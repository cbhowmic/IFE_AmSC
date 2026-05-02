from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn


class SurrogateMLP(nn.Module):
    def __init__(self, input_dim: int, output_dim: int, hidden_dim: int = 100, num_hidden_layers: int = 3):
        super().__init__()

        layers = []
        in_dim = input_dim

        for _ in range(num_hidden_layers):
            layers.append(nn.Linear(in_dim, hidden_dim))
            layers.append(nn.ReLU())
            in_dim = hidden_dim

        layers.append(nn.Linear(hidden_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class RhinoSurrogate:
    def __init__(self, artifact_path: str, device: str = "cpu"):
        self.device = torch.device(device)

        artifact = torch.load(artifact_path, map_location=self.device)

        self.model_config = artifact["model_config"]
        self.input_specs = artifact["input_specs"]
        self.output_specs = artifact["output_specs"]

        self.model = SurrogateMLP(
            input_dim=self.model_config["input_dim"],
            output_dim=self.model_config["output_dim"],
            hidden_dim=self.model_config["hidden_dim"],
            num_hidden_layers=self.model_config["num_hidden_layers"],
        ).to(self.device)

        self.model.load_state_dict(artifact["model_state_dict"])
        self.model.eval()

        self.x_mean = artifact["x_mean"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)
        self.x_std  = artifact["x_std"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)
        self.y_mean = artifact["y_mean"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)
        self.y_std  = artifact["y_std"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)

        self.input_keys = [spec["key"] for spec in self.input_specs]
        self.input_display_names = [spec["display_name"] for spec in self.input_specs]
        self.output_keys = [spec["key"] for spec in self.output_specs]
        self.output_display_names = [spec["display_name"] for spec in self.output_specs]

    def predict_array(self, x_raw: np.ndarray) -> np.ndarray:
        x_raw = np.asarray(x_raw, dtype=np.float32).reshape(1, -1)

        if x_raw.shape[1] != len(self.input_specs):
            raise ValueError(
                f"Expected {len(self.input_specs)} inputs, got {x_raw.shape[1]}"
            )

        x_norm = (x_raw - self.x_mean) / self.x_std
        x_tensor = torch.tensor(x_norm, dtype=torch.float32, device=self.device)

        with torch.no_grad():
            y_norm = self.model(x_tensor).detach().cpu().numpy()

        y_raw = y_norm * self.y_std + self.y_mean
        return y_raw.reshape(-1)

    def predict_from_dict(self, inputs: dict[str, float]) -> dict[str, float]:
        missing = [k for k in self.input_keys if k not in inputs]
        if missing:
            raise ValueError(f"Missing input keys: {missing}")

        x_raw = np.array([[inputs[k] for k in self.input_keys]], dtype=np.float32)
        y_raw = self.predict_array(x_raw)

        return {
            spec["key"]: float(value)
            for spec, value in zip(self.output_specs, y_raw)
        }

    def predict_from_named_args(
        self,
        I0_SD: float,
        Ndotminus: float,
        beta: float,
    ) -> dict[str, float]:
        return self.predict_from_dict(
            {
                "I0_SD": I0_SD,
                "Ndotminus": Ndotminus,
                "beta": beta,
            }
        )

    def predict_with_display_names(self, inputs: dict[str, float]) -> dict[str, float]:
        y_by_key = self.predict_from_dict(inputs)

        return {
            spec["display_name"]: y_by_key[spec["key"]]
            for spec in self.output_specs
        }