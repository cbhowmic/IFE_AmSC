"""MLflow serving wrapper for a normalized RHINO PyTorch surrogate."""

from __future__ import annotations

from typing import Any, Sequence

import mlflow.pyfunc
import numpy as np
import pandas as pd
import torch


class RhinoSurrogatePyFunc(mlflow.pyfunc.PythonModel):
    """Serve raw named RHINO inputs and return denormalized named outputs."""

    def __init__(
        self,
        model: torch.nn.Module,
        input_columns: Sequence[str],
        output_columns: Sequence[str],
        x_mean: np.ndarray,
        x_std: np.ndarray,
        y_mean: np.ndarray,
        y_std: np.ndarray,
    ) -> None:
        """Store the trained network and its complete inference contract."""
        self.model = model.to("cpu").eval()
        self.input_columns = list(input_columns)
        self.output_columns = list(output_columns)
        self.x_mean = np.asarray(x_mean, dtype=np.float32).reshape(-1)
        self.x_std = np.asarray(x_std, dtype=np.float32).reshape(-1)
        self.y_mean = np.asarray(y_mean, dtype=np.float32).reshape(-1)
        self.y_std = np.asarray(y_std, dtype=np.float32).reshape(-1)

    def predict(
        self,
        context: Any,
        model_input: pd.DataFrame,
        params: dict[str, Any] | None = None,
    ) -> pd.DataFrame:
        """Normalize raw inputs, run inference, and restore physical outputs."""
        del context, params
        if not isinstance(model_input, pd.DataFrame):
            model_input = pd.DataFrame(model_input)

        missing = [
            column for column in self.input_columns if column not in model_input
        ]
        if missing:
            raise ValueError(
                "Inference input is missing required column(s): "
                + ", ".join(missing)
            )

        try:
            values = model_input[self.input_columns].apply(
                pd.to_numeric,
                errors="raise",
            )
        except (TypeError, ValueError) as error:
            raise ValueError("Inference input columns must be numeric") from error

        x = values.to_numpy(dtype=np.float32)
        if not np.isfinite(x).all():
            raise ValueError("Inference inputs must not contain NaN or infinity")
        normalized_x = (x - self.x_mean) / self.x_std

        with torch.no_grad():
            normalized_y = self.model(
                torch.tensor(normalized_x, dtype=torch.float32)
            ).numpy()
        y = normalized_y * self.y_std + self.y_mean
        return pd.DataFrame(y, columns=self.output_columns, index=model_input.index)
