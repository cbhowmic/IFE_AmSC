from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import torch


DEFAULT_DATA_ROOT = Path(__file__).resolve().parents[1] / "data" / "surrogate_bp_output"


class SimulationLookup:
    def __init__(
        self,
        index_path: str | Path,
        artifact_path: str | Path,
        data_root: str | Path | None = None,
    ):
        self.index_path = Path(index_path).expanduser().resolve()
        env_data_root = os.environ.get("RHINO_SIM_DATA_ROOT")
        configured_data_root = data_root if data_root is not None else env_data_root
        self.data_root = (
            Path(configured_data_root).expanduser()
            if configured_data_root
            else DEFAULT_DATA_ROOT
        )
        if not self.data_root.is_absolute():
            self.data_root = (self.index_path.parent / self.data_root).resolve()
        else:
            self.data_root = self.data_root.resolve()

        with self.index_path.open() as f:
            payload = json.load(f)

        self.input_specs = payload["input_specs"]
        self.output_specs = payload["output_specs"]
        self.records = payload["records"]

        artifact = torch.load(artifact_path, map_location="cpu")
        self.x_mean = artifact["x_mean"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)
        self.x_std = artifact["x_std"].detach().cpu().numpy().reshape(1, -1).astype(np.float32)

        self.input_keys = [spec["key"] for spec in self.input_specs]

        # Build one matrix of all indexed inputs in the same order as input_specs
        self.X = np.array(
            [
                [record["inputs"][key] for key in self.input_keys]
                for record in self.records
            ],
            dtype=np.float32,
        )

        self.X_norm = (self.X - self.x_mean) / self.x_std

    def resolve_simulation_path(self, indexed_path: str | Path) -> Path:
        """Resolve archived index paths against the local portable data folder."""
        path = Path(indexed_path).expanduser()

        if path.exists():
            return path.resolve()

        parts = path.parts
        if "surrogate_bp_output" in parts:
            data_root_index = parts.index("surrogate_bp_output")
            suffix = parts[data_root_index + 1 :]
            if suffix:
                return (self.data_root.joinpath(*suffix)).resolve()

        if path.is_absolute():
            return path

        data_root_candidate = self.data_root / path
        if data_root_candidate.exists():
            return data_root_candidate.resolve()

        return (self.index_path.parent / path).resolve()

    def _query_to_array(self, query_inputs: dict[str, float]) -> np.ndarray:
        missing = [k for k in self.input_keys if k not in query_inputs]
        if missing:
            raise ValueError(f"Missing input keys: {missing}")

        return np.array(
            [[query_inputs[key] for key in self.input_keys]],
            dtype=np.float32,
        )

    def find_nearest(self, query_inputs: dict[str, float]) -> dict:
        x_query = self._query_to_array(query_inputs)
        x_query_norm = (x_query - self.x_mean) / self.x_std

        distances = np.linalg.norm(self.X_norm - x_query_norm, axis=1)
        idx = int(np.argmin(distances))

        nearest_record = self.records[idx]
        resolved_path = self.resolve_simulation_path(nearest_record["path"])

        return {
            "query_inputs": {k: float(query_inputs[k]) for k in self.input_keys},
            "nearest_index": idx,
            "distance_normalized": float(distances[idx]),
            "simulation_id": nearest_record["simulation_id"],
            "path": str(resolved_path),
            "indexed_path": nearest_record["path"],
            "path_exists": resolved_path.exists(),
            "scenario": nearest_record["scenario"],
            "inputs": nearest_record["inputs"],
            "outputs": nearest_record["outputs"],
        }

    def find_k_nearest(self, query_inputs: dict[str, float], k: int = 2) -> list[dict]:
        x_query = self._query_to_array(query_inputs)
        x_query_norm = (x_query - self.x_mean) / self.x_std

        distances = np.linalg.norm(self.X_norm - x_query_norm, axis=1)
        order = np.argsort(distances)[:k]

        results = []
        for idx in order:
            record = self.records[int(idx)]
            resolved_path = self.resolve_simulation_path(record["path"])
            results.append(
                {
                    "nearest_index": int(idx),
                    "distance_normalized": float(distances[idx]),
                    "simulation_id": record["simulation_id"],
                    "path": str(resolved_path),
                    "indexed_path": record["path"],
                    "path_exists": resolved_path.exists(),
                    "scenario": record["scenario"],
                    "inputs": record["inputs"],
                    "outputs": record["outputs"],
                }
            )

        return results
