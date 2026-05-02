from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch


class SimulationLookup:
    def __init__(self, index_path: str, artifact_path: str):
        with open(index_path) as f:
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

        return {
            "query_inputs": {k: float(query_inputs[k]) for k in self.input_keys},
            "nearest_index": idx,
            "distance_normalized": float(distances[idx]),
            "simulation_id": nearest_record["simulation_id"],
            "path": nearest_record["path"],
            "scenario": nearest_record["scenario"],
            "inputs": nearest_record["inputs"],
            "outputs": nearest_record["outputs"],
        }

    def find_k_nearest(self, query_inputs: dict[str, float], k: int = 3) -> list[dict]:
        x_query = self._query_to_array(query_inputs)
        x_query_norm = (x_query - self.x_mean) / self.x_std

        distances = np.linalg.norm(self.X_norm - x_query_norm, axis=1)
        order = np.argsort(distances)[:k]

        results = []
        for idx in order:
            record = self.records[int(idx)]
            results.append(
                {
                    "nearest_index": int(idx),
                    "distance_normalized": float(distances[idx]),
                    "simulation_id": record["simulation_id"],
                    "path": record["path"],
                    "scenario": record["scenario"],
                    "inputs": record["inputs"],
                    "outputs": record["outputs"],
                }
            )

        return results