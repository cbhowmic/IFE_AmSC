#!/usr/bin/env python3
"""Build React Flow graph payloads from RHINO/openPMD simulations.

This module:
1. loads RHINO/openPMD data for one simulation
2. extracts subsystem names, injectors, masses, and fractional inflows
3. downsamples the mass and time series
4. builds a React Flow-compatible payload
5. optionally injects surrogate prediction results for frontend display
6. treats missing fractional inflows as intentional wall-loss behavior

The returned payload is intended to be consumed by the React frontend.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
import openpmd_api as io


N_COMPONENTS = 21
DOWNSAMPLE_STEP = 1
SPECIES = "Tritium"
LAYOUT_SCALE_X = 0.86
LAYOUT_SCALE_Y = 0.86

# Keep your SAVED_POSITIONS block exactly as you already have it.
SAVED_POSITIONS = [
    {
      "id": "0",
      "position": {"x": 1946, "y": -133}
    },
    {
      "id": "1",
      "position": {"x": 71.64952667055849, "y": -133}
    },
    {
      "id": "2",
      "position": {"x": 558.3820898538551, "y": -133}
    },
    {
      "id": "3",
      "position": {"x": 1360, "y": -133}
    },
    {
      "id": "4",
      "position": {"x": 1656.8719706654426, "y": -133}
    },
    {
      "id": "5",
      "position": {"x": 1360, "y": 133}
    },
    {
      "id": "6",
      "position": {"x": 1360, "y": 762}
    },
    {
      "id": "7",
      "position": {"x": 1946, "y": 762}
    },
    {
      "id": "8",
      "position": {"x": 1360, "y": 524.2087731367648}
    },
    {
      "id": "9",
      "position": {"x": 2288, "y": 372}
    },
    {
      "id": "10",
      "position": {"x": 1946, "y": 372}
    },
    {
      "id": "11",
      "position": {"x": 1012, "y": 372}
    },
    {
      "id": "12",
      "position": {"x": 1012, "y": 762}
    },
    {
      "id": "13",
      "position": {"x": 1012, "y": 133}
    },
    {
      "id": "14",
      "position": {"x": 1360, "y": 1018}
    },
    {
      "id": "15",
      "position": {"x": 512, "y": 762}
    },
    {
      "id": "16",
      "position": {"x": 54.363138197245405, "y": 383.64780531095334}
    },
    {
      "id": "17",
      "position": {"x": 2572.9101085153334, "y": 372}
    },
    {
      "id": "18",
      "position": {"x": 559.1109837515487, "y": -394.34766136899736}
    },
    {
      "id": "19",
      "position": {"x": 1946, "y": 1018}
    },
    {
      "id": "20",
      "position": {"x": 693.8872040292795, "y": 1035.8068817822084}
    },
]

WALL_LOSS_NODE = {
    "id": "wall_loss",
    "position": {"x": 1011, "y": 1635},
    "data": {
        "label": "Wall Loss",
        "initialMass": 0.0,
        "massSeries": [],
        "isLossNode": True,
        "isDiagnostic": False,
        "hasMeaningfulSteadyState": False,
    },
}

SURROGATE_PREDICTION_DEFAULTS = {
    "isotopeSeparationTritium": {
        "label": "Isotope separation T",
        "value": None,
        "unit": "g",
    },
    "plantDoublingTime": {
        "label": "Plant doubling time",
        "value": None,
        "unit": "d",
    },
    "minimalStartupInventory": {
        "label": "Minimal startup inventory",
        "value": None,
        "unit": "g",
    },
}


def is_diagnostic_subsystem(name: str) -> bool:
    return "box" in name.lower()


def open_series(bp5_path: Path) -> io.Series:
    return io.Series(
        str(bp5_path),
        io.Access.read_only,
        '{"verify_homogeneous_extents": false}',
    )


def _get_first_snapshot(series: io.Series):
    snapshots = series.snapshots()
    if not snapshots:
        raise ValueError("The series does not contain any snapshots.")
    return snapshots[0]


def _validate_species(snapshot: Any, species: str) -> None:
    available = [sp for sp in snapshot.particles if sp != "Times"]
    if species == "Times" or species not in available:
        raise ValueError(
            f"Invalid species {species!r}. Available species are: {available}"
        )


def get_mass_inventory(series: io.Series, species: str = SPECIES) -> np.ndarray:
    snapshot = _get_first_snapshot(series)
    _validate_species(snapshot, species)
    mass = snapshot.particles[species]["mass"][io.Record_Component.SCALAR]
    data = mass.load_chunk()
    series.flush()
    return data


def get_time_series(series: io.Series) -> np.ndarray:
    snapshot = _get_first_snapshot(series)
    times = snapshot.particles["Times"]["data"][io.Record_Component.SCALAR]
    data = times.load_chunk()
    series.flush()
    return data


def get_subsystems_dict(series: io.Series, species: str = SPECIES) -> dict[str, int]:
    snapshot = _get_first_snapshot(series)
    _validate_species(snapshot, species)
    return {
        name: subsystem.get_attribute("id")
        for name, subsystem in snapshot.particles[species]["subsystems"].items()
    }


def get_injectors_map(series: io.Series, species: str = SPECIES) -> dict[int, Any]:
    snapshot = _get_first_snapshot(series)
    _validate_species(snapshot, species)
    return {
        subsystem.get_attribute("id"): subsystem.get_attribute("injectors")
        for _, subsystem in snapshot.particles[species]["subsystems"].items()
    }


def get_flows_map(series: io.Series, species: str = SPECIES) -> dict[int, Any]:
    snapshot = _get_first_snapshot(series)
    _validate_species(snapshot, species)
    return {
        subsystem.get_attribute("id"): subsystem.get_attribute("fractional inflows")
        for _, subsystem in snapshot.particles[species]["subsystems"].items()
    }


def normalize_to_list(value: Any) -> list[Any]:
    if value is None or value == "None":
        return []
    if isinstance(value, list):
        return value
    return [value]


def compact_positions(
    positions: list[dict[str, Any]],
    scale_x: float = LAYOUT_SCALE_X,
    scale_y: float = LAYOUT_SCALE_Y,
) -> list[dict[str, Any]]:
    xs = [float(item["position"]["x"]) for item in positions]
    ys = [float(item["position"]["y"]) for item in positions]
    center_x = 0.5 * (min(xs) + max(xs))
    center_y = 0.5 * (min(ys) + max(ys))

    compacted: list[dict[str, Any]] = []
    for item in positions:
        x = float(item["position"]["x"])
        y = float(item["position"]["y"])
        compacted.append(
            {
                "id": item["id"],
                "position": {
                    "x": center_x + (x - center_x) * scale_x,
                    "y": center_y + (y - center_y) * scale_y,
                },
            }
        )

    return compacted


def compact_single_position(
    position: dict[str, float],
    reference_positions: list[dict[str, Any]],
    scale_x: float = LAYOUT_SCALE_X,
    scale_y: float = LAYOUT_SCALE_Y,
) -> dict[str, float]:
    xs = [float(item["position"]["x"]) for item in reference_positions]
    ys = [float(item["position"]["y"]) for item in reference_positions]
    center_x = 0.5 * (min(xs) + max(xs))
    center_y = 0.5 * (min(ys) + max(ys))

    x = float(position["x"])
    y = float(position["y"])
    return {
        "x": center_x + (x - center_x) * scale_x,
        "y": center_y + (y - center_y) * scale_y,
    }


def build_react_prediction_payload(
    prediction_display: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Convert surrogate display-name outputs into the React payload schema."""
    payload = copy.deepcopy(SURROGATE_PREDICTION_DEFAULTS)

    if prediction_display is None:
        return payload

    payload["isotopeSeparationTritium"]["value"] = prediction_display.get(
        "Tritium in isotope separation [g]"
    )
    payload["plantDoublingTime"]["value"] = prediction_display.get(
        "Plant doubling time [d]"
    )
    payload["minimalStartupInventory"]["value"] = prediction_display.get(
        "Minimum startup inventory [g]"
    )

    return payload


def with_payload_defaults(payload: dict[str, Any]) -> dict[str, Any]:
    """Backfill optional fields so older payloads remain frontend-compatible."""
    for node in payload.get("nodes", []):
        data = node.setdefault("data", {})
        label = str(data.get("label", ""))
        is_diagnostic = bool(data.get("isDiagnostic", is_diagnostic_subsystem(label)))
        data.setdefault("isDiagnostic", is_diagnostic)
        data.setdefault("hasMeaningfulSteadyState", not is_diagnostic)

    payload.setdefault("time", {})
    payload["time"].setdefault("timeSeries", [])
    payload["time"].setdefault("timeUnit", "days")
    payload["time"].setdefault("downsampleStep", None)

    payload.setdefault(
        "surrogatePredictions",
        copy.deepcopy(SURROGATE_PREDICTION_DEFAULTS),
    )
    payload.setdefault("simulation_path", None)

    return payload


def build_nodes_rf(
    names_map: dict[int, str],
    mass_ds: np.ndarray,
    saved_positions: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    compacted_positions = compact_positions(saved_positions)
    pos_map = {int(item["id"]): item["position"] for item in compacted_positions}

    nodes_rf: list[dict[str, Any]] = []
    for component_id in range(mass_ds.shape[0]):
        if component_id not in pos_map:
            raise KeyError(f"Missing saved position for component {component_id}")

        label = names_map.get(component_id, f"Subsystem {component_id}")
        is_diagnostic = is_diagnostic_subsystem(label)
        p = pos_map[component_id]

        nodes_rf.append(
            {
                "id": str(component_id),
                "position": {
                    "x": float(p["x"]),
                    "y": float(p["y"]),
                },
                "data": {
                    "label": label,
                    "initialMass": float(mass_ds[component_id, 0]),
                    "massSeries": mass_ds[component_id, :].tolist(),
                    "isDiagnostic": is_diagnostic,
                    "hasMeaningfulSteadyState": not is_diagnostic,
                },
            }
        )

    wall_loss_node = copy.deepcopy(WALL_LOSS_NODE)
    wall_loss_node["position"] = compact_single_position(
        wall_loss_node["position"],
        saved_positions,
    )
    wall_loss_node["data"]["massSeries"] = [0.0] * mass_ds.shape[1]
    nodes_rf.append(wall_loss_node)

    return nodes_rf


def build_edges_rf(
    injectors_map: dict[int, Any],
    flows_map: dict[int, Any],
) -> list[dict[str, Any]]:
    edges_rf: list[dict[str, Any]] = []

    for target, raw_injectors in injectors_map.items():
        injectors = normalize_to_list(raw_injectors)
        if not injectors:
            continue

        raw_flows = flows_map.get(target, None)
        flow_missing = raw_flows is None or raw_flows == "None"

        if flow_missing:
            for source in injectors:
                edges_rf.append(
                    {
                        "id": f"{source}-{target}",
                        "source": str(source),
                        "target": str(target),
                        "data": {
                            "isLoss": False,
                            "isUnknownFlow": True,
                            "baseFlow": None,
                        },
                    }
                )

            edges_rf.append(
                {
                    "id": f"{target}-wall_loss",
                    "source": str(target),
                    "target": "wall_loss",
                    "data": {
                        "isLoss": True,
                        "isUnknownFlow": False,
                        "baseFlow": None,
                    },
                }
            )
            continue

        flows = normalize_to_list(raw_flows)
        if len(injectors) != len(flows):
            print(
                f"Skipping component {target}: mismatch "
                f"injectors={len(injectors)} flows={len(flows)}"
            )
            continue

        for source, flow_value in zip(injectors, flows):
            edges_rf.append(
                {
                    "id": f"{source}-{target}",
                    "source": str(source),
                    "target": str(target),
                    "data": {
                        "isLoss": False,
                        "isUnknownFlow": False,
                        "baseFlow": float(flow_value),
                    },
                }
            )

    return edges_rf


def build_graph_from_simulation(
    sim_path: str | Path,
    downsample_step: int = DOWNSAMPLE_STEP,
    surrogate_predictions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the React payload for one archived RHINO simulation."""
    if downsample_step < 1:
        raise ValueError("downsample_step must be >= 1")

    sim_path = Path(sim_path)
    series = open_series(sim_path)

    try:
        mass = get_mass_inventory(series)
        time = get_time_series(series)

        if mass.shape[0] != N_COMPONENTS:
            print(
                f"Warning: expected {N_COMPONENTS} components, found {mass.shape[0]}"
            )

        mass_ds = mass[:, ::downsample_step]
        time_ds = time[::downsample_step]

        subsystems = get_subsystems_dict(series)
        names_map = {value: key for key, value in subsystems.items()}
        injectors_map = get_injectors_map(series)
        flows_map = get_flows_map(series)

        nodes_rf = build_nodes_rf(names_map, mass_ds, SAVED_POSITIONS)
        edges_rf = build_edges_rf(injectors_map, flows_map)

        payload = {
            "simulation_path": str(sim_path),
            "nodes": nodes_rf,
            "edges": edges_rf,
            "time": {
                "timeSeries": time_ds.tolist(),
                "timeUnit": "days",
                "downsampleStep": downsample_step,
            },
            "surrogatePredictions": (
                copy.deepcopy(surrogate_predictions)
                if surrogate_predictions is not None
                else copy.deepcopy(SURROGATE_PREDICTION_DEFAULTS)
            ),
        }

        return with_payload_defaults(payload)

    finally:
        series.close()


def write_graph_payload(
    payload: dict[str, Any],
    payload_out: str | None = None,
    nodes_out: str = "nodes_rf.json",
    edges_out: str = "edges_rf.json",
    time_out: str = "time_rf.json",
) -> None:
    """Write either a full payload JSON or split files for frontend debugging."""
    payload = with_payload_defaults(payload)

    if payload_out is not None:
        Path(payload_out).write_text(json.dumps(payload))

    Path(nodes_out).write_text(json.dumps(payload["nodes"]))
    Path(edges_out).write_text(json.dumps(payload["edges"]))
    Path(time_out).write_text(
        json.dumps(
            {
                "time": payload["time"],
                "surrogatePredictions": payload["surrogatePredictions"],
                "simulation_path": payload["simulation_path"],
            }
        )
    )