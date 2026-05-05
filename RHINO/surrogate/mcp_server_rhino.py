from __future__ import annotations

import os
import socket
import subprocess
import time
from pathlib import Path
from urllib.parse import quote

from fastmcp import FastMCP
from starlette.responses import JSONResponse
from rhino_surrogate_runtime import RhinoSurrogate
from simulation_lookup import SimulationLookup 
from graph_builder import (
    build_graph_from_simulation,
    build_react_prediction_payload,
    write_graph_payload,
)

mcp = FastMCP("RHINO Surrogate")

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
MODEL_PATH = DATA_DIR / "rhino_surrogate.pt"
INDEX_PATH = DATA_DIR / "simulation_index.json"
VIZ_DIR = BASE_DIR / "react_flow_viz"
VIZ_PUBLIC_DIR = VIZ_DIR / "public"

_visualization_process: subprocess.Popen | None = None
_latest_graph: dict | None = None
_graph_revision = 0


# Load once when the server starts
surrogate = RhinoSurrogate(str(MODEL_PATH), device="cpu")

lookup = SimulationLookup(
    index_path=INDEX_PATH,
    artifact_path=MODEL_PATH,
)


def _is_port_open(host: str, port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.25)
        return sock.connect_ex((host, port)) == 0


def _write_visualization_files(graph: dict) -> None:
    VIZ_PUBLIC_DIR.mkdir(parents=True, exist_ok=True)
    write_graph_payload(
        graph,
        payload_out=str(VIZ_PUBLIC_DIR / "graph_payload.json"),
        nodes_out=str(VIZ_PUBLIC_DIR / "nodes_rf.json"),
        edges_out=str(VIZ_PUBLIC_DIR / "edges_rf.json"),
        time_out=str(VIZ_PUBLIC_DIR / "time_rf.json"),
    )


def _json_response(payload: dict | list) -> JSONResponse:
    return JSONResponse(
        payload,
        headers={
            "Access-Control-Allow-Origin": "*",
            "Cache-Control": "no-store",
        },
    )


def _mcp_base_url() -> str:
    host = os.environ.get("RHINO_MCP_HOST", "127.0.0.1")
    port = int(os.environ.get("RHINO_MCP_PORT", "8000"))
    return f"http://{host}:{port}"


def _read_fallback_json(file_name: str) -> dict | list:
    fallback_path = VIZ_PUBLIC_DIR / file_name
    if fallback_path.exists():
        import json

        return json.loads(fallback_path.read_text())

    raise FileNotFoundError(f"No fallback JSON found at {fallback_path}")


def _current_graph_payload() -> dict:
    if _latest_graph is not None:
        return {**_latest_graph, "graphRevision": _graph_revision}

    fallback = _read_fallback_json("graph_payload.json")
    if isinstance(fallback, dict):
        return {**fallback, "graphRevision": 0}

    raise ValueError("Fallback graph_payload.json is not a JSON object")


def _start_visualization_server(host: str = "127.0.0.1", port: int = 5173) -> dict:
    global _visualization_process

    data_base_url = _mcp_base_url()
    url = f"http://{host}:{port}/?dataBaseUrl={quote(data_base_url, safe='')}"
    if _is_port_open(host, port):
        return {"status": "already_running", "url": url, "port": port}

    if not (VIZ_DIR / "package.json").exists():
        return {
            "status": "not_started",
            "reason": f"Missing visualization app at {VIZ_DIR}",
            "url": None,
            "port": port,
        }

    if not (VIZ_DIR / "node_modules").exists():
        return {
            "status": "not_started",
            "reason": f"Missing node_modules in {VIZ_DIR}; run `npm install` there before the demo.",
            "url": None,
            "port": port,
        }

    _visualization_process = subprocess.Popen(
        ["npm", "run", "dev", "--", "--host", host, "--port", str(port)],
        cwd=VIZ_DIR,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    time.sleep(1.0)

    if _visualization_process.poll() is not None:
        return {
            "status": "failed",
            "reason": "Vite exited immediately.",
            "url": None,
            "port": port,
        }

    return {
        "status": "started",
        "url": url,
        "port": port,
        "pid": _visualization_process.pid,
    }


def _add_nearest_simulation_context(
    graph: dict,
    nearest: dict,
    prediction_by_key: dict[str, float],
    prediction_payload: dict,
) -> dict:
    for spec in lookup.output_specs:
        key = spec["key"]
        display_name = spec["display_name"]
        pred_val = float(prediction_by_key[key])
        sim_val = float(nearest["outputs"][key])
        abs_err = abs(pred_val - sim_val)
        rel_err = abs_err / max(abs(sim_val), 1e-12)

        if display_name == "Tritium in isotope separation [g]":
            payload_key = "isotopeSeparationTritium"
        elif display_name == "Plant doubling time [d]":
            payload_key = "plantDoublingTime"
        elif display_name == "Minimum startup inventory [g]":
            payload_key = "minimalStartupInventory"
        else:
            continue

        prediction_payload[payload_key].update(
            {
                "nearestSimulation": sim_val,
                "absoluteError": abs_err,
                "relativeError": rel_err,
            }
        )

    simulation_path = Path(nearest["path"])
    graph["surrogatePredictions"] = prediction_payload
    graph["simulationSummary"] = {
        "simulationId": nearest["simulation_id"],
        "runFolder": simulation_path.name,
        "scenario": nearest["scenario"],
        "path": nearest["path"],
        "indexedPath": nearest.get("indexed_path"),
        "pathExists": nearest.get("path_exists"),
        "distanceNormalized": nearest["distance_normalized"],
        "inputs": nearest["inputs"],
        "outputs": nearest["outputs"],
    }

    return graph
    
@mcp.tool
def predict_rhino_surrogate(
    Ndotminus: float,
    beta: float,
) -> dict[str, float]:
    """
    Predict RHINO surrogate outputs from 2 physical input values.

    Inputs:
    - Ndotminus: tritium burning rate in grams per day [g/d]
    - beta: burn fraction as a unitless fraction [-]

    The inputs must be physical values, not normalized values.
    Normalization is handled internally by the surrogate runtime.

    Returns:
    A dictionary mapping human-readable output names to predicted physical values.
    """
    return surrogate.predict_with_display_names(
        {
            "Ndotminus": Ndotminus,
            "beta": beta,
        }
    )

@mcp.tool
def find_nearest_simulation(
    Ndotminus: float,
    beta: float,
) -> dict:
    """
    Find the nearest archived simulation in normalized input space.

    Inputs are physical values:
    - Ndotminus in grams per day [g/d]
    - beta as a unitless fraction [-]

    The lookup uses internally normalized inputs to compute distance, but callers
    should provide physical values.

    Returns:
    Metadata for the nearest archived simulation, including its path, stored inputs,
    stored outputs, and the normalized input-space distance from the query point.
    """
    return lookup.find_nearest(
        {
            "Ndotminus": Ndotminus,
            "beta": beta,
        }
    )

@mcp.tool
def predict_and_compare_to_nearest_simulation(
    Ndotminus: float,
    beta: float,
) -> dict:
    """
    Predict RHINO surrogate outputs from physical inputs, find the nearest archived
    simulation in normalized input space, and compare the prediction to that
    simulation's outputs.
    
    Inputs:
    - Ndotminus: tritium burning rate [g/d]
    - beta: burn fraction [-]
    
    Inputs must be given in physical units, not normalized units. Normalization is
    handled internally by the surrogate runtime.
    
    Returns:
    A dictionary with the following entries:
    
    - prediction:
        Dictionary mapping human-readable output names to surrogate-predicted
        physical values. Current outputs are:
        * Tritium in isotope separation [g]
        * Plant doubling time [d]
        * Minimum startup inventory [g]
    
    - nearest_simulation:
        Dictionary describing the nearest archived simulation in normalized input
        space. It contains:
        * simulation_id: identifier of the archived simulation
        * path: filesystem path to the simulation data
        * scenario: scenario or campaign label
        * distance_normalized: Euclidean distance from the query point in normalized
          input space
        * inputs: dictionary of the archived simulation input values
        * outputs: dictionary of the archived simulation output values
    
    - comparison:
        Dictionary keyed by human-readable output name. For each output, it stores:
        * predicted: surrogate-predicted value
        * nearest_simulation: value from the nearest archived simulation
        * absolute_error: absolute difference between prediction and simulation
        * relative_error: relative difference between prediction and simulation
    
    Notes:
    The nearest simulation is selected by distance in normalized input space, not by
    matching outputs. If the nearest-neighbor distance is large, the comparison may
    be less representative.

    This tool requires access to a local simulation index and archived simulation metadata.
    """
    query_inputs = {
        "Ndotminus": Ndotminus,
        "beta": beta,
    }

    prediction = surrogate.predict_from_dict(query_inputs)
    
    prediction_display = surrogate.predict_with_display_names(query_inputs)
    nearest = lookup.find_nearest(query_inputs)

    comparison = {}
    for spec in lookup.output_specs:
        key = spec["key"]
        name = spec["display_name"]

        pred_val = float(prediction[key])
        sim_val = float(nearest["outputs"][key])

        abs_err = abs(pred_val - sim_val)
        rel_err = abs_err / max(abs(sim_val), 1e-12)

        comparison[name] = {
            "predicted": pred_val,
            "nearest_simulation": sim_val,
            "absolute_error": abs_err,
            "relative_error": rel_err,
        }

    return {
        "prediction": prediction_display,
        "nearest_simulation": nearest,
        "comparison": comparison,
    }


@mcp.tool
def build_graph_for_nearest_simulation(
    Ndotminus: float,
    beta: float,
    write_visualization_files: bool = False,
    start_visualization: bool = True,
    visualization_port: int = 5173,
) -> dict:
    """
    Find the nearest archived simulation for the given physical inputs and
    return React Flow graph data for that simulation.

    By default this stores the latest graph in memory, exposes it through MCP
    HTTP JSON routes, and starts the local Vite visualization server. Existing
    JSON files remain available as fallback data.
    """
    global _latest_graph
    global _graph_revision

    query_inputs = {
        "Ndotminus": Ndotminus,
        "beta": beta,
    }

    nearest = lookup.find_nearest(query_inputs)
    prediction_by_key = surrogate.predict_from_dict(query_inputs)
    prediction_display = surrogate.predict_with_display_names(query_inputs)
    react_predictions = build_react_prediction_payload(prediction_display)

    graph = build_graph_from_simulation(
        nearest["path"],
        surrogate_predictions=react_predictions,
    )
    graph = _add_nearest_simulation_context(
        graph,
        nearest,
        prediction_by_key,
        react_predictions,
    )
    _latest_graph = graph
    _graph_revision += 1

    visualization = {
        "status": "not_requested",
        "url": None,
        "port": visualization_port,
        "data_source": f"{_mcp_base_url()}/graph_payload.json",
    }
    if write_visualization_files:
        _write_visualization_files(graph)
        visualization["files_written"] = {
            "payload": str(VIZ_PUBLIC_DIR / "graph_payload.json"),
            "nodes": str(VIZ_PUBLIC_DIR / "nodes_rf.json"),
            "edges": str(VIZ_PUBLIC_DIR / "edges_rf.json"),
            "time": str(VIZ_PUBLIC_DIR / "time_rf.json"),
        }

    if start_visualization:
        visualization.update(
            _start_visualization_server(port=int(visualization_port))
        )

    return {
        "nearest_simulation": nearest,
        "graph": graph,
        "visualization": visualization,
    }

@mcp.custom_route("/health", methods=["GET"])
async def health_check(request):
    return JSONResponse({"status": "healthy", "service": "rhino-mcp"})


@mcp.custom_route("/graph_payload.json", methods=["GET"])
async def graph_payload(request):
    return _json_response(_current_graph_payload())


@mcp.custom_route("/graph_version.json", methods=["GET"])
async def graph_version(request):
    return _json_response(
        {
            "graphRevision": _graph_revision,
            "hasLiveGraph": _latest_graph is not None,
        }
    )


@mcp.custom_route("/nodes_rf.json", methods=["GET"])
async def nodes_payload(request):
    if _latest_graph is not None:
        return _json_response(_latest_graph["nodes"])
    return _json_response(_read_fallback_json("nodes_rf.json"))


@mcp.custom_route("/edges_rf.json", methods=["GET"])
async def edges_payload(request):
    if _latest_graph is not None:
        return _json_response(_latest_graph["edges"])
    return _json_response(_read_fallback_json("edges_rf.json"))


@mcp.custom_route("/time_rf.json", methods=["GET"])
async def time_payload(request):
    if _latest_graph is not None:
        return _json_response(
            {
                "time": _latest_graph["time"],
                "surrogatePredictions": _latest_graph["surrogatePredictions"],
                "simulation_path": _latest_graph["simulation_path"],
            }
        )
    return _json_response(_read_fallback_json("time_rf.json"))

if __name__ == "__main__":
    mcp.run(
        transport="http",
        host=os.environ.get("RHINO_MCP_HOST", "127.0.0.1"),
        port=int(os.environ.get("RHINO_MCP_PORT", "8000")),
    )
