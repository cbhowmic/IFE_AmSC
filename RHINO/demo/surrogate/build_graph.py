from pathlib import Path

from RHINO.demo.surrogate.rhino_surrogate_runtime import RhinoSurrogate
from RHINO.demo.surrogate.graph_builder import (
    build_graph_from_simulation,
    write_graph_payload,
    build_react_prediction_payload,
)
from RHINO.demo.surrogate.simulation_lookup import SimulationLookup


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
VIZ_PUBLIC_DIR = BASE_DIR / "react_flow_viz" / "public"


def main():
    # Physical inputs
    inputs = {
        "Ndotminus": 500.0,
        "beta": 0.1,
    }

    # Load surrogate
    surrogate = RhinoSurrogate(DATA_DIR / "rhino_surrogate.pt", device="cpu")

    # Predict surrogate outputs
    prediction_display = surrogate.predict_with_display_names(inputs)
    print("\nSurrogate prediction:")
    print(prediction_display)

    # Load nearest-simulation lookup
    lookup = SimulationLookup(
        index_path=DATA_DIR / "simulation_index.json",
        artifact_path=DATA_DIR / "rhino_surrogate.pt",
    )

    # Find nearest archived simulation
    nearest = lookup.find_nearest(inputs)
    print("\nNearest simulation:")
    print("simulation_id:", nearest["simulation_id"])
    print("distance_normalized:", nearest["distance_normalized"])
    print("inputs:", nearest["inputs"])
    print("outputs:", nearest["outputs"])

    # Convert surrogate outputs to React-friendly payload
    prediction_payload = build_react_prediction_payload(prediction_display)

    # Build graph from nearest simulation
    sim_path = nearest["path"]
    graph = build_graph_from_simulation(
        sim_path,
        surrogate_predictions=prediction_payload,
    )

    print("\nGraph payload:")
    print(graph.keys())
    print(graph["time"].keys())
    print(graph["surrogatePredictions"])

    # Write the bundled visualization fallback JSON files.
    VIZ_PUBLIC_DIR.mkdir(parents=True, exist_ok=True)
    write_graph_payload(
        graph,
        payload_out=str(VIZ_PUBLIC_DIR / "graph_payload.json"),
        nodes_out=str(VIZ_PUBLIC_DIR / "nodes_rf.json"),
        edges_out=str(VIZ_PUBLIC_DIR / "edges_rf.json"),
        time_out=str(VIZ_PUBLIC_DIR / "time_rf.json"),
    )
    print("\nWrote graph payload files into react_flow_viz/public")


if __name__ == "__main__":
    main()
