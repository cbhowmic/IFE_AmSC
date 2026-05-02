from rhino_surrogate_runtime import RhinoSurrogate
from graph_builder import (
    build_graph_from_simulation,
    write_graph_payload,
    build_react_prediction_payload,
)
from simulation_lookup import SimulationLookup


def main():
    # Physical inputs
    inputs = {
        "I0_SD": 100.0,
        "Ndotminus": 500.0,
        "beta": 0.1,
    }

    # Load surrogate
    surrogate = RhinoSurrogate("rhino_surrogate.pt", device="cpu")

    # Predict surrogate outputs
    prediction_display = surrogate.predict_with_display_names(inputs)
    print("\nSurrogate prediction:")
    print(prediction_display)

    # Load nearest-simulation lookup
    lookup = SimulationLookup(
        index_path="simulation_index.json",
        artifact_path="rhino_surrogate.pt",
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

    # Write debug JSON files locally
    write_graph_payload(graph, payload_out="graph_payload.json")
    print("\nWrote graph_payload.json, nodes_rf.json, edges_rf.json, time_rf.json")


if __name__ == "__main__":
    main()