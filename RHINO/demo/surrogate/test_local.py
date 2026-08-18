from pathlib import Path

from RHINO.demo.surrogate.rhino_surrogate_runtime import RhinoSurrogate

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

print("Testing surrogate")

surrogate = RhinoSurrogate(DATA_DIR / "rhino_surrogate.pt", device="cpu")

result = surrogate.predict_from_named_args(
    Ndotminus=500.0,
    beta=0.1,
)

print("Surrogate prediction from named args")
print(result)

result = surrogate.predict_with_display_names({
    "Ndotminus": 500.0,
    "beta": 0.1,
})

print("Surrogate prediction from display names")
print(result)

print("#################")
print("Testing simulation lookup")

from RHINO.demo.surrogate.simulation_lookup import SimulationLookup

lookup = SimulationLookup(
    index_path=DATA_DIR / "simulation_index.json",
    artifact_path=DATA_DIR / "rhino_surrogate.pt",
)

result = lookup.find_nearest(
    {
        "Ndotminus": 273.,
        "beta": 0.16,
    }
)
print("Testing simulation lookup")
print(result["simulation_id"])
print(result["distance_normalized"])
print(result["inputs"])
print(result["outputs"])

result = lookup.find_nearest(
    {
        "Ndotminus": 272.,
        "beta": 0.15,
    }
)
print(result["simulation_id"])
print(result["distance_normalized"])
print(result["inputs"])
print(result["outputs"])


record0 = lookup.records[0]

result0 = lookup.find_nearest(record0["inputs"])

print("expected:", record0["simulation_id"])
print("found:   ", result0["simulation_id"])
print("distance:", result0["distance_normalized"])

print("#################")
print("Testing graph builder")


from RHINO.demo.surrogate.graph_builder import build_graph_from_simulation

sim_path = DATA_DIR / "surrogate_bp_output" / "2026-04-30" / "10-39-07.bp5"
graph = build_graph_from_simulation(sim_path)

print(graph.keys())
