from fastmcp import FastMCP
from starlette.responses import JSONResponse
from rhino_surrogate_runtime import RhinoSurrogate
from simulation_lookup import SimulationLookup 
from graph_builder import build_graph_from_simulation, build_react_prediction_payload

mcp = FastMCP("RHINO Surrogate")

# Load once when the server starts
surrogate = RhinoSurrogate("rhino_surrogate.pt", device="cpu")

lookup = SimulationLookup(
    index_path="simulation_index.json",
    artifact_path="rhino_surrogate.pt",
)
    
@mcp.tool
def predict_rhino_surrogate(
    I0_SD: float,
    Ndotminus: float,
    beta: float,
) -> dict[str, float]:
    """
    Predict RHINO surrogate outputs from 3 physical input values.

    Inputs:
    - I0_SD: startup inventory in grams [g]
    - Ndotminus: tritium burning rate in grams per day [g/d]
    - beta: burn fraction as a unitless fraction [-]

    The inputs must be physical values, not normalized values.
    Normalization is handled internally by the surrogate runtime.

    Returns:
    A dictionary mapping human-readable output names to predicted physical values.
    """
    return surrogate.predict_with_display_names(
        {
            "I0_SD": I0_SD,
            "Ndotminus": Ndotminus,
            "beta": beta,
        }
    )

@mcp.tool
def find_nearest_simulation(
    I0_SD: float,
    Ndotminus: float,
    beta: float,
) -> dict:
    """
    Find the nearest archived simulation in normalized input space.

    Inputs are physical values:
    - I0_SD in grams [g]
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
            "I0_SD": I0_SD,
            "Ndotminus": Ndotminus,
            "beta": beta,
        }
    )

@mcp.tool
def predict_and_compare_to_nearest_simulation(
    I0_SD: float,
    Ndotminus: float,
    beta: float,
) -> dict:
    """
    Predict RHINO surrogate outputs from physical inputs, find the nearest archived
    simulation in normalized input space, and compare the prediction to that
    simulation's outputs.
    
    Inputs:
    - I0_SD: startup inventory [g]
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
        "I0_SD": I0_SD,
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
    I0_SD: float,
    Ndotminus: float,
    beta: float,
) -> dict:
    """
    Find the nearest archived simulation for the given physical inputs and
    return React Flow graph data for that simulation, including surrogate
    prediction values formatted for the frontend.
    """
    query_inputs = {
        "I0_SD": I0_SD,
        "Ndotminus": Ndotminus,
        "beta": beta,
    }

    nearest = lookup.find_nearest(query_inputs)
    prediction_display = surrogate.predict_with_display_names(query_inputs)
    react_predictions = build_react_prediction_payload(prediction_display)

    graph = build_graph_from_simulation(
        nearest["path"],
        surrogate_predictions=react_predictions,
    )

    return {
        "nearest_simulation": nearest,
        "graph": graph,
    }

@mcp.custom_route("/health", methods=["GET"])
async def health_check(request):
    return JSONResponse({"status": "healthy", "service": "rhino-mcp"})

if __name__ == "__main__":
    mcp.run(transport="http", host="0.0.0.0", port=8000)