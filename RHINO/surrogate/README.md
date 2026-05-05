# RHINO Surrogate Demo

This folder contains the live demo workflow:

1. take two physical inputs, `Ndotminus` and `beta`
2. run `rhino_surrogate.pt` through `rhino_surrogate_runtime.py`
3. find the nearest archived BP5 simulation
4. build the React Flow graph payload
5. optionally start the local visualization app

The old startup inventory input, `I0_SD`, is no longer used by the surrogate or
MCP tools.

## Demo Setup Notes

There are two separate dependency systems in this demo.

Python dependencies are listed in `RHINO/surrogate/requirements.txt`:

```bash
conda activate scspdemo
pip install -r RHINO/surrogate/requirements.txt
```

The visualization dependencies are JavaScript/React dependencies. They do not go
in `requirements.txt`; they are listed in
`RHINO/surrogate/react_flow_viz/package.json` and pinned by
`RHINO/surrogate/react_flow_viz/package-lock.json`.

Install them from the visualization folder:

```bash
cd RHINO/surrogate/react_flow_viz
npm ci
```

Use `npm ci` on the demo computer because it installs exactly the versions in
`package-lock.json`, which is more reproducible than a fresh `npm install`.

## Portable Data Layout

Keep the downloaded BP5 simulation data in:

```text
RHINO/data/surrogate_bp_output/
```

The current `simulation_index.json` may contain old absolute paths from another
machine. That is okay: `SimulationLookup` resolves those indexed paths into the
local `RHINO/data/surrogate_bp_output` folder at runtime. If the demo computer
uses a different data location, set:

```bash
export RHINO_SIM_DATA_ROOT=/path/to/surrogate_bp_output
```

## Smoke Tests

From `RHINO/surrogate`:

```bash
conda activate scspdemo
python test_local.py
python test_client.py
```

`test_client.py` uses an in-process FastMCP client, so it does not need port
8000 and will not fail if another local process is already using that port.

## Demo Launch

Use two terminals.

Terminal 1 starts the MCP server:

```bash
cd RHINO/surrogate
conda activate scspdemo
python mcp_server_rhino.py
```

The MCP server defaults to `http://127.0.0.1:8000/mcp`. If port 8000 is already
busy, use another port:

```bash
RHINO_MCP_PORT=8001 python mcp_server_rhino.py
```

Terminal 2 sends demo inputs to the MCP server:

```bash
cd RHINO/surrogate
conda activate scspdemo
python demo_client.py --Ndotminus 500 --beta 0.1
```

Calling `build_graph_for_nearest_simulation` keeps the latest graph payload in
the MCP server and exposes it through these HTTP routes:

```text
http://127.0.0.1:8000/graph_payload.json
http://127.0.0.1:8000/nodes_rf.json
http://127.0.0.1:8000/edges_rf.json
http://127.0.0.1:8000/time_rf.json
```

The React Flow app starts at `http://127.0.0.1:5173` with a `dataBaseUrl`
query parameter pointing back to the MCP server, so it reads the live graph
directly from the server. If no live graph has been generated yet, the MCP
routes use the existing JSON files as fallback data.

The visualization panel shows the archived simulation being visualized using
the BP5 run folder name, for example `00-02-22.bp5`. The surrogate panel compares
each surrogate prediction against that nearest simulation's stored output, with
absolute and relative error shown beside the two values.

When new inputs are sent through `demo_client.py` or another MCP client, the
server increments a graph revision counter. The React Flow app checks
`http://127.0.0.1:8000/graph_version.json` periodically and reloads the live
graph automatically when that revision changes.

To also write fresh JSON files into `react_flow_viz/public`, pass:

```bash
python demo_client.py --Ndotminus 500 --beta 0.1 --write-json-files
```

---

# RHINO Module

This module processes RHINO simulation outputs and converts them into AI-ready ADIOS datasets following the project schema.

This repository contains Python scripts to:
1) **Write** RHINO fuel-cycle inventory into **OpenPMD** format using **ADIOS2 BP5**, and  
2) **Read**, validate, and **plot** time-series and steady-state inventories from the produced `.bp5` file.

The scripts are designed for RHINO inventory outputs stored as Pandas pickles (`.pkl`) and follow a fixed canonical list of plant subsystems. A given run may only provide a subset of those subsystems; the writer maps provided subsystems into the canonical ordering, and the reader plots only the subsystems included in the run.

---

## Directory Structure

```text
IFE_AmSC/
├── RHINO/
    ├── scripts
        ├── rhinoWriteADIOS.py        # Write RHINO inventory data to ADIOS2 
        ├── rhinoReadADIOS.py         # Read, validate, and plot RHINO inventory data
    ├── README.md
    ├── Data/                         # User places RHINO data here
    │   └── ...
    ├── output/                       # ADIOS2 BP5 files written here
    └── plots/                        # Generated plots saved here
```

## 1) Obtain RHINO data
RHINO simulation data are available from the following SharePoint location:

https://doellnl.sharepoint.com/teams/FESAmSC-IFEPPDesignDig.Twin/Shared%20Documents/Forms/AllItems.aspx?id=%2Fteams%2FFESAmSC%2DIFEPPDesignDig%2ETwin%2FShared%20Documents%2FGeneral%2FData&viewid=5d1dddd2%2D1bbe%2D44ff%2D818c%2D39dd2133c599

The directory contains data from multiple RHINO runs. Each individual run should include the following five required files:

### Time-series inventories
- `*_Generic_FuelCycle_T.pkl`
- `*_Generic_FuelCycle_D.pkl`

These files contain time-dependent inventories for Tritium (T) and Deuterium (D) across RHINO plant subsystems.

### Steady-state inventories
- `*_Generic_FuelCycle_T_SteadyState.pkl`
- `*_Generic_FuelCycle_D_SteadyState.pkl`

These files represent the steady-state (final) inventories corresponding to the time-series data.

### Metadata
- `*_meta.pkl`

This file contains run-level metadata describing the RHINO simulation configuration.


## 2) Place RHINO data
Place the downloaded data inside `RHINO/Data/`.

Each RHINO run should be placed in its own subdirectory containing the five required files.


## 3) Environment Setup
This module uses the project-wide Conda environment:
```bash
conda activate IFE_AmSC
```
If the environment has not yet been created, build it from the repository root:
```bash
conda env create -f environment.yml
conda activate IFE_AmSC
```

## 4) Running the Scripts

### 4.1) Write ADIOS file from RHINO data
From inside `RHINO/scripts/`:
```bash
python rhinoWriteADIOS.py
```
This script:
- Reads `.pkl` inventory files from `RHINO/Data/`
- Maps subsystems into canonical ordering
- Writes an OpenPMD-compliant ADIOS2 BP5 file
- Stores the output in `RHINO/output/`.


### 4.2) Read, Validate, and Plot
From inside `RHINO/scripts/`:
```bash
python rhinoReadADIOS.py
```
This script:
- Reads the generated BP5 file from RHINO/output/
- Validates data consistency
- Generates time-series and steady-state plots
- Saves plots into `RHINO/plots/`.
