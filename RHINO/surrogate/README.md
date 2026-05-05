# RHINO Surrogate Demo

## TLDR

Run the RHINO surrogate for `Ndotminus` and `beta`, match the nearest archived BP5 simulation, and visualize the result with the local React Flow demo. The demo assumes the directory structure below and requires two terminals.

Installation:
```bash 
git clone https://github.com/cbhowmic/IFE_AmSC.git
cd IFE_AmSC
git checkout surrogate
conda create -y --name scspdemo -c conda-forge python nodejs
conda activate scspdemo
pip install -r RHINO/surrogate/requirements.txt
npm --prefix RHINO/surrogate/react_flow_viz ci
unzip RHINO/surrogate/data/surrogate_bp_output.zip -d RHINO/surrogate/data/
```

From the same terminal:

```bash
RHINO_MCP_PORT=8001 python RHINO/surrogate/mcp_server_rhino.py
```

Open a second terminal:

```bash
conda activate scspdemo
python RHINO/surrogate/demo_client.py --Ndotminus 500 --beta 0.1 --mcp-url http://127.0.0.1:8001/mcp
```

Click on the link to open the visualization in the browser. 


## Details 
This folder contains the workflow for the live demo:

1. take two physical inputs, `Ndotminus` and `beta`
2. run `rhino_surrogate.pt` through `rhino_surrogate_runtime.py`
3. find the nearest archived BP5 simulation from **local folder** `data/`
4. build the React Flow graph payload
5. start the local visualization app

## Directory Structure

From `RHINO/surrogate`:

```text
$ tree -L 3 -I 'node_modules|__pycache__|dist'
.
├── build_graph.py
├── build_simulation_index.ipynb
├── build_simulation_index.py
├── data
│   ├── rhino_surrogate.pt
│   ├── simulation_index.json
│   ├── surrogate_bp_output
│   │   ├── 2026-04-29
│   │   ├── 2026-04-30
│   │   └── 2026-05-01
│   └── surrogate_bp_output.zip
├── demo_client.py
├── graph_builder.py
├── mcp_server_rhino.py
├── react_flow_viz
│   ├── index.html
│   ├── package.json
│   ├── package-lock.json
│   ├── public
│   │   ├── edges_rf.json
│   │   ├── favicon.svg
│   │   ├── graph_payload.json
│   │   ├── icons.svg
│   │   ├── nodes_rf.json
│   │   └── time_rf.json
│   ├── src
│   │   ├── App.css
│   │   ├── App.jsx
│   │   ├── assets
│   │   ├── index.css
│   │   ├── main.jsx
│   │   └── README.md
│   └── vite.config.js
├── README.md
├── requirements.txt
├── rhinoSurrogate.ipynb
├── rhinoSurrogate.py
├── rhino_surrogate_runtime.py
├── simulation_lookup.py
├── test_client.py
└── test_local.py
```

## Demo Setup

From the repository root, create one Conda environment with Python plus Node.js.
The `nodejs` Conda package also provides `npm`:

```bash
conda create -y --name scspdemo -c conda-forge python nodejs
conda activate scspdemo
```

If the environment already exists, add Node.js to it with:

```bash
conda activate scspdemo
conda install -y -c conda-forge nodejs
```

Install the Python dependencies:

```bash
pip install -r RHINO/surrogate/requirements.txt
```

Check that Node.js and `npm` are available:

```bash
node --version
npm --version
```

Install the React Flow visualization dependencies:

```bash
npm --prefix RHINO/surrogate/react_flow_viz ci
```

Use `npm ci` on the demo computer because it installs exactly the versions in
`package-lock.json`, which is more reproducible than a fresh `npm install`.

Keep the downloaded BP5 simulation data in:

```text
RHINO/surrogate/data/surrogate_bp_output/
```

If you have the zipped BP5 data folder, unzip it from the repository root:

```bash
unzip RHINO/surrogate/data/surrogate_bp_output.zip -d RHINO/surrogate/data/
```

The current `simulation_index.json` contains absolute paths from NERSC.
That is okay: `SimulationLookup` resolves those indexed paths into the
local `RHINO/surrogate/data/surrogate_bp_output` folder at runtime. If the demo
computer uses a different data location, set:

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

## Input Parameters

Pass the two physical inputs with the demo client flags below. Use physical
values, not normalized values; the surrogate runtime handles normalization
internally.


| Physical input | CLI flag      | Units                |
| -------------- | ------------- | -------------------- |
| Burning rate   | `--Ndotminus` | grams per day, `g/d` |
| Burn fraction  | `--beta`      | unitless fraction    |


Example:

```bash
python RHINO/surrogate/demo_client.py --Ndotminus 500 --beta 0.1
```

This command requires the MCP server to be running first; see the two-terminal
launch instructions below.

The same parameter names, `Ndotminus` and `beta`, are used by the MCP tools
`predict_rhino_surrogate`, `find_nearest_simulation`,
`predict_and_compare_to_nearest_simulation`, and
`build_graph_for_nearest_simulation`.

## Demo Launch

Use two terminals from the repository root.

Terminal 1 starts the MCP server:

```bash
conda activate scspdemo
python RHINO/surrogate/mcp_server_rhino.py
```

The MCP server defaults to `http://127.0.0.1:8000/mcp`. If port 8000 is already
busy, use another port:

```bash
RHINO_MCP_PORT=8001 python RHINO/surrogate/mcp_server_rhino.py
```

Terminal 2 sends demo inputs to the MCP server:

```bash
conda activate scspdemo
python RHINO/surrogate/demo_client.py --Ndotminus 500 --beta 0.1
```

If the server is running on port 8001, point the client to that port:

```bash
python RHINO/surrogate/demo_client.py --Ndotminus 500 --beta 0.1 --mcp-url http://127.0.0.1:8001/mcp
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
routes use the bundled fallback JSON files in `react_flow_viz/public`.

The visualization panel shows the archived simulation being visualized using
the BP5 run folder name, for example `00-02-22.bp5`. The surrogate panel compares
each surrogate prediction against that nearest simulation's stored output, with
absolute and relative error shown beside the two values.

When new inputs are sent through `demo_client.py` or another MCP client, the
server increments a graph revision counter. The React Flow app checks
`http://127.0.0.1:8000/graph_version.json` periodically and reloads the live
graph automatically when that revision changes.

To also refresh the bundled fallback JSON files in `react_flow_viz/public`, pass:

```bash
python RHINO/surrogate/demo_client.py --Ndotminus 500 --beta 0.1 --write-json-files
```
