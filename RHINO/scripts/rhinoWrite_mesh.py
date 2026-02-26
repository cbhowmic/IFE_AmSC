# This file reads RHINO output data (pkl files) and writes it in an openPMD series
# Uses ADIOS2 as a backend (.bp5) with variable-based iteration encoding
# Stores data in mesh records (i.e. DataFrame-like)
# Data is stored in a shared folder on NERSC
# /global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino

import openpmd_api as io
import numpy as np
import pandas as pd
import json
import os

# Configuration
OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhino_mesh.bp5"
SECONDS_PER_DAY = 86400.0


# Load RHINO Data
T_ts_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T.pkl")
D_ts_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D.pkl")
meta_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_meta.pkl")
T_ss_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T_SteadyState.pkl")
D_ss_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D_SteadyState.pkl")

meta = meta_df[0].to_dict()
for k, v in list(meta.items()):
    if isinstance(v, np.generic):
        meta[k] = v.item()

# Canonical subsystem ordering
labels_subsystem = [
    "Storage_Delivery", "Fueling", "Fusion_Chamber_Pump", "Pd_Cleanup",
    "Protium_Removal", "Exhaust_Processing", "Gas_Detrit", "Water_Detrit",
    "Glovebox", "Stack", "Isotope_Seperation", "Blanket_Extraction",
    "Heat_Exchanger", "Power_Conv_Loop", "Vent_Detrit", "Blanket",
    "Decay_Box", "Stack_Box", "Burn_Box", "Gen_Box", "Uptake_Box"
]
canon = {name: i for i, name in enumerate(labels_subsystem)}

# Build canonical arrays matching with subsystem ordering
T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64))
nSubsystems, Nt = T_ts.shape

D_ts = np.zeros((nSubsystems, Nt), dtype=np.float64)
for name, row in D_ts_df.iterrows():
    D_ts[canon[str(name)], :] = row.to_numpy(dtype=np.float64)

T_ss = T_ss_df.iloc[:, 0].to_numpy(dtype=np.float64)

D_ss = np.zeros((nSubsystems,), dtype=np.float64)
for name, row in D_ss_df.iterrows():
    D_ss[canon[str(name)]] = float(row.iloc[0])

dt = float(meta["dt"]) if "dt" in meta else None


# Create openPMD Series
adios2_cfg = r'''
{
  "iteration_encoding": "variable_based",
  "adios2": {
    "engine": {
      "type": "bp5",
      "parameters": {
        "StatsLevel": "0"
      }
    }
  }
}
'''
series = io.Series(OUTPUT_PATH, io.Access_Type.create_linear, adios2_cfg)
print("Writing RHINO data in openPMD meshes...")


# Series-level metadata
series.author = "Holly Flynn <Holly.Flynn@srnl.doe.gov>"
series.date = "2025-12-16"

series.set_attribute("software", "RHINO")
series.set_attribute("softwareDescription", "RHINO: Fusion Pilot Plant fuel cycle")
series.set_attribute("softwareVersion", "1.0")
series.set_attribute("schema", "OpenPMD")
series.set_attribute("schemaVersion", "0.17.0")

series.set_attribute("subsystemLabels", labels_subsystem)

if dt is not None:
    series.set_attribute("timeStep", float(dt))
    series.set_attribute("timeUnit", "day")
    series.set_attribute("startTime", 0.0)
    series.set_attribute("endTime", float((Nt - 1) * dt))

for k, v in meta.items():
    if isinstance(v, (int, float, str, bool)):
        series.set_attribute(f"meta/{k}", v)

series.flush()

# Iterations
it = series.iterations[0]
it.set_attribute("description", "Fuel-cycle inventory time series")

if dt is not None:
    it.time = 0.0
    it.dt = float(dt)
    it.time_unit_SI = SECONDS_PER_DAY

series.meshes_path = "inventory"


# Tritium time-series mesh
meshT = it.meshes["Tritium"]
compT = meshT[io.Mesh_Record_Component.SCALAR]

compT.reset_dataset(io.Dataset(T_ts.dtype, T_ts.shape))
compT.store_chunk(T_ts)

meshT.unit_dimension = {io.Unit_Dimension.M: 1}
compT.unit_SI = 1e-3  

meshT.set_attribute("axisLabels", ["subsystem"])


# Deuterium time-series mesh
meshD = it.meshes["Deuterium"]
compD = meshD[io.Mesh_Record_Component.SCALAR]

compD.reset_dataset(io.Dataset(D_ts.dtype, D_ts.shape))
compD.store_chunk(D_ts)

meshD.unit_dimension = {io.Unit_Dimension.M: 1}
compD.unit_SI = 1e-3

meshD.set_attribute("axisLabels", ["subsystem"])


# Steady-state meshes 
meshTss = it.meshes["Tritium_steady"]
compTss = meshTss[io.Mesh_Record_Component.SCALAR]
compTss.reset_dataset(io.Dataset(T_ss.dtype, T_ss.shape))
compTss.store_chunk(T_ss)

meshTss.unit_dimension = {io.Unit_Dimension.M: 1}
compTss.unit_SI = 1e-3
meshTss.set_attribute("axisLabels", ["subsystem"])

meshDss = it.meshes["Deuterium_steady"]
compDss = meshDss[io.Mesh_Record_Component.SCALAR]
compDss.reset_dataset(io.Dataset(D_ss.dtype, D_ss.shape))
compDss.store_chunk(D_ss)

meshDss.unit_dimension = {io.Unit_Dimension.M: 1}
compDss.unit_SI = 1e-3
meshDss.set_attribute("axisLabels", ["subsystem"])


it.close()
series.close()

print("RHINO mesh data written successfully.")
print("Output:", OUTPUT_PATH)
