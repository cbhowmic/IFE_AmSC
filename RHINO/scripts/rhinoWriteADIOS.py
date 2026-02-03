import openpmd_api as io
import numpy as np
import os
import pandas as pd
import json
import math

# Load the RHINO timeseries data
T_ts_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T.pkl")
D_ts_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D.pkl")
meta_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_meta.pkl")     
T_ss_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T_SteadyState.pkl")   
D_ss_df = pd.read_pickle("Data/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D_SteadyState.pkl")   

# Extract metadata as dictionary
meta = meta_df[0].to_dict()
for k, v in list(meta.items()):
    if isinstance(v, np.generic):
        meta[k] = v.item()

# Time-series arrays
T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64, copy=True))
D_ts = np.ascontiguousarray(D_ts_df.to_numpy(dtype=np.float64, copy=True))
nT, Nt1 = T_ts.shape           # nT=21
nD, Nt2 = D_ts.shape          # nD=15
assert Nt1 == Nt2, "Time dimension mismatch between T and D"

# Steady-state arrays
T_ss = T_ss_df.iloc[:, 0].to_numpy(dtype=np.float64)
D_ss = D_ss_df.iloc[:, 0].to_numpy(dtype=np.float64)
assert T_ss.shape[0] == nT, "Steady-state T size mismatch"
assert D_ss.shape[0] == nD, "Steady-state D size mismatch"

# Last columns of T and D should match steady-state arrays
last_T = T_ts[:, -1]
allclose_T = np.allclose(last_T, T_ss, atol=0, rtol=0)
max_abs_T = float(np.max(np.abs(last_T - T_ss)))
print('Exact match for T:', allclose_T)
print('Max abs diff for T:', max_abs_T)
last_D = D_ts[:, -1]
allclose_D = np.allclose(last_D, D_ss, atol=0, rtol=0)
max_abs_D = float(np.max(np.abs(last_D - D_ss)))
print('Exact match for D:', allclose_D)
print('Max abs diff for D:', max_abs_D)

# Make T and D of the same size by padding D with zeros
T_labels = [str(x) for x in T_ts_df.index.tolist()]        
D_labels = [str(x) for x in D_ts_df.index.tolist()]    
D_ts_expanded = np.zeros((nT, Nt1), dtype=np.float64)
D_ts_expanded[:nD, :] = D_ts
D_ts = np.ascontiguousarray(D_ts_expanded)

D_ss_expanded = np.zeros((nT,), dtype=np.float64)
D_ss_expanded[:nD] = D_ss
D_ss = np.ascontiguousarray(D_ss_expanded)

# Metadata
dt = float(meta["dt"]) if "dt" in meta else None
print("dt:", dt)

# Create an openPMD series
OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhinoADIOS.bp5"

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
print("Creating RHINO data in openPMD format (ADIOS)...")

# Top-level metadata
series.author = "Chandreyee Bhowmick <ccb@ornl.gov>"
series.date = "2026-02-03"
series.set_attribute("software", "RHINO: Fusion Pilot Plant fuel cycle")
series.set_attribute("softwareVersion", "1.0")   # TODO: update version
series.set_attribute("schema", "rhino-1.0")
series.set_attribute("OpenPMD_version", io.__version__)

labels_subsystem = ["Storage_Delivery", "Fueling", "Fusion_Chamber_Pump", "Pd_Cleanup", "Protium_Removal", "Exhaust_Processing", "Gas_Detrit", "Water_Detrit", "Glovebox", "Stack", "Isotope_Seperation", "Blanket_Extraction", "Heat_Exchanger", "Power_Conv_Loop", "Vent_Detrit", "Blanket", "Decay_Box", "Stack_Box", "Burn_Box", "Gen_Box", "Uptake_Box"]
series.set_attribute("subsystem labels", labels_subsystem)
if dt is not None:
    t_start = 0.0
    t_end = float((Nt1 - 1) * dt)
    t_duration = float(Nt1 * dt)   

    series.set_attribute("time/start time", float(t_start))
    series.set_attribute("time/end time", float(t_end))
    series.set_attribute("time/step", float(dt))
    series.set_attribute("time/unit", "day")

for k, v in meta.items():
    if isinstance(v, (int, float, str, bool, np.integer, np.floating)):
        series.set_attribute(f"meta/{k}", v)

series.flush()

# Record time-series data in blocks
SECONDS_PER_DAY = 86400.0
it = series.iterations[0]  
it.set_attribute("description", "Inventory data of various gas species in pilot power plant")

# Time info
if dt is not None:
    it.time = float(0 * dt)
    it.dt = float(dt)
    it.time_unit_SI = SECONDS_PER_DAY
    
# Tritium inventory
meshT = it.meshes["Tritium"]
compT = meshT[io.Mesh_Record_Component.SCALAR]
dataT = np.ascontiguousarray(T_ts)

compT.reset_dataset(io.Dataset(dataT.dtype, dataT.shape))
compT.store_chunk(dataT)

meshT.unit_dimension = {io.Unit_Dimension.M: 1}
compT.unit_SI = 1e-3
    
# Deuterium inventory
meshD = it.meshes["Deuterium"]
compD = meshD[io.Mesh_Record_Component.SCALAR]
dataD = np.ascontiguousarray(D_ts)

compD.reset_dataset(io.Dataset(dataD.dtype, dataD.shape))
compD.store_chunk(dataD)
    
meshD.unit_dimension = {io.Unit_Dimension.M: 1}
compD.unit_SI = 1e-3

# Meash metadata
meshT.set_attribute("axisLabels", "subsystem")
meshD.set_attribute("axisLabels", "subsystem")
    
if dt is not None:
    meshT.grid_spacing = [1.0, float(dt)]
    meshT.grid_global_offset = [0.0, 0.0]
    meshD.grid_spacing = [1.0, float(dt)]
    meshD.grid_global_offset = [0.0, 0.0]

# Tritium steady-state (1D)
meshTss = it.meshes["Tritium_steady"]
compTss = meshTss[io.Mesh_Record_Component.SCALAR]
dataTss = np.ascontiguousarray(T_ss)   
compTss.reset_dataset(io.Dataset(dataTss.dtype, dataTss.shape))
compTss.store_chunk(dataTss)
meshTss.unit_dimension = {io.Unit_Dimension.M: 1}
compTss.unit_SI = 1e-3  

# Deuterium steady-state (1D)
meshDss = it.meshes["Deuterium_steady"]
compDss = meshDss[io.Mesh_Record_Component.SCALAR]
dataDss = np.ascontiguousarray(D_ss)   
compDss.reset_dataset(io.Dataset(dataDss.dtype, dataDss.shape))
compDss.store_chunk(dataDss)
meshDss.unit_dimension = {io.Unit_Dimension.M: 1}
compDss.unit_SI = 1e-3  

# Steady-state mesh metadata
meshTss.set_attribute("axisLabels", "subsystem")  
meshDss.set_attribute("axisLabels", "subsystem")

it.close()

# Close the series
series.flush()
series.close()
print("RHINO data written in openPMD block-iteration format (ADIOS).")
print("Output:", OUTPUT_PATH)