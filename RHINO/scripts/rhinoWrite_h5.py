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

# Labels for T and D
labels_T = json.dumps([str(x) for x in T_ts_df.index.tolist()])
labels_D = json.dumps([str(x) for x in D_ts_df.index.tolist()])

# Metadata
dt = float(meta["dt"]) if "dt" in meta else None

# Block parameters
block_len = 100
nblocks = math.ceil(Nt1 / block_len)

# Create an openPMD series
OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhino_block.h5"
series = io.Series(OUTPUT_PATH, io.Access.create)
print("Creating RHINO data in openPMD format (HDF5)...")
print("Block length:", block_len, "Number of blocks:", nblocks)

# Top-level metadata
series.author = "Chandreyee Bhowmick <ccb@ornl.gov>"
series.date = "2026-01-12"
series.set_attribute("schema", "rhino-1.0")
series.set_attribute("application", "RHINO: Fusion Pilot Plant fuel cycle")

# RHINO specific metadata
series.set_attribute("layout", "block_iterations")
series.set_attribute("block_len", int(block_len))
series.set_attribute("Nt_total", int(Nt1))

series.set_attribute("labels/T_rows", labels_T)
series.set_attribute("labels/D_rows", labels_D)
series.set_attribute("meta/dt_unit", "day")
for k, v in meta.items():
    series.set_attribute(f"meta/{k}", v)


# Record time-series data in blocks
SECONDS_PER_DAY = 86400.0
for b in range(nblocks):
    t0 = b * block_len
    t1 = min((b + 1) * block_len, Nt1)
    this_len = t1 - t0

    it = series.iterations[b]

    # Iteration metadata describing this block
    it.set_attribute("description", "Time-series inventory data")
    it.set_attribute("block_start_step", int(t0))
    it.set_attribute("block_len", int(this_len))

    # Time info
    if dt is not None:
        it.time = float(t0 * dt)
        it.dt = float(dt)
        it.time_unit_SI = SECONDS_PER_DAY

    # Tritium block info
    meshT = it.meshes["T"]
    compT = meshT[io.Mesh_Record_Component.SCALAR]
    dataT = np.ascontiguousarray(T_ts[:, t0:t1])

    compT.reset_dataset(io.Dataset(dataT.dtype, dataT.shape))
    compT.store_chunk(dataT)
    
    meshT.unit_dimension = {io.Unit_Dimension.M: 1}
    compT.unit_SI = 1e-3

    # Deuterium block info
    meshD = it.meshes["D"]
    compD = meshD[io.Mesh_Record_Component.SCALAR]
    dataD = np.ascontiguousarray(D_ts[:, t0:t1])
    
    compD.reset_dataset(io.Dataset(dataD.dtype, dataD.shape))
    compD.store_chunk(dataD)
    
    meshD.unit_dimension = {io.Unit_Dimension.M: 1}
    compD.unit_SI = 1e-3

    # Meash metadata
    meshT.set_attribute("axisLabels", ["subsystem", "time_in_block"])
    meshD.set_attribute("axisLabels", ["subsystem", "time_in_block"])
    meshT.set_attribute("axisUnits", ["subsystem index", "day"])
    meshD.set_attribute("axisUnits", ["subsystem index", "day"])
    meshT.set_attribute("rowLabels", labels_T)
    meshD.set_attribute("rowLabels", labels_D)

    if dt is not None:
        meshT.grid_spacing = [1.0, float(dt)]
        meshT.grid_global_offset = [0.0, 0.0]
        meshD.grid_spacing = [1.0, float(dt)]
        meshD.grid_global_offset = [0.0, 0.0]
    
    it.close()

# Record steady-state data in the last iteration
it_ss = series.iterations[nblocks]
it_ss.set_attribute("description", "Steady-state inventory data")
if dt is not None:
    it_ss.time = float((Nt1 - 1) * dt)
    it_ss.dt = float(dt)
    it_ss.time_unit_SI = SECONDS_PER_DAY

# Tritium steady-state (1D)
meshTss = it_ss.meshes["T_steady"]
compTss = meshTss[io.Mesh_Record_Component.SCALAR]
dataTss = np.ascontiguousarray(T_ss)   
compTss.reset_dataset(io.Dataset(dataTss.dtype, dataTss.shape))
compTss.store_chunk(dataTss)
meshTss.unit_dimension = {io.Unit_Dimension.M: 1}
compTss.unit_SI = 1e-3  

# Deuterium steady-state (1D)
meshDss = it_ss.meshes["D_steady"]
compDss = meshDss[io.Mesh_Record_Component.SCALAR]
dataDss = np.ascontiguousarray(D_ss)   
compDss.reset_dataset(io.Dataset(dataDss.dtype, dataDss.shape))
compDss.store_chunk(dataDss)
meshDss.unit_dimension = {io.Unit_Dimension.M: 1}
compDss.unit_SI = 1e-3  

# Steady-state mesh metadata
meshTss.set_attribute("axisLabels", ["subsystem"])  
meshDss.set_attribute("axisLabels", ["subsystem"])
meshTss.set_attribute("axisUnits", ["subsystem index"])
meshDss.set_attribute("axisUnits", ["subsystem index"])
meshTss.set_attribute("rowLabels", labels_T)
meshDss.set_attribute("rowLabels", labels_D)

it_ss.close()

# Close the series
series.flush()
series.close()
print("RHINO data written in openPMD block-iteration format (HDF5).")
print("Output:", OUTPUT_PATH)