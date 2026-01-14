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

# Extract metadata as dictionary
meta = meta_df[0].to_dict()
for k, v in list(meta.items()):
    if isinstance(v, np.generic):
        meta[k] = v.item()

# Time-series arrays
T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64, copy=True))
D_ts = np.ascontiguousarray(D_ts_df.to_numpy(dtype=np.float64, copy=True))
nT, Nt = T_ts.shape           # nT=21
nD, Nt2 = D_ts.shape          # nD=15
assert Nt == Nt2, "Time dimension mismatch between T and D"

# Labels for T and D
labels_T = json.dumps([str(x) for x in T_ts_df.index.tolist()])
labels_D = json.dumps([str(x) for x in D_ts_df.index.tolist()])

# Metadata
dt = float(meta["dt"]) if "dt" in meta else None

# Block parameters
block_len = 100
nblocks = math.ceil(Nt / block_len)

# Create an openPMD series
OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/rhino_block.h5"
series = io.Series(OUTPUT_PATH, io.Access.create)
print("Creating RHINO data in openPMD format (HDF5)...")
print("Block length:", block_len, "Number of blocks:", nblocks)

# Series attributes
series.author = "Chandreyee Bhowmick <ccb@ornl.gov>"
series.date = "2026-01-12"
series.set_attribute("schema", "rhino-1.0")
series.set_attribute("application", "RHINO: Fusion Pilot Plant fuel cycle")

# Metadata
series.set_attribute("rhino/layout", "block_iterations")
series.set_attribute("rhino/block_len", int(block_len))
series.set_attribute("rhino/Nt_total", int(Nt))

series.set_attribute("rhino/labels/T_rows_json", labels_T)
series.set_attribute("rhino/labels/D_rows_json", labels_D)
for k, v in meta.items():
    series.set_attribute(f"rhino/meta/{k}", v)


# Record time-series data in blocks
for b in range(nblocks):
    t0 = b * block_len
    t1 = min((b + 1) * block_len, Nt)
    this_len = t1 - t0

    it = series.iterations[b]

    # Iteration metadata describing this block
    it.set_attribute("rhino/block_start_step", int(t0))
    it.set_attribute("rhino/block_len", int(this_len))

    # Time info
    if dt is not None:
        it.time = float(t0 * dt)
        it.dt = float(dt)
        it.time_unit_SI = 1.0

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

    it.close()


# Close the series
series.flush()
series.close()
print("RHINO data written in openPMD block-iteration format (HDF5).")
print("Output:", OUTPUT_PATH)