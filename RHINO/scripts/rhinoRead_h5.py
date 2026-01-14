import openpmd_api as io
import json
import numpy as np
import matplotlib.pyplot as plt
import os
import math

H5_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/rhino_block.h5"
OUTDIR = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/Plots"
os.makedirs(OUTDIR, exist_ok=True)

# Block parameters
nblocks = len(series.iterations)
downsample = 100

series = io.Series(H5_PATH, io.Access.read_only)

# Print metadata
block_len = int(series.get_attribute("rhino/block_len")) if series.contains_attribute("rhino/block_len") else None
Nt_total = int(series.get_attribute("rhino/Nt_total")) if series.contains_attribute("rhino/Nt_total") else None

dt = None
if series.contains_attribute("rhino/meta/dt"):
    dt = float(series.get_attribute("rhino/meta/dt"))

print("Opened:", H5_PATH)
print("openPMD version:", series.openPMD)
print("block_len:", block_len, "Nt_total:", Nt_total, "dt:", dt)

# Labels (optional)
labels_T = None
if series.contains_attribute("rhino/labels/T_rows_json"):
    labels_T = json.loads(series.get_attribute("rhino/labels/T_rows_json"))
    print("T labels:", len(labels_T), "(first 5)", labels_T[:5])

# Choose blocks to read
iter_keys = sorted(series.iterations)
blocks = iter_keys[:min(nblocks, len(iter_keys))]
print(f"Reading {len(blocks)} blocks: {blocks[0]} .. {blocks[-1]}")

def get_attr(obj, name):
    try:
        return obj.get_attribute(name)
    except Exception:
        return None

print("\n=== Series metadata ===")
for key in [
    "schema",
    "schema/version",
    "application",
    "run/generated_utc",
    "rhino/meta/dt",
    "rhino/meta/units",
]:
    val = get_attr(series, key)
    if val is not None:
        print(f"{key}: {val}")


# Read blocks and build timeseries
times = []
T_total = []
D_total = []

for b in blocks:
    it = series.iterations[b]

    # block start step + length (stored on iteration)
    start_step = it.get_attribute("rhino/block_start_step") if it.contains_attribute("rhino/block_start_step") else b * (block_len or 0)
    this_len = it.get_attribute("rhino/block_len") if it.contains_attribute("rhino/block_len") else (block_len or 0)
    start_step = int(start_step)
    this_len = int(this_len)

    # Read 2D blocks: (nSubsystem, this_len)
    T_block = it.meshes["T"][io.Mesh_Record_Component.SCALAR].load_chunk()
    D_block = it.meshes["D"][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()

    # Compute totals across subsystem dimension -> 1D over time-in-block
    T_sum = np.sum(T_block, axis=0)
    D_sum = np.sum(D_block, axis=0)

    # Make time vector for this block
    if dt is not None:
        tvec = (start_step + np.arange(this_len)) * dt
    else:
        tvec = start_step + np.arange(this_len)

    # Append
    if downsample > 1:
        tvec = tvec[::downsample]
        T_sum = T_sum[::downsample]
        D_sum = D_sum[::downsample]

    times.append(tvec)
    T_total.append(T_sum)
    D_total.append(D_sum)

    it.close()

# Concatenate blocks into one long array
times = np.concatenate(times) if len(times) else np.array([])
T_total = np.concatenate(T_total) if len(T_total) else np.array([])
D_total = np.concatenate(D_total) if len(D_total) else np.array([])

series.close()

print("Built timeseries points:", len(times))
print("Time range:", (times[0], times[-1]) if len(times) else None)

# Plot timeseries
plt.figure()
plt.plot(times, T_total, label="Total Tritium")
plt.plot(times, D_total, label="Total Deuterium")
plt.xlabel("Time [s]" if dt is not None else "Step")
plt.ylabel("Total Inventory [kg]")
plt.title(f"RHINO totals from blocks (first {len(blocks)} blocks, downsample={downsample})")
plt.grid(True)
plt.legend()

png_path = os.path.join(OUTDIR, f"rhino_first_{len(blocks)}_blocks.png")
pdf_path = os.path.join(OUTDIR, f"rhino_first_{len(blocks)}_blocks.pdf")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.show()

print("Saved:")
print(" ", png_path)
print(" ", pdf_path)