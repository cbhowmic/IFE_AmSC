import openpmd_api as io
import json
import numpy as np
import matplotlib.pyplot as plt
import os

BP5_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhino_block.bp5"
OUTDIR = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/plots"
os.makedirs(OUTDIR, exist_ok=True)

downsample = 100

series = io.Series(BP5_PATH, io.Access.read_only)

print("\n=== File opened ===")
print("Path:", BP5_PATH)
print("openPMD version:", series.openPMD)

# --------------------------------------------------
# Metadata 
# --------------------------------------------------

# File-level metadata
block_len = int(series.get_attribute("block_len"))
Nt_total = int(series.get_attribute("Nt_total"))
dt = float(series.get_attribute("meta_dt")) if series.contains_attribute("meta_dt") else None
dt_unit = series.get_attribute("meta_dt_unit") if series.contains_attribute("meta_dt_unit") else None  
metadata_json = series.get_attribute("meta_json") if series.contains_attribute("meta_json") else None

print("\n=== File-level metadata ===")
print("block_len:", block_len)
print("Nt_total:", Nt_total)
print("dt:", dt)
print("dt_unit:", dt_unit)
print("metadata_json:", metadata_json)


# Labels 
labels_T = None
labels_D = None
if series.contains_attribute("labels_T_rows"):
    labels_T = json.loads(series.get_attribute("labels_T_rows"))
if series.contains_attribute("labels_D_rows"):
    labels_D = json.loads(series.get_attribute("labels_D_rows"))

print("T labels:", len(labels_T) if labels_T is not None else None)
print("D labels:", len(labels_D) if labels_D is not None else None)


# --------------------------------------------------
# Steady-state inventories 
# --------------------------------------------------
iter_keys = sorted(series.iterations)
steady_it = iter_keys[-1]
it_ss = series.iterations[steady_it]
T_ss = it_ss.meshes["T_steady"][io.Mesh_Record_Component.SCALAR].load_chunk()
D_ss = it_ss.meshes["D_steady"][io.Mesh_Record_Component.SCALAR].load_chunk()
series.flush()  

print("T_ss shape:", T_ss.shape, "| total inventory [g]:", float(np.sum(T_ss)))
print("D_ss shape:", D_ss.shape, "| total inventory [g]:", float(np.sum(D_ss)))


# --------------------------------------------------
# Time-series inventories
# --------------------------------------------------
block_info = []
for k in iter_keys:
    it = series.iterations[k]
    if ("T" in it.meshes) and ("D" in it.meshes):
        start_step = int(it.get_attribute("block_start_step")) if it.contains_attribute("block_start_step") else k * (block_len or 0)
        block_info.append((start_step, k))
block_info.sort(key=lambda x: x[0])
blocks = [k for _, k in block_info]
starts = [s for s, _ in block_info] 
print(f"Reading {len(blocks)} blocks: {blocks[0]} .. {blocks[-1]}")

# Read blocks and build timeseries
times = []
T_total = []
D_total = []
T_subsys = []   
D_subsys = []

for start_step, k in zip(starts, blocks):
    it = series.iterations[k]
    this_len = int(it.get_attribute("block_len"))

    T_block = it.meshes["T"][io.Mesh_Record_Component.SCALAR].load_chunk()
    D_block = it.meshes["D"][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()

    if downsample > 1:
        idx = np.arange(0, this_len, downsample)
        T_block = T_block[:, idx]
        D_block = D_block[:, idx]
    else:
        idx = np.arange(this_len)

    tvec = (start_step + idx).astype(np.float64)
    if dt is not None:
        tvec *= dt

    times.append(tvec)
    T_total.append(np.sum(T_block, axis=0))
    D_total.append(np.sum(D_block, axis=0))
    T_subsys.append(T_block)
    D_subsys.append(D_block)

# Concatenate
times = np.concatenate(times)
T_total = np.concatenate(T_total)
D_total = np.concatenate(D_total)
T_subsys = np.concatenate(T_subsys, axis=1)
D_subsys = np.concatenate(D_subsys, axis=1)

print("\n=== Time-series summary ===")
if len(times):
    print("Time range:", float(times[0]), "to", float(times[-1]), f"[{dt_unit or 'unknown'}]")


# --------------------------------------------------
# Plot inventories
# --------------------------------------------------

# Total inventories
xlabel = f"Time [{dt_unit}]" if dt is not None else "Step"

plt.figure()
plt.plot(times, T_total, label="Tritium")
plt.plot(times, D_total, label="Deuterium")
plt.xlabel(xlabel)
plt.ylabel("Total Inventory [g]")
plt.title("Total inventory of Tritium and Deuterium across the subsystems")
plt.grid(True)
plt.legend()

png_path = os.path.join(OUTDIR, f"rhino_total_blocks_down{downsample}_ADIOS.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", png_path)


# Inventories per subsystem
plt.figure()
for i in range(T_subsys.shape[0]):
    plt.plot(times, T_subsys[i, :], label=labels_T[i])
plt.xlabel(xlabel)
plt.ylabel("Inventory [g]")
plt.title("Tritium inventory by subsystem")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)
png_path = os.path.join(OUTDIR, f"rhino_T_subsystems_down{downsample}_ADIOS.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", png_path)

plt.figure()
for i in range(D_subsys.shape[0]):
    plt.plot(times, D_subsys[i, :], label=labels_D[i])
plt.xlabel(xlabel)
plt.ylabel("Inventory [g]")
plt.title("Deuterium inventory by subsystem")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)
png_path = os.path.join(OUTDIR, f"rhino_D_subsystems_down{downsample}_ADIOS.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", png_path)