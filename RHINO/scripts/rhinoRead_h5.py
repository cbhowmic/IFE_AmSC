import openpmd_api as io
import json
import numpy as np
import matplotlib.pyplot as plt
import os
import math

H5_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhino_block.h5"
OUTDIR = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/plots"
os.makedirs(OUTDIR, exist_ok=True)

series = io.Series(H5_PATH, io.Access.read_only)
downsample = 100

print("\n=== File opened ===")
print("Path:", H5_PATH)
print("openPMD version:", series.openPMD)

# --------------------------------------------------
# Metadata 
# --------------------------------------------------

# File-level metadata
block_len = int(series.get_attribute("block_len")) if series.contains_attribute("block_len") else None
Nt_total = int(series.get_attribute("Nt_total")) if series.contains_attribute("Nt_total") else None
dt = float(series.get_attribute("meta/dt")) if series.contains_attribute("meta/dt") else None
dt_unit = series.get_attribute("meta/dt_unit") if series.contains_attribute("meta/dt_unit") else None  

print("\n=== File-level metadata ===")
print("block_len:", block_len)
print("Nt_total:", Nt_total)
print("dt:", dt)
print("dt_unit:", dt_unit)


# Labels 
labels_T = json.loads(series.get_attribute("labels/T_rows")) if series.contains_attribute("labels/T_rows") else None
labels_D = json.loads(series.get_attribute("labels/D_rows")) if series.contains_attribute("labels/D_rows") else None
print("\n=== Subsystem labels ===")
if labels_T is not None:
    print("T labels:", len(labels_T), labels_T)
if labels_D is not None:
    print("D labels:", len(labels_D), labels_D)


# Series-level attributes
def get_attr(obj, name):
    try:
        return obj.get_attribute(name)
    except Exception:
        return None

print("\n=== Series-level metadata ===")
for key in [
    "schema",
    "schema/version",
    "application",
    "run/generated_utc",
    "meta/dt",
    "meta/dt_unit",
    "meta/units",
]:
    val = get_attr(series, key)
    if val is not None:
        print(f"{key}: {val}")



# --------------------------------------------------
# Steady-state inventories 
# --------------------------------------------------
iter_keys = sorted(series.iterations)
steady_it = iter_keys[-1]
it = series.iterations[steady_it]
if ("T_steady" not in it.meshes) and ("D_steady" not in it.meshes):
    steady_it = None

T_ss = None
D_ss = None
if steady_it is None:
    print("\n=== Steady-state ===")
    print("No steady-state iteration found (no T_steady/D_steady meshes).")
else:
    it = series.iterations[steady_it]
    print("\n=== Steady-state ===")
    print("Iteration index:", steady_it)

    if "T_steady" in it.meshes:
        T_ss = it.meshes["T_steady"][io.Mesh_Record_Component.SCALAR].load_chunk()
    if "D_steady" in it.meshes:
        D_ss = it.meshes["D_steady"][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()

    if T_ss is not None:
        print("T_ss shape:", T_ss.shape, "| total [g]:", float(np.sum(T_ss)))
    if D_ss is not None:
        print("D_ss shape:", D_ss.shape, "| total [g]:", float(np.sum(D_ss)))



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

for start_step, b in zip(starts, blocks):
    it = series.iterations[b]
    this_len = int(it.get_attribute("block_len")) if it.contains_attribute("block_len") else (block_len or 0)

    T_block = it.meshes["T"][io.Mesh_Record_Component.SCALAR].load_chunk()
    D_block = it.meshes["D"][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()

    if downsample > 1:
        T_block = T_block[:, ::downsample]
        D_block = D_block[:, ::downsample]
        t_index = np.arange(0, this_len, downsample)
    else:
        t_index = np.arange(this_len)

    if dt is not None:
        tvec = (start_step + t_index) * dt
    else:
        tvec = (start_step + t_index).astype(np.float64)

    # Compute total inventories across subsystems
    T_sum = np.sum(T_block, axis=0)
    D_sum = np.sum(D_block, axis=0)

    times.append(tvec)
    T_total.append(T_sum)
    D_total.append(D_sum)
    T_subsys.append(T_block)
    D_subsys.append(D_block)

# Concatenate
times = np.concatenate(times) if len(times) else np.array([], dtype=np.float64)
T_total = np.concatenate(T_total) if len(T_total) else np.array([], dtype=np.float64)
D_total = np.concatenate(D_total) if len(D_total) else np.array([], dtype=np.float64)
T_subsys = np.concatenate(T_subsys, axis=1) if len(T_subsys) else None
D_subsys = np.concatenate(D_subsys, axis=1) if len(D_subsys) else None

print("\n=== Time-series summary ===")
if len(times):
    print("Time range:", float(times[0]), "to", float(times[-1]), f"[{dt_unit or 'unknown'}]")


# --------------------------------------------------
# Plot inventories
# --------------------------------------------------

# Total inventories
plt.figure()
plt.plot(times, T_total, label="Tritium")
plt.plot(times, D_total, label="Deuterium")
plt.xlabel(f"Time [{dt_unit}]" if dt is not None and dt_unit is not None else ("Time [day]" if dt is not None else "Step"))
plt.ylabel("Total Inventory [g]")
plt.title("Total inventory of Tritium and Deuterium across the subsystems")
plt.grid(True)
plt.legend()

png_path = os.path.join(OUTDIR, f"rhino_total_blocks_down{downsample}.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", png_path)


# Inventories per subsystem
if T_subsys is not None and labels_T is not None:
    plt.figure()
    for i in range(T_subsys.shape[0]):
        lbl = labels_T[i] if i < len(labels_T) else f"T_row_{i}"
        plt.plot(times, T_subsys[i, :], label=lbl)
    plt.xlabel(f"Time [{dt_unit}]" if dt is not None and dt_unit is not None else ("Time [day]" if dt is not None else "Step"))
    plt.ylabel("Inventory [g]")
    plt.title("Tritium inventory by subsystem")
    plt.grid(True)
    plt.legend(fontsize=7, ncol=2)

    png_path = os.path.join(OUTDIR, f"rhino_T_subsystems_down{downsample}.png")
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.show()
    print("Saved:", png_path)

if D_subsys is not None and labels_D is not None:
    plt.figure()
    for i in range(D_subsys.shape[0]):
        lbl = labels_D[i] if i < len(labels_D) else f"D_row_{i}"
        plt.plot(times, D_subsys[i, :], label=lbl)
    plt.xlabel(f"Time [{dt_unit}]" if dt is not None and dt_unit is not None else ("Time [day]" if dt is not None else "Step"))
    plt.ylabel("Inventory [g]")
    plt.title("Deuterium inventory by subsystem")
    plt.grid(True)
    plt.legend(fontsize=7, ncol=2)

    png_path = os.path.join(OUTDIR, f"rhino_D_subsystems_down{downsample}.png")
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.show()
    print("Saved:", png_path)