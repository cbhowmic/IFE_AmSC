import openpmd_api as io
import json
import numpy as np
import matplotlib.pyplot as plt
import os

BP5_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhinoADIOS.bp5"
OUTDIR = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/plots"
os.makedirs(OUTDIR, exist_ok=True)

downsample = 100

series = io.Series(BP5_PATH, io.Access.read_only)

print("\n=== File opened ===")
print("Path:", BP5_PATH)

# --------------------------------------------------
# Metadata 
# --------------------------------------------------

print("Series attributes:", list(series.attributes))
print("Iter0 attributes:", list(series.iterations[0].attributes))


# File-level metadata
def get_attr(path, default=None):
    return series.get_attribute(path) if series.contains_attribute(path) else default

dt = get_attr("time step", None)
dt_unit = get_attr("time unit", None)
t_start = get_attr("start time", None)
t_end = get_attr("end time", None)

labels_subsystem = get_attr("subsystem labels", None)
if labels_subsystem is not None:
    labels_subsystem = [str(x) for x in labels_subsystem]

print("\n=== File-level metadata ===")
print("dt:", dt)
print("dt_unit:", dt_unit)
print("t_start:", t_start)
print("t_end:", t_end)

canon = {name: i for i, name in enumerate(labels_subsystem)} if labels_subsystem else {}
labels_T_rows = json.loads(get_attr("Tritium subsystems", "[]"))
labels_D_rows = json.loads(get_attr("Deuterium subsystems", "[]"))
map_T = json.loads(get_attr("Tritium to all subsystems_index", "{}"))
map_D = json.loads(get_attr("Deuterium to all subsystems_index", "{}"))
print("n_canonical_subsystems:", len(labels_subsystem) if labels_subsystem else None)
print("n_T_subsystems_provided:", len(labels_T_rows))
print("n_D_subsystems_provided:", len(labels_D_rows))


it = series.iterations[0]

def load_scalar_mesh(mesh_name):
    if mesh_name not in it.meshes:
        raise KeyError(f"Mesh '{mesh_name}' not found. Available meshes: {list(it.meshes.keys())}")
    arr = it.meshes[mesh_name][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()
    return np.asarray(arr)

# Load data meshes
T_ts = load_scalar_mesh("Tritium")
D_ts = load_scalar_mesh("Deuterium")
T_ss = load_scalar_mesh("Tritium_steady")
D_ss = load_scalar_mesh("Deuterium_steady")

print("\n=== Data shapes ===")
print("T_ts:", T_ts.shape)
print("D_ts:", D_ts.shape)
print("T_ss:", T_ss.shape)
print("D_ss:", D_ss.shape)


# -------------------Basic check------------------------------
if T_ts.ndim != 2 or D_ts.ndim != 2:
    raise ValueError("Expected Tritium/Deuterium time-series meshes to be 2D (nSubsystem, Nt).")
if T_ss.ndim != 1 or D_ss.ndim != 1:
    raise ValueError("Expected Tritium_steady/Deuterium_steady to be 1D (nSubsystem,).")
nsub, Nt = T_ts.shape
if D_ts.shape != (nsub, Nt):
    raise ValueError("Tritium and Deuterium time-series shapes do not match.")
if T_ss.shape[0] != nsub or D_ss.shape[0] != nsub:
    raise ValueError("Steady-state length does not match nSubsystem.")

# ---------------------------------------------------
#  Build time vectors and downsample
# ---------------------------------------------------
if downsample is None or downsample < 1:
    downsample = 1
idx = np.arange(0, Nt, downsample)

if dt is not None:
    times = idx.astype(np.float64) * float(dt)
    xlabel = f"Time [{dt_unit or 'unknown'}]"
else:
    times = idx.astype(np.float64)
    xlabel = "Step"

T_ds = T_ts[:, idx]
D_ds = D_ts[:, idx]

def rows_from_mapping(labels_rows, mapping, fallback_canon):
    rows = []
    names = []
    for name in labels_rows:
        if name in mapping:
            rows.append(int(mapping[name]))
            names.append(name)
        elif name in fallback_canon:
            rows.append(int(fallback_canon[name]))
            names.append(name)
        else:
            # silently skip unknown labels (or raise if you prefer)
            print(f"Warning: subsystem label '{name}' not found in mapping/canonical list; skipping.")
    return np.array(rows, dtype=int), names

rows_T, names_T = rows_from_mapping(labels_T_rows, map_T, canon)
rows_D, names_D = rows_from_mapping(labels_D_rows, map_D, canon)

T_sel = T_ds[rows_T, :] if len(rows_T) else np.empty((0, len(times)))
D_sel = D_ds[rows_D, :] if len(rows_D) else np.empty((0, len(times)))

T_total = np.sum(T_sel, axis=0) if T_sel.size else np.zeros_like(times)
D_total = np.sum(D_sel, axis=0) if D_sel.size else np.zeros_like(times)

print("\n=== Time-series summary ===")
if len(times):
    print("Time range:", float(times[0]), "to", float(times[-1]), f"[{dt_unit or 'unknown'}]")


# -------------------Steady-state values vs last entries in time-series-------------------------------
if len(rows_T):
    last_T = T_ts[rows_T, -1]
    T_ss_sel = T_ss[rows_T]
    print("\n=== SS check (Tritium)===")
    print("Exact match:", np.allclose(last_T, T_ss_sel, atol=0, rtol=0),
          "| max abs diff:", float(np.max(np.abs(last_T - T_ss_sel))))
if len(rows_D):
    last_D = D_ts[rows_D, -1]
    D_ss_sel = D_ss[rows_D]
    print("\n=== SS check (Deuterium)===")
    print("Exact match:", np.allclose(last_D, D_ss_sel, atol=0, rtol=0),
          "| max abs diff:", float(np.max(np.abs(last_D - D_ss_sel))))

# --------------------------------------------------
# Plot inventories
# --------------------------------------------------
xlabel = f"Time [{dt_unit}]" if dt is not None else "Step"

# Total inventories
plt.figure()
plt.plot(times, T_total, label=f"Tritium (n={len(rows_T)})")
plt.plot(times, D_total, label=f"Deuterium (n={len(rows_D)})")
plt.xlabel(xlabel)
plt.ylabel("Total Inventory [g]")
plt.title("Total inventory of Tritium and Deuterium across subsystems")
plt.grid(True)
plt.legend()
png_path = os.path.join(OUTDIR, f"rhino_total_down{downsample}_ADIOS.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", png_path)

# Tritium by subsystem
if len(rows_T):
    plt.figure()
    for j, name in enumerate(names_T):
        plt.plot(times, T_sel[j, :], label=name)
    plt.xlabel(xlabel)
    plt.ylabel("Inventory [g]")
    plt.title("Tritium inventory by subsystem")
    plt.grid(True)
    plt.legend(fontsize=7, ncol=2)

    png_path = os.path.join(OUTDIR, f"rhino_T_subsystems_down{downsample}_ADIOS.png")
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.show()
    print("Saved:", png_path)
else:
    print("No Tritium subsystem labels found in metadata; skipping T-by-subsystem plot.")


# Deuterium by subsystem
if len(rows_D):
    plt.figure()
    for j, name in enumerate(names_D):
        plt.plot(times, D_sel[j, :], label=name)
    plt.xlabel(xlabel)
    plt.ylabel("Inventory [g]")
    plt.title("Deuterium inventory by subsystem")
    plt.grid(True)
    plt.legend(fontsize=7, ncol=2)

    png_path = os.path.join(OUTDIR, f"rhino_D_subsystems_down{downsample}_ADIOS.png")
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.show()
    print("Saved:", png_path)
else:
    print("No Deuterium subsystem labels found in metadata; skipping D-by-subsystem plot.")

series.close()