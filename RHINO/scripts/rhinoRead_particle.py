import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
import os


# Configuration
BP5_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhinoADIOS_particles.bp5"
OUTDIR = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/plots"
os.makedirs(OUTDIR, exist_ok=True)

downsample = 100

series = io.Series(
    BP5_PATH,
    io.Access.read_only,
    '{"verify_homogeneous_extents": false}'
)

print("\n=== File opened ===")
print("Path:", BP5_PATH)
print("Series attributes:", list(series.attributes))


# Metadata
def get_attr(name, default=None):
    return series.get_attribute(name) if series.contains_attribute(name) else default

dt = get_attr("timeStep", None)
dt_unit = get_attr("timeUnit", None)
Nt_meta = get_attr("timeSteps", None)

labels_subsystem = get_attr("subsystemLabels", None)
if labels_subsystem is None:
    raise RuntimeError("Missing 'subsystemLabels' attribute.")

labels_subsystem = [str(x) for x in labels_subsystem]

print("\n=== File-level metadata ===")
print("dt:", dt)
print("dt_unit:", dt_unit)
print("timeSteps:", Nt_meta)
print("n_subsystems:", len(labels_subsystem))


# Load Particle Data 
it = series.iterations[0]

def load_species(name):
    if name not in it.particles:
        raise KeyError(f"Particle species '{name}' not found.")

    pt = it.particles[name]

    inv = pt["inventory"][io.Record_Component.SCALAR].load_chunk()
    ss  = pt["inventory_steady"][io.Record_Component.SCALAR].load_chunk()

    series.flush()

    return np.asarray(inv), np.asarray(ss)

T_ts, T_ss = load_species("Tritium")
D_ts, D_ss = load_species("Deuterium")

print("\n=== Data shapes ===")
print("T_ts:", T_ts.shape)
print("D_ts:", D_ts.shape)
print("T_ss:", T_ss.shape)
print("D_ss:", D_ss.shape)


# Sanity Checks
if T_ts.ndim != 2 or D_ts.ndim != 2:
    raise ValueError("Expected inventory to be 2D (subsystem × time).")

if T_ss.ndim != 1 or D_ss.ndim != 1:
    raise ValueError("Expected steady-state to be 1D (subsystem).")

nsub, Nt = T_ts.shape

if D_ts.shape != (nsub, Nt):
    raise ValueError("Tritium and Deuterium shapes do not match.")

if T_ss.shape[0] != nsub or D_ss.shape[0] != nsub:
    raise ValueError("Steady-state length mismatch.")


# Build Time Vector
if downsample is None or downsample < 1:
    downsample = 1

idx = np.arange(0, Nt, downsample)

if dt is not None:
    times = idx.astype(np.float64) * float(dt)
    xlabel = f"Time [{dt_unit}]"
else:
    times = idx.astype(np.float64)
    xlabel = "Step"

T_ds = T_ts[:, idx]
D_ds = D_ts[:, idx]


# Total inventories
T_total = np.sum(T_ds, axis=0)
D_total = np.sum(D_ds, axis=0)

print("\n=== Time-series summary ===")
if len(times):
    print("Time range:", float(times[0]), "to", float(times[-1]), f"[{dt_unit}]")

# Steady-State Check
last_T = T_ts[:, -1]
last_D = D_ts[:, -1]

print("\n=== Steady-state consistency ===")
print("T exact match:", np.allclose(last_T, T_ss, atol=0, rtol=0),
      "| max abs diff:", float(np.max(np.abs(last_T - T_ss))))
print("D exact match:", np.allclose(last_D, D_ss, atol=0, rtol=0),
      "| max abs diff:", float(np.max(np.abs(last_D - D_ss))))


# Plot total inventory
plt.figure()
plt.plot(times, T_total, label="Tritium")
plt.plot(times, D_total, label="Deuterium")

plt.xlabel(xlabel)
plt.ylabel("Total Inventory [g]")
plt.title("Total inventory across subsystems")
plt.grid(True)
plt.legend()

png_path = os.path.join(OUTDIR, f"rhino_particles_total.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()

print("Saved:", png_path)


# Plot inventory per subsystem
plt.figure()
for i, name in enumerate(labels_subsystem):
    plt.plot(times, T_ds[i, :], label=name)

plt.xlabel(xlabel)
plt.ylabel("Inventory [g]")
plt.title("Tritium inventory by subsystem")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)

png_path = os.path.join(OUTDIR, f"rhino_particles_T_subsystems.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()

print("Saved:", png_path)

plt.figure()
for i, name in enumerate(labels_subsystem):
    plt.plot(times, D_ds[i, :], label=name)

plt.xlabel(xlabel)
plt.ylabel("Inventory [g]")
plt.title("Deuterium inventory by subsystem (particle)")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)

png_path = os.path.join(OUTDIR, f"rhino_particles_D_subsystems.png")
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()

print("Saved:", png_path)

series.close()
