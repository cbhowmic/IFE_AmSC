import openpmd_api as io
import numpy as np
import pandas as pd

# Configuration
OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/output/rhinoADIOS_particles.bp5"
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

labels_subsystem = [
    "Storage_Delivery", "Fueling", "Fusion_Chamber_Pump", "Pd_Cleanup",
    "Protium_Removal", "Exhaust_Processing", "Gas_Detrit", "Water_Detrit",
    "Glovebox", "Stack", "Isotope_Seperation", "Blanket_Extraction",
    "Heat_Exchanger", "Power_Conv_Loop", "Vent_Detrit", "Blanket",
    "Decay_Box", "Stack_Box", "Burn_Box", "Gen_Box", "Uptake_Box"
]

canon = {name: i for i, name in enumerate(labels_subsystem)}

# Build canonical arrays
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

# Create Series
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

print("Writing RHINO data in particle representation...")

# Series metadata
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
    series.set_attribute("timeSteps", int(Nt))

series.flush()

# Iteration
series.particles_path = "inventory"
it = series.iterations[0]

if dt is not None:
    it.time = 0.0
    it.dt = float(dt)
    it.time_unit_SI = SECONDS_PER_DAY

ids = np.arange(nSubsystems, dtype=np.uint64)
pos_data = np.arange(nSubsystems, dtype=np.float64)

def write_species(name, data_ts, data_ss):

    pt = it.particles[name]

    pt.set_attribute("particleShape", "point")
    pt.set_attribute("particleType", "subsystem")
    pt.set_attribute("timeAxis", 1)  

    # id 
    id_rec = pt["id"][io.Record_Component.SCALAR]
    id_rec.reset_dataset(io.Dataset(ids.dtype, ids.shape))
    id_rec.store_chunk(ids)

    # logical position
    pos_rec = pt["position"]["x"]
    pos_rec.reset_dataset(io.Dataset(pos_data.dtype, pos_data.shape))
    pos_rec.store_chunk(pos_data)
    pos_rec.unit_SI = 1.0

    # time-series inventory data
    data_arr = np.ascontiguousarray(data_ts)
    inv_rec = pt["inventory"][io.Record_Component.SCALAR]
    inv_rec.reset_dataset(io.Dataset(data_arr.dtype, data_arr.shape))
    inv_rec.store_chunk(data_arr)

    pt["inventory"].unit_dimension = {io.Unit_Dimension.M: 1}
    inv_rec.unit_SI = 1e-3

    # steady-state inventory data
    ss_arr = np.ascontiguousarray(data_ss)
    ss_rec = pt["inventory_steady"][io.Record_Component.SCALAR]
    ss_rec.reset_dataset(io.Dataset(ss_arr.dtype, ss_arr.shape))
    ss_rec.store_chunk(ss_arr)

    pt["inventory_steady"].unit_dimension = {io.Unit_Dimension.M: 1}
    ss_rec.unit_SI = 1e-3

write_species("Tritium", T_ts, T_ss)
write_species("Deuterium", D_ts, D_ss)

it.close()
series.close()

print("FAST particle file written.")
print("Output:", OUTPUT_PATH)
