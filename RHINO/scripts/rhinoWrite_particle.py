# This file reads RHINO output data (pkl files) and writes it in an openPMD series
# Uses ADIOS2 as a backend (.bp5) with variable-based iteration encoding
# Stores data in particle records (i.e. DataFrame-like)
# Data is stored in a shared folder on NERSC
# /global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino

import openpmd_api as io
import numpy as np
import pandas as pd

# Configuration
OUTPUT_PATH = "/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/bp_output/rhino_particles.bp5"
DATA_PATH = "/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/older_data/Data"
# OUTPUT_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/bp_output/rhino_particles.bp5"
# DATA_PATH = "/home/ccb/ccb/Projects/IFE_AmSC/RHINO/Data"
SECONDS_PER_DAY = 86400.0


# Load RHINO Data
T_ts_df = pd.read_pickle(f"{DATA_PATH}/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T.pkl")
D_ts_df = pd.read_pickle(f"{DATA_PATH}/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D.pkl")
T_ss_df = pd.read_pickle(f"{DATA_PATH}/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_T_SteadyState.pkl")
D_ss_df = pd.read_pickle(f"{DATA_PATH}/2025-12-16/06-57-07_AmSC_Generic_FuelCycle_D_SteadyState.pkl")
meta_df = pd.read_pickle(f"{DATA_PATH}/2025-12-16/06-57-07_AmSC_meta.pkl")

# metadata
meta = meta_df[0].to_dict()
for k, v in list(meta.items()):
    if isinstance(v, np.generic):
        meta[k] = v.item()

labels_subsystem = [
    "Storage_Delivery", 
    "Fueling", 
    "Fusion_Chamber_Pump", 
    "Pd_Cleanup",
    "Protium_Removal", 
    "Exhaust_Processing", 
    "Gas_Detrit", 
    "Water_Detrit",
    "Glovebox", 
    "Stack", 
    "Isotope_Seperation", 
    "Blanket_Extraction",
    "Heat_Exchanger", 
    "Power_Conv_Loop", 
    "Vent_Detrit", 
    "Blanket",
    "Decay_Box", 
    "Stack_Box", 
    "Burn_Box", 
    "Gen_Box", 
    "Uptake_Box"
]
canon = {name: i for i, name in enumerate(labels_subsystem)}

T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64))
nSubsystems, Nt = T_ts.shape
D_ts = np.zeros((nSubsystems, Nt), dtype=np.float64)
for name, row in D_ts_df.iterrows():
    D_ts[canon[str(name)], :] = row.to_numpy(dtype=np.float64)
T_ss = T_ss_df.iloc[:, 0].to_numpy(dtype=np.float64)
D_ss = np.zeros((nSubsystems,), dtype=np.float64)
for name, row in D_ss_df.iterrows():
    D_ss[canon[str(name)]] = float(row.iloc[0])

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

# ==============================
# Root Schema Metadata
series.set_attribute("schema", "OpenPMD+X")
series.set_attribute("schemaVersion", "0.17.0")
series.set_attribute("basePath", "/data")
series.set_attribute("particlesPath", "inventory")
series.set_attribute("iterationEncoding", "variableBased")
series.set_attribute("iterationFormat", "bp5")


# ==============================
# Software Metadata
series.set_attribute("software/softwareName", "RHINO")
series.set_attribute("software/softwareDescription", "RHINO: Fusion Pilot Plant fuel cycle simulation")
series.set_attribute("software/softwareVersion", "1.0")
# Placeholder for additional attributes under software/
series.set_attribute("software/versionControlSoftware", "")
series.set_attribute("software/softwareCommit", "")
series.set_attribute("software/softwareDocumentation", "")


# ==============================
# Provenance metadata
series.set_attribute("provenance/author", "Holly Flynn")
series.set_attribute("provenance/authorAffiliation", "Savannah River National Laboratory")
series.set_attribute("provenance/authorEmail", "Holly.Flynn@srnl.doe.gov")
series.set_attribute("provenance/creationDate", "2025-12-16")
# Placeholder for additional attributes under provenance/
series.set_attribute("provenance/inputDirectory", "")
series.set_attribute("provenance/inputFiles", "")
series.set_attribute("provenance/originalDataDirectory", "")
series.set_attribute("provenance/originalDataFiles", "")


# ==============================
# System metadata (placeholders)
series.set_attribute("system/systemIP", "")
series.set_attribute("system/systemDescription", "")


# ==============================
# Application metadata (RHINO-specific)
series.set_attribute("metadata/subsystems/description", "Subsystems of the pilot plant fuel cycle")
ids = np.arange(nSubsystems, dtype=np.uint64)
series.set_attribute("metadata/subsystems/id", ids.tolist())
series.set_attribute("metadata/subsystems/labels", labels_subsystem)
# series.set_attribute("metadata/subsystems/mapping", canon)        # dict as an attribute is not supported
# Placeholder for additional attributes under metadata/subsystems/
series.set_attribute("metadata/subsystems/connections", "")
series.set_attribute("metadata/subsystems/connectionType", "")

series.set_attribute("metadata/species/description", "Gas species in the fuel cycle")
series.set_attribute("metadata/species/names", ["Tritium", "Deuterium"])
# Placeholder for additional attributes under metadata/species/
series.set_attribute("metadata/species/Tritium/subsystems", "")
series.set_attribute("metadata/species/Tritium/subsystemsConnection", "")
series.set_attribute("metadata/species/Deuterium/subsystems", "")
series.set_attribute("metadata/species/Deuterium/subsystemsConnection", "")

series.flush()


# ==============================
# Data
# ==============================
# Iteration
dt = float(meta["dt"]) if "dt" in meta else None
series.particles_path = "inventory"
it = series.snapshots()[0]
if dt is not None:
    it.time = 0.0
    it.dt = float(dt)
    it.time_unit_SI = SECONDS_PER_DAY
    it.set_attribute("timeUnitLabel", "day")

pos_data = np.arange(nSubsystems, dtype=np.float64)

def write_species(name, data_ts, data_ss):

    pt = it.particles[name]

    pt.set_attribute("description", "Inventory across subsystems for species " + name)
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

    # time-series inventory data (copy to ensure writable for store_chunk)
    data_arr = np.ascontiguousarray(data_ts).copy()
    inv_rec = pt["mass"][io.Record_Component.SCALAR]
    inv_rec.reset_dataset(io.Dataset(data_arr.dtype, data_arr.shape))
    inv_rec.store_chunk(data_arr)

    pt["mass"].unit_dimension = {io.Unit_Dimension.M: 1}
    inv_rec.unit_SI = 1e-3

    # steady-state inventory data (copy to ensure writable for store_chunk)
    ss_arr = np.ascontiguousarray(data_ss).copy()
    ss_rec = pt["mass_steady"][io.Record_Component.SCALAR]
    ss_rec.reset_dataset(io.Dataset(ss_arr.dtype, ss_arr.shape))
    ss_rec.store_chunk(ss_arr)

    pt["mass_steady"].unit_dimension = {io.Unit_Dimension.M: 1}
    ss_rec.unit_SI = 1e-3

write_species("Tritium", T_ts, T_ss)
write_species("Deuterium", D_ts, D_ss)

it.close()
series.close()

print("RHINO data written to ADIOS-OpenPMD in particle representation.")
print("Output:", OUTPUT_PATH)
