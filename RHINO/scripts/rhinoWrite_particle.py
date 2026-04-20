"""
This file reads RHINO data (pkl files) 
and writes it into an openPMD series
using ADIOS2 as a backend (.bp5)
"""

import sys
import numpy as np
import pandas as pd
import openpmd_api as io
 
# Notes for Holly: 
# It would be nice if the folder structure were DATA_PATH/PREFIX/T.pkl, etc.
# It would be nice to have the individual input files in each folder 
# The metafile name does not follow the same convention as the others (i.e. no infix)
# Attributes = Metadata
# Note that the subsystems ids are not the same for T and D: correct?

#####################
### Configuration ###
#####################
# Folder where the bp file will be saved 
OUTPUT_PATH = "/global/cfs/cdirs/m3239/aforment/IFE_AmSC/RHINO/scripts/tmp/rhino_particles.bp5"

# Paths where RHINO data is 
RHINO_PATH = "/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/more_data_runs" 
DATA_PATH  = f"{RHINO_PATH}/Extraction/Data" 
INPUT_PATH = f"{RHINO_PATH}/makeJSON.py" 
PREFIX="12-34-31" 
INFIX ="IFE_AmSC_500MW_FuelCycle"  

# Import input file 
sys.path.append(RHINO_PATH)
from makeJSON import InputFile

#######################
### Load RHINO data ###
#######################
# Time-Series and Steady-State data
T_ts_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_T.pkl")
D_ts_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_D.pkl")
T_ss_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_T_SteadyState.pkl")
D_ss_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_D_SteadyState.pkl")
# Metafile
meta_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_IFE_meta.pkl")

########################
### Extract metadata ###
########################
# Convert metadata to dictionary
meta = meta_df[0].to_dict()
for k, v in list(meta.items()):
    if isinstance(v, np.generic):
        meta[k] = v.item()

# Get timestep
dt = float(meta["dt"]) if "dt" in meta else None
if dt is None:
    raise ValueError("dt is required for time-series data")   

# Get simulation time 
endtime = float(meta["calc_length"])
times = np.arange(0, endtime+dt, dt)
nt = len(times)

# Useful constant 
SECONDS_PER_DAY = 86400.0

######################
### Extract inputs ###
######################
my_inputs = {}
# Tritium 
my_inputs["Tritium"] = {}
for k,v in InputFile['Systems_T'].items():
    my_inputs["Tritium"][v[0]] = {"id": k, 
                                  "processing time": v[1], 
                                  "nonradioactive loss fraction": v[2], 
                                  "fractional inflows": v[3], 
                                  "initial mass": v[4], 
                                  "source": v[5], 
                                  "injectors": v[6], 
                                  "label": v[7]}
# Deuterium 
my_inputs["Deuterium"] = {}
for k,v in InputFile['Systems_D'].items():
    my_inputs["Deuterium"][v[0]] = {"id": k, 
                                  "processing time": v[1], 
                                  "nonradioactive loss fraction": v[2], 
                                  "fractional inflows": v[3], 
                                  "initial mass": v[4], 
                                  "source": v[5], 
                                  "injectors": v[6], 
                                  "label": v[7]}
# Superset of all subsystem names 
all_subsystems_names = np.union1d(list(my_inputs["Tritium"]), list(my_inputs["Deuterium"]))

###########################
### Extract inventories ###
###########################
my_inventory = {}
# Tritium 
T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64))
nSubsystemsT, Nt = T_ts.shape
T_ss = T_ss_df.iloc[:, 1].to_numpy(dtype=np.float64)
my_inventory["Tritium"]   = {"data_ts": T_ts, "data_ss": T_ss}
# Deuterium 
D_ts = np.ascontiguousarray(D_ts_df.to_numpy(dtype=np.float64))
nSubsystemsD, Nt = D_ts.shape
D_ss = D_ss_df.iloc[:, 1].to_numpy(dtype=np.float64)
my_inventory["Deuterium"] = {"data_ts": D_ts, "data_ss": D_ss}

#############################
### Create openPMD series ###
#############################
adios2_cfg = r'''
{
  "iteration_encoding": "variable_based",
  "adios2": {
    "modifiable_attributes": false,
    "use_group_table": false,
    "engine": {
      "type": "bp5",
      "parameters": {
        "StatsLevel": "1",
        "AsyncWrite": "guided"
      }
    }
  }
}
'''
series = io.Series(OUTPUT_PATH, io.Access_Type.create_linear, adios2_cfg)
print("Converting RHINO data into openPMD/ADIOS2 format...")
print(f"Input: {DATA_PATH}")

#########################
### Series attributes ###
#########################
# openPMD-ready attributes
series.particles_path = "inventory"
series.set_attribute("software", "RHINO")
series.set_attribute("softwareVersion", "1.0")
series.set_attribute("softwareDescription", "RHINO: Fusion Pilot Plant fuel cycle simulation")
series.set_attribute("author", "Holly Flynn")
series.set_attribute("authorAffiliation", "Savannah River National Laboratory")
series.set_attribute("authorEmail", "Holly.Flynn@srnl.doe.gov")
series.set_attribute("date", "2025-12-16")
series.set_attribute("machine", "")
series.set_attribute("comment", f"Provenance: data path is {DATA_PATH}, input file is {INPUT_PATH}")
# General inputs (common to D and T) 
series.set_attribute("input:TBR:Tritium Breeding Ratio", InputFile["System Inputs"]["TBR"])
series.set_attribute("input:TBRr:Required Tritium Breeding Ratio", InputFile["System Inputs"]["TBRr"])
series.set_attribute("input:beta:Burn fraction", InputFile["System Inputs"]["beta"])
series.set_attribute("input:eta:Fueling efficiency", InputFile["System Inputs"]["eta"])
series.set_attribute("input:Ndotminus:Tritium burned per day", InputFile["System Inputs"]["Ndotminus"])
series.set_attribute("input:MW:Power output in MW", InputFile["System Inputs"]["MW"])
series.set_attribute("input:I0_SD:Starting inventory", InputFile["System Inputs"]["I0_SD"])

##########################
### Create iteration 0 ###
##########################
it = series.snapshots()[0]
it.time = 0.0
it.dt = float(dt)
it.time_unit_SI = SECONDS_PER_DAY

#######################
### Save time array ###
#######################
species = it.particles["Times"]  
species.set_attribute("description", "Times")
record = species["data"]
record.unit_dimension =  {io.Unit_Dimension.T: 1}
record.unit_SI = SECONDS_PER_DAY
data = np.ascontiguousarray(times).copy()
dataset = io.Dataset(times.dtype, times.shape)
component = record[io.Record_Component.SCALAR]
component = record[io.Record_Component.SCALAR]
component.reset_dataset(dataset)
component.store_chunk(data)

#############################
### Save species: T and D ###
#############################
def write_species(name, data_ts, data_ss):

    pt = it.particles[name]
    pt.set_attribute("description", "Inventory across subsystems for species " + name)
    pt.set_attribute("timeAxis", 1)  
    pt.set_attribute("subsystemsAxis", 0)

    # subsystems data
    record = pt["subsystems"]
    for k,v in my_inputs[name].items():
        component = record[k]
        component.make_empty
        for kk, vv in my_inputs[name][k].items():
            component.set_attribute(kk, vv)
        component.reset_dataset(io.Dataset(data.dtype, data.shape))
        component.store_chunk(data)

    # time-series inventory data
    data_arr = np.ascontiguousarray(data_ts).copy()
    inv_rec = pt["mass"][io.Record_Component.SCALAR]
    inv_rec.reset_dataset(io.Dataset(data_arr.dtype, data_arr.shape))
    inv_rec.store_chunk(data_arr)

    pt["mass"].unit_dimension = {io.Unit_Dimension.M: 1}
    inv_rec.unit_SI = 1e-3
    
    # steady-state inventory data
    ss_arr = np.ascontiguousarray(data_ss).copy()
    ss_rec = pt["mass_steady"][io.Record_Component.SCALAR]
    ss_rec.reset_dataset(io.Dataset(ss_arr.dtype, ss_arr.shape))
    ss_rec.store_chunk(ss_arr)

    pt["mass_steady"].unit_dimension = {io.Unit_Dimension.M: 1}
    ss_rec.unit_SI = 1e-3

write_species("Tritium", T_ts, T_ss)
write_species("Deuterium", D_ts, D_ss)

######################
### Close and save ###
######################
it.close()
series.close()
print("RHINO data written to ADIOS-OpenPMD in particle representation.")
print("Output:", OUTPUT_PATH)
