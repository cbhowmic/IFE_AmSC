"""
This file reads RHINO data (pkl files) 
and writes it into an openPMD series
using ADIOS2 as a backend (.bp5)
This file only reads T ss data.
"""


"""
get_estimate(
tritium_mass_per_target_mg : float,
pulse_rep_rate_Hz : float,
burn_fraction: float [=beta in input file],
thermal_power_output_MW : float [double check with Holly, suspect it's a constant =500, should be in the intput file but it's not],
blanket_flow_rate_kg_s : float [],
) ->
T_startup_inventory_g : float [Holly's script],
T_doubling_time_hrs : float [Holly's script],
T_steady_state_hrs : float [max across subs of time to eq. do it in the write script],
T_steady_state_isotope_separation_inventory_g : float [ss at isotope box, write script]

"""

import sys
import pickle 
import numpy as np
import pandas as pd
import openpmd_api as io
 
    # Notes for Holly: 
    # It would be nice if the folder structure were DATA_PATH/PREFIX/T.pkl, etc.
    # It would be nice to have the individual input files in each folder 
    # The metafile name does not follow the same convention as the others (i.e. no infix)
    # Attributes = Metadata
    # Note that the subsystems ids are not the same for T and D: correct?

def save_bp(DATA_PATH, PREFIX, INFIX, OUTPUT_PATH):
    #####################
    ### Configuration ###
    #####################
    # Folder where the bp file will be saved 
    
    
    
    #PREFIX="22-21-28" 
    #INFIX ="IFE_AmSC_500MW_FuelCycle"  
    
    # Paths where RHINO data is 
    #RHINO_PATH = "/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data" 
    #DATA_PATH  = f"{RHINO_PATH}/Power&BurnFractionScan_Daily_Reduced1" 
    INPUT_PATH = f"{DATA_PATH}/{PREFIX}_IFE_input.pkl" 
    POSTPROC_PATH = f"{DATA_PATH}/{PREFIX}_IFE_processed.pkl" 
    
    # Import input file 
    #sys.path.append(RHINO_PATH)
    #from makeJSON import InputFile
    
    #######################
    ### Load RHINO data ###
    #######################
    # Time-Series and Steady-State data
    T_ts_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_T_reduced.pkl")
    #D_ts_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_D.pkl")
    T_ss_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_T_SteadyState.pkl")
    # 0 processing time
    # 1 ss
    # 2 flow 
    #D_ss_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_{INFIX}_D_SteadyState.pkl")
    # Metafile
    meta_df = pd.read_pickle(f"{DATA_PATH}/{PREFIX}_IFE_meta.pkl")
    # Input file
    InputFile = pd.read_pickle(INPUT_PATH)
    PostProcData = pd.read_pickle(POSTPROC_PATH)
    
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
    #times = np.arange(0, endtime+dt, dt)
    #nt = len(times)
    
    # Useful constant 
    SECONDS_PER_DAY = 86400.0
    
    ######################
    ### Extract inputs ###
    ######################
    my_inputs = {}
    # Tritium 
    my_inputs["Tritium"] = {}
    for k,v in InputFile['Systems_T'].items():
        if isinstance(v, list):
            my_inputs["Tritium"][v[0]] = {"id": int(k), 
                                          "processing time": v[1], 
                                          "nonradioactive loss fraction": v[2], 
                                          "fractional inflows": v[3], 
                                          "initial mass": v[4], 
                                          "source": v[5], 
                                          "injectors": v[6], 
                                          "label": v[7]}
    
    # Superset of all subsystem names 
    all_subsystems_names = list(my_inputs["Tritium"])
    
    ###########################
    ### Extract inventories ###
    ###########################
    my_inventory = {}
    # Tritium 
    T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64))
    nSubsystemsT, Nt = T_ts.shape
    T_ss = T_ss_df.iloc[:, 1].to_numpy(dtype=np.float64)
    my_inventory["Tritium"]   = {"data_ss": T_ss, "data_ts": T_ts}
    
    # Save times
    times = np.linspace(0, endtime, Nt)

    # Time to steady state 
    # if "box" not in substystem and subs != Storage_delivery -> compute time to ss
<<<<<<< HEAD
    SS_time = ...
=======
    ss_tol = 0.02  # tolerance band for ss

    def time_to_steady_state(times, ts, ss, tol=0.01):
        if np.isclose(ss, 0.0):
            err = np.abs(ts - ss)
        else:
            err = np.abs((ts - ss) / ss)
        inside = err <= tol
        i = len(times) - 1
        while i >= 0 and inside[i]:
            i -= 1
        if i == len(times) - 1:
            return np.nan
        return float(times[i + 1])
    t_ss_by_subsystem = {}
    for i, name in enumerate(all_subsystems_names):
        lname = name.lower()
        if "box" in lname or "storage" in lname or "delivery" in lname:
            continue
        t_ss_by_subsystem[name] = time_to_steady_state(
            times=times,
            ts=T_ts[i, :],
            ss=T_ss[i],
            tol=ss_tol,
        )

    valid_t_ss = [t for t in t_ss_by_subsystem.values() if not np.isnan(t)]
    if len(valid_t_ss) > 0:
        ss_time = float(np.max(valid_t_ss))
    else:
        ss_time = np.nan
>>>>>>> baeb9b7 (Computing time to  steady-state)

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

    # Post-processes outputs computed by Holly 
    # I0 (g)	Imin (g)	I_startup (g)	I_subtract (g)	reserve_time (days)	Iops (g)	plant_doubling_time (days)
    for k,v in PostProcData.items():
        series.set_attribute(f"output:{k}", v)
<<<<<<< HEAD
    series.set_attribute("output:Steady state time (days)", SS_time)
=======
    series.set_attribute("output:Steady state time (days)", ss_time)
>>>>>>> baeb9b7 (Computing time to  steady-state)
    
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
    def write_species(name, data_ss, data_ts):
    
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
    
    write_species("Tritium", T_ss, T_ts)
    
    ######################
    ### Close and save ###
    ######################
    it.close()
    series.close()
    print("RHINO data written to ADIOS-OpenPMD in particle representation.")
    print("Output:", OUTPUT_PATH)