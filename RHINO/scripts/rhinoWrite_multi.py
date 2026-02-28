# This file reads iterates over multiple RHINO output data (pkl files)
# Writes the data using openPMD series and ADIOS2 as a backend (.bp5) with variable-based iteration encoding
# Stores data in particle records (i.e. DataFrame-like)
# Each run outputs a separate ADIOS directory, all the outputs are stored in output/ directory with the run ID as the identifier
# RHINO data is stored in a shared folder on NERSC
# /global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino

import openpmd_api as io
import numpy as np
import pandas as pd
import re
from pathlib import Path
from datetime import datetime, timezone
import socket
import json

ROOT_PATH = Path("/home/ccb/ccb/Dropbox/IFE_AmSC/RHINO/Data")
SCENARIOS = ["Burn Fraction Changes", "TBR", "Protium Removal", "Extraction", "Power Scan"]    
OUTPUT_ROOT = Path("/home/ccb/ccb/Projects/IFE_AmSC/RHINO/bp_output")
SECONDS_PER_DAY = 86400.0

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


# ADIOS2 configuration (variable-based iteration encoding, bp5 engine)
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

# canonical subsystem labels 
labels_subsystem = [
    "Storage_Delivery", "Fueling", "Fusion_Chamber_Pump", "Pd_Cleanup",
    "Protium_Removal", "Exhaust_Processing", "Gas_Detrit", "Water_Detrit",
    "Glovebox", "Stack", "Isotope_Seperation", "Blanket_Extraction",
    "Heat_Exchanger", "Power_Conv_Loop", "Vent_Detrit", "Blanket",
    "Decay_Box", "Stack_Box", "Burn_Box", "Gen_Box", "Uptake_Box"
]
canon = {name: i for i, name in enumerate(labels_subsystem)}


# helper function to parse structured fields from log file
def parse_log_text(log_text):
    md = {}

    m = re.search(r"RHINO(?:\.py)?:\s*(.+)", log_text)
    if m:
        md["softwareDescription"] = m.group(1).strip()

    m = re.search(r"Author:\s*(.+)", log_text)
    if m:
        md["author"] = m.group(1).strip()

    m = re.search(r"Contributors:\s*(.+)", log_text)
    if m:
        md["contributors"] = [c.strip() for c in m.group(1).split(",") if c.strip()]

    m = re.search(r"Status:\s*(.+)", log_text)
    if m:
        md["status"] = m.group(1).strip()

    m = re.search(r"Version\s*:? (\d+)", log_text)
    if m:
        md["softwareVersion"] = m.group(1).strip()
    else:
        m2 = re.search(r"Version\s*(\S+)", log_text)
        if m2:
            md["softwareVersion"] = m2.group(1).strip()

    m = re.search(r"Latest version as of:\s*([0-9/ -]+)", log_text)
    if m:
        md["releaseDate"] = m.group(1).strip()

    m = re.search(r"Date of current run:\s*([0-9\-\/]+)", log_text)
    if m:
        md["runDate"] = m.group(1).strip()

    m = re.search(r"Time of current run:\s*([0-9:\-]+)", log_text)
    if m:
        md["runTime"] = m.group(1).strip()

    m = re.search(r"from '(.*?)'", log_text)
    if m:
        md["executablePath"] = m.group(1).strip()

    m = re.search(r"\(?([A-Za-z]:?[/\\].*?\.json)\)?", log_text)
    if m:
        md["inputFilePath"] = m.group(1).strip()
        md["inputFiles"] = [Path(md["inputFilePath"]).name]
    else:
        m2 = re.search(r"(\S+\.json)", log_text)
        if m2:
            md["inputFiles"] = [m2.group(1)]
    return md



# =====================================
# Main loop: iterate over parameter folders and runs, write each run to a separate ADIOS2 series
# =====================================
for scenario in SCENARIOS:
    try:
        scenario_path = ROOT_PATH / scenario
        data_root = scenario_path / "Data"
        log_root = scenario_path / "Logs"

        if not data_root.exists():
            print(f"Scenario folder missing or no Data/: {data_root} — skipping")
            continue

        if not log_root.exists():
            print(f"Scenario folder missing Logs/: {log_root} — continuing (logs may be missing)")
        print(f"\nProcessing scenario: '{scenario}' at {scenario_path}")
        
        # iterate runs under Data/
        for tfile in sorted(data_root.glob("*_FuelCycle_T.pkl")):
            
            # Extract correct prefixes
            run_prefix_data = tfile.name.split("_FuelCycle_T.pkl")[0]
            run_prefix_meta = "_".join(run_prefix_data.split("_")[:2])
            run_time_prefix = run_prefix_data.split("_")[0]
            print(f"Found run: {run_prefix_data}")

            t_path = data_root / f"{run_prefix_data}_FuelCycle_T.pkl"
            d_path = data_root / f"{run_prefix_data}_FuelCycle_D.pkl"
            t_ss_path = data_root / f"{run_prefix_data}_FuelCycle_T_SteadyState.pkl"
            d_ss_path = data_root / f"{run_prefix_data}_FuelCycle_D_SteadyState.pkl"
            meta_path = data_root / f"{run_prefix_meta}_meta.pkl"

            missing = [p for p in (t_path, d_path, t_ss_path, d_ss_path, meta_path) if not p.exists()]
            if missing:
                print(f"WARNING: missing data files for run {run_prefix_data}")
                print("Missing:", [m.name for m in missing])
                continue
            
            # find log file 
            log_file_candidate = log_root / run_time_prefix
            log_text = ""
            parsed_log = {}

            if log_file_candidate.exists():
                with open(log_file_candidate, "r", encoding="utf-8", errors="replace") as f:
                    log_text = f.read()

                version_match = re.search(r"Version (\d+)", log_text)
                date_match = re.search(r"Date of current run:\s*(.*)", log_text)
                time_match = re.search(r"Time of current run:\s*(.*)", log_text)
                author_match = re.search(r"Author:\s*(.*)", log_text)
                contrib_match = re.search(r"Contributors:\s*(.*)", log_text)

                if version_match:
                    parsed_log["softwareVersion"] = version_match.group(1)

                if date_match:
                    parsed_log["runDate"] = date_match.group(1)

                if time_match:
                    parsed_log["runTime"] = time_match.group(1)

                if author_match:
                    parsed_log["author"] = author_match.group(1)

                if contrib_match:
                    parsed_log["contributors"] = [
                        c.strip() for c in contrib_match.group(1).split(",")
                    ]
            else:
                print(f"WARNING: log file not found for {run_time_prefix}")

            # --- load RHINO data from .pkl files ---
            T_ts_df = pd.read_pickle(t_path)
            D_ts_df = pd.read_pickle(d_path)
            T_ss_df = pd.read_pickle(t_ss_path)
            D_ss_df = pd.read_pickle(d_ss_path)
            meta_df = pd.read_pickle(meta_path)

            T_subsystems_present = [str(x) for x in list(T_ts_df.index)]  
            D_subsystems_present = [str(x) for x in list(D_ts_df.index)] 
            
            # convert metadata
            meta = meta_df[0].to_dict() if hasattr(meta_df, "__len__") else dict(meta_df)
            for k, v in list(meta.items()):
                if isinstance(v, np.generic):
                    meta[k] = v.item()

            # canonical arrays
            T_ts = np.ascontiguousarray(T_ts_df.to_numpy(dtype=np.float64))
            nSubsystems, Nt = T_ts.shape
            D_ts = np.zeros((nSubsystems, Nt), dtype=np.float64)
            for name, row in D_ts_df.iterrows():
                D_ts[canon[str(name)], :] = row.to_numpy(dtype=np.float64)
            T_ss = T_ss_df.iloc[:, 0].to_numpy(dtype=np.float64)
            D_ss = np.zeros((nSubsystems,), dtype=np.float64)
            for name, row in D_ss_df.iterrows():
                D_ss[canon[str(name)]] = float(row.iloc[0])

            # build output path 
            safe_param = scenario.replace(" ", "_")
            out_name = f"{run_prefix_data}.bp5"
            OUTPUT_PATH = OUTPUT_ROOT / out_name

            # create OpenPMD series (one per run)
            series = io.Series(str(OUTPUT_PATH), io.Access_Type.create_linear, adios2_cfg)
            
            # root metadata
            series.set_attribute("schema", "OpenPMD+X")
            series.set_attribute("schemaVersion", "0.17.0")
            series.set_attribute("basePath", "/data")
            series.set_attribute("particlesPath", "inventory")
            series.set_attribute("iterationEncoding", "variableBased")
            
            # software metadata (some from parsed log)
            series.set_attribute("software/softwareName", "RHINO")            
            if "softwareDescription" in parsed_log:
                series.set_attribute("software/softwareDescription", parsed_log["softwareDescription"])
            else:
                series.set_attribute("software/softwareDescription", "Reduced Hydrogen INventory Optimization model for Fusion Fuel Cycle")
            if "softwareVersion" in parsed_log:
                series.set_attribute("software/softwareVersion", str(parsed_log["softwareVersion"]))
            else:
                series.set_attribute("software/softwareVersion", str(meta.get("version", "unknown")))
            if "status" in parsed_log:
                series.set_attribute("software/status", parsed_log["status"])
            if "releaseDate" in parsed_log:
                series.set_attribute("software/releaseDate", parsed_log["releaseDate"])
            series.set_attribute("software/versionControlSoftware", "")
            series.set_attribute("software/softwareCommit", "")
            series.set_attribute("software/softwareDocumentation", "")
            
            # provenance metadata
            series.set_attribute("provenance/author", parsed_log.get("author", "Holly B. Flynn"))
            if "contributors" in parsed_log:
                series.set_attribute("provenance/contributors", parsed_log["contributors"])
            series.set_attribute("provenance/authorAffiliation", parsed_log.get("authorAffiliation", "Savannah River National Laboratory"))
            series.set_attribute("provenance/authorEmail", "Holly.Flynn@srnl.doe.gov")
            if "runDate" in parsed_log:
                series.set_attribute("provenance/creationDate", parsed_log["runDate"])
            else:
                series.set_attribute("provenance/creationDate", run_prefix_data.split("_")[0])
            if "runTime" in parsed_log:
                series.set_attribute("provenance/creationTime", parsed_log["runTime"])
            else:
                series.set_attribute("provenance/creationTime", run_time_prefix)
            if "executablePath" in parsed_log:
                series.set_attribute("provenance/executablePath", parsed_log["executablePath"])
            if "inputFiles" in parsed_log:
                series.set_attribute("provenance/inputFiles", parsed_log["inputFiles"])
            elif "inputFilePath" in parsed_log:
                series.set_attribute("provenance/inputFiles", [Path(parsed_log["inputFilePath"]).name])
            series.set_attribute("provenance/inputDirectory", "")
            series.set_attribute("provenance/originalDataDirectory", str(data_root))
            series.set_attribute("provenance/originalDataFiles", [p.name for p in (t_path, d_path, t_ss_path, d_ss_path, meta_path)])

            # system metadata
            try:
                ip = socket.gethostbyname(socket.gethostname())
            except Exception:
                    ip = ""
            series.set_attribute("system/systemIP", ip)
            series.set_attribute("system/systemDescription", f"Scenario: {scenario}; run: {run_prefix_data}")
        
            # metadata/subsystems
            ids = np.arange(nSubsystems, dtype=np.uint64)
            series.set_attribute("metadata/subsystems/description", "Subsystems of the pilot plant fuel cycle")
            series.set_attribute("metadata/subsystems/id", ids.tolist())
            series.set_attribute("metadata/subsystems/labels", labels_subsystem)
            series.set_attribute("metadata/subsystems/connections", [])
            series.set_attribute("metadata/subsystems/connectionType", "directed_nonweighted")
            
            # metadata/species
            series.set_attribute("metadata/species/description", "Gas species in the fuel cycle")
            series.set_attribute("metadata/species/names", ["Tritium", "Deuterium"])
            series.set_attribute("metadata/species/Tritium/subsystems", T_subsystems_present)
            series.set_attribute("metadata/species/Deuterium/subsystems", D_subsystems_present)
            series.set_attribute("metadata/species/Tritium/subsystemsConnection", "")
            series.set_attribute("metadata/species/Deuterium/subsystemsConnection", "")

            # additional code-specific metadata
            series.set_attribute("scenario", scenario)
            series.set_attribute("runID", run_prefix_data)
            
            series.flush()

            series.particles_path = "inventory"
            it = series.snapshots()[0]

            # time info 
            dt_val = float(meta["dt"]) if "dt" in meta else None
            if dt_val is not None:
                it.time = 0.0
                it.dt = float(dt_val)
                it.time_unit_SI = SECONDS_PER_DAY
                it.set_attribute("timeUnitLabel", "day")
                it.set_attribute("timeSteps", int(Nt))            

            pos_data = np.arange(nSubsystems, dtype=np.float64)

            # helper function to write each species
            def write_species(it_obj, ids_arr, pos_arr, name, data_ts, data_ss):
                pt = it_obj.particles[name]
                pt.set_attribute("description", f"Inventory across subsystems for species {name}")
                pt.set_attribute("particleShape", "point")
                pt.set_attribute("particleType", "subsystem")
                
                id_rec = pt["id"][io.Record_Component.SCALAR]
                id_rec.reset_dataset(io.Dataset(ids_arr.dtype, ids_arr.shape))
                id_rec.store_chunk(ids_arr)
                
                pos_rec = pt["position"]["x"]
                pos_rec.reset_dataset(io.Dataset(pos_arr.dtype, pos_arr.shape))
                pos_rec.store_chunk(pos_arr)
                pos_rec.unit_SI = 1.0
                
                # time-series inventory
                data_arr = np.ascontiguousarray(data_ts).copy()
                inv_rec = pt["mass"][io.Record_Component.SCALAR]
                inv_rec.reset_dataset(io.Dataset(data_arr.dtype, data_arr.shape))
                inv_rec.store_chunk(data_arr)
                pt["mass"].unit_dimension = {io.Unit_Dimension.M: 1}
                inv_rec.unit_SI = 1e-3
                
                # steady-state inventory
                ss_arr = np.ascontiguousarray(data_ss).copy()
                ss_rec = pt["mass_steady"][io.Record_Component.SCALAR]
                ss_rec.reset_dataset(io.Dataset(ss_arr.dtype, ss_arr.shape))
                ss_rec.store_chunk(ss_arr)
                pt["mass_steady"].unit_dimension = {io.Unit_Dimension.M: 1}
                ss_rec.unit_SI = 1e-3
                pt["mass"].set_attribute("timeAxis", 1)

            write_species(it, ids, pos_data, "Tritium", T_ts, T_ss)
            write_species(it, ids, pos_data, "Deuterium", D_ts, D_ss)

            it.close()
            series.close()
            print(f"Wrote: {OUTPUT_PATH}")

    except Exception as e:
        print(f"ERROR processing scenario '{scenario}': {e}")
        continue

print("\nAll parameter studies processed.")