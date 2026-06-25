from pathlib import Path
import sys
import numpy as np
import pandas as pd
import openpmd_api as io
from rhinoWrite import rhino_to_adios

ROOT_PATH=Path("/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data/")
SCENARIOS=["2026-04-29", "2026-04-30", "2026-05-01"]

OUTPUT_ROOT=Path("/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/surrogate_bp_output/")


for scenario in SCENARIOS:
    print(scenario)
    try:
        scenario_path = ROOT_PATH / scenario
        data_root = scenario_path

        if not data_root.exists():
            print(f"Scenario folder missing or no Data/: {data_root} — skipping")
            continue
        
        for tfile in sorted(data_root.glob("*_T_reduced.pkl")):
            try:
                # Extract correct prefixes
                run_prefix_data = tfile.name.split("_T_reduced.pkl")[0]
                run_prefix_meta = "_".join(run_prefix_data.split("_")[:2])
                run_time_prefix = run_prefix_data.split("_")[0]
                
                prefix = run_time_prefix 
                infix = "_".join(run_prefix_data.split("_")[1:])
                #print(f"Found run: {infix}")
                if scenario==SCENARIOS[1] and prefix=="03-51-04":
                    continue
                    
            except Exception as e:
                print(f"ERROR processing run '{run_prefix_data}' in scenario '{scenario}': {e}")
                continue

            # build output path 
            safe_param = scenario.replace(" ", "_")
            safe_param = safe_param.replace("&", "And")
            
            out_name = f"{safe_param}/{run_time_prefix}.bp5"
            OUTPUT_PATH = OUTPUT_ROOT / out_name
        
            rhino_to_adios(scenario_path, prefix, infix, OUTPUT_PATH)
    
    except Exception as e:
        print(f"ERROR processing scenario '{scenario}': {e}")
        continue