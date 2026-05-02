# RHINO Module

This module processes RHINO simulation outputs and converts them into AI-ready ADIOS datasets following the project schema.

This repository contains Python scripts to:
1) **Write** RHINO fuel-cycle inventory into **OpenPMD** format using **ADIOS2 BP5**, and  
2) **Read**, validate, and **plot** time-series and steady-state inventories from the produced `.bp5` file.

The scripts are designed for RHINO inventory outputs stored as Pandas pickles (`.pkl`) and follow a fixed canonical list of plant subsystems. A given run may only provide a subset of those subsystems; the writer maps provided subsystems into the canonical ordering, and the reader plots only the subsystems included in the run.

---

## Directory Structure

```text
IFE_AmSC/
├── RHINO/
    ├── scripts
        ├── rhinoWriteADIOS.py        # Write RHINO inventory data to ADIOS2 
        ├── rhinoReadADIOS.py         # Read, validate, and plot RHINO inventory data
    ├── README.md
    ├── Data/                         # User places RHINO data here
    │   └── ...
    ├── output/                       # ADIOS2 BP5 files written here
    └── plots/                        # Generated plots saved here
```

## 1) Obtain RHINO data
RHINO simulation data are available from the following SharePoint location:

https://doellnl.sharepoint.com/teams/FESAmSC-IFEPPDesignDig.Twin/Shared%20Documents/Forms/AllItems.aspx?id=%2Fteams%2FFESAmSC%2DIFEPPDesignDig%2ETwin%2FShared%20Documents%2FGeneral%2FData&viewid=5d1dddd2%2D1bbe%2D44ff%2D818c%2D39dd2133c599

The directory contains data from multiple RHINO runs. Each individual run should include the following five required files:

### Time-series inventories
- `*_Generic_FuelCycle_T.pkl`
- `*_Generic_FuelCycle_D.pkl`

These files contain time-dependent inventories for Tritium (T) and Deuterium (D) across RHINO plant subsystems.

### Steady-state inventories
- `*_Generic_FuelCycle_T_SteadyState.pkl`
- `*_Generic_FuelCycle_D_SteadyState.pkl`

These files represent the steady-state (final) inventories corresponding to the time-series data.

### Metadata
- `*_meta.pkl`

This file contains run-level metadata describing the RHINO simulation configuration.


## 2) Place RHINO data
Place the downloaded data inside `RHINO/Data/`.

Each RHINO run should be placed in its own subdirectory containing the five required files.


## 3) Environment Setup
This module uses the project-wide Conda environment:
```bash
conda activate IFE_AmSC
```
If the environment has not yet been created, build it from the repository root:
```bash
conda env create -f environment.yml
conda activate IFE_AmSC
```

## 4) Running the Scripts

### 4.1) Write ADIOS file from RHINO data
From inside `RHINO/scripts/`:
```bash
python rhinoWriteADIOS.py
```
This script:
- Reads `.pkl` inventory files from `RHINO/Data/`
- Maps subsystems into canonical ordering
- Writes an OpenPMD-compliant ADIOS2 BP5 file
- Stores the output in `RHINO/output/`.


### 4.2) Read, Validate, and Plot
From inside `RHINO/scripts/`:
```bash
python rhinoReadADIOS.py
```
This script:
- Reads the generated BP5 file from RHINO/output/
- Validates data consistency
- Generates time-series and steady-state plots
- Saves plots into `RHINO/plots/`.
