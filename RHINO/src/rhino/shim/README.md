# RHINO Shim Layer

This layer converts raw RHINO simulation outputs into
openPMD-compliant ADIOS2 BP5 series.

## Responsibilities

- Standardize RHINO outputs
- Write openPMD metadata
- Store run attributes
- Store subsystem time-series
- Produce AI-ready campaign data

## Outputs

One `.bp5` series per RHINO run.

These outputs are consumed by:
- campaign management
- campaign indexing
- feature extraction
- surrogate training

## Usage
```
PREFIX=14-24-04
INFIX=IFE_AmSC_500MW_FuelCycle
DATA_PATH="/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data/2026-04-29/"
OUTPUT_PATH=$HOME/my.bp5
rhino-write --data-path "$DATA_PATH" --prefix $PREFIX --infix $INFIX --output-path $OUTPUT_PATH
```

```
OUTPUT_ROOT=$HOME
ROOT_PATH="/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data/"
SCENARIOS=2026-04-29
rhino-write-multiple --root-path "$ROOT_PATH" --scenarios $SCENARIOS --output-root $OUTPUT_ROOT 
```
