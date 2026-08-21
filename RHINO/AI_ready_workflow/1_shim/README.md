# RHINO Shim Layer

This layer converts raw RHINO simulation outputs into
openPMD-compliant ADIOS2 BP5 series.

## Responsibilities

- Standardize RHINO outputs
- Ingest metadata from the RHINO simulation
- Store run attributes, inputs/outputs
- Capture scientific context for later analysis, visualization, troubleshooting

## Outputs

One `ADIOS (.bp5)` data product per RHINO run.

The outputs of the shim layer are then consumed by the campaign management layer.

## Files and Usage

### `rhinoWrite.py`

Contains the conversion logic for one RHINO run. It reads the run's pickle
files, extracts simulation inputs, derived outputs, subsystem metadata,
steady-state inventories, and time-series inventories, and writes them as an
openPMD/ADIOS2 BP5 series. The `rhino-write` command invokes this file's
command-line entry point.

Convert one run:

```bash
PREFIX=14-24-04
INFIX=IFE_AmSC_500MW_FuelCycle
DATA_PATH="/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data/2026-04-29/"
OUTPUT_PATH="$HOME/my.bp5"

rhino-write \
    --data-path "$DATA_PATH" \
    --prefix "$PREFIX" \
    --infix "$INFIX" \
    --output-path "$OUTPUT_PATH"
```

The prefix and infix identify the related RHINO files within `DATA_PATH`. The
result is a single BP5 data product at `OUTPUT_PATH`.

The converter also derives the simulation datetime from the source naming
convention. `DATA_PATH` must end with a scenario date in `YYYY-MM-DD` form, and
`PREFIX` must be the run time in `HH-MM-SS` form. For example, scenario
`2026-04-29` and prefix `14-24-04` are stored in the series-level `date`
attribute as `2026-04-29T14:24:04`. The source data does not encode a timezone,
so the stored ISO 8601 datetime has no timezone offset. Invalid names stop the
conversion instead of silently writing incorrect provenance.

### `rhinoWrite_multiple.py`

Provides batch conversion for one or more scenario directories. It discovers
files ending in `_T_reduced.pkl`, derives each run's prefix and infix from its
filename, and calls the single-run converter. The `rhino-write-multiple`
command invokes this file's command-line entry point.

Convert every discovered run in one or more scenarios:

```bash
OUTPUT_ROOT="$HOME/rhino_bp5"
ROOT_PATH="/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/Surrogate Data/"

rhino-write-multiple \
    --root-path "$ROOT_PATH" \
    --scenarios 2026-04-29 2026-04-30 \
    --output-root "$OUTPUT_ROOT"
```

To omit a known run, add `--skip-run SCENARIO:PREFIX`. The option can be
specified more than once:

```bash
rhino-write-multiple \
    --root-path "$ROOT_PATH" \
    --scenarios 2026-04-29 \
    --output-root "$OUTPUT_ROOT" \
    --skip-run 2026-04-29:14-24-04
```

If a discovered run is incomplete or cannot be converted, the batch command
reports that run's filename, scenario, and error, then continues with the next
run. Review the command output for `ERROR processing run` before treating a
batch conversion as complete.

### `examples/rhinoWrite_multiple.ipynb`

An interactive example of the multi-run conversion workflow. Use it to inspect
the source data, adjust scenario selections, and run batch conversion from a
Jupyter environment.

```bash
jupyter lab examples/rhinoWrite_multiple.ipynb
```

### `schema.md`

Documents the schema considered by the shim,
including software and provenance metadata, subsystem descriptions, and data
records. It is a reference document and is not executed.
