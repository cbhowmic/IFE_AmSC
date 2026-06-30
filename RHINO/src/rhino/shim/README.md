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