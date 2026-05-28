# RHINO Campaign Management Layer

This layer organizes openPMD/ADIOS BP5 files produced by the shim layer
into campaign archives and campaign indexes.

## Responsibilities

- Create campaign archives from BP5 outputs
- Add RHINO runs to campaign archives
- Build or refresh campaign indexes
- Query archived runs using SQL
- Expose reusable query functions for feature extraction and surrogate workflows

## Inputs

- `.bp5` openPMD/ADIOS series from `RHINO/shim`

## Outputs

- Campaign archive, e.g. `rhino1.aca`
- Campaign index database
- Queryable run metadata and attributes