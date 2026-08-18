# RHINO AI-Ready Data Workflow

This workflow turns RHINO fuel-cycle simulation results into standardized,
searchable, and machine-learning-ready data. The overall idea is to preserve
the scientific content and provenance of each RHINO run while progressively
adding the structure needed to find runs, compare them, and assemble training
datasets.

The workflow separates format conversion, campaign organization, and feature
construction into three layers. Each layer has a well-defined output that can
be inspected or reused independently.

## Workflow at a Glance

```text
RHINO simulation outputs
        |
        v
1. Shim: standardize each run as openPMD/ADIOS2 BP5
        |
        v
2. Campaign: package archival replicas, group runs, and build an index
        |
        v
3. Feature extraction: select and join inputs, outputs, and array data
        |
        v
ML-ready feature table
        |
        v
Surrogate training and downstream analysis
```

## Workflow Layers

### 1. Shim and Data Standardization

The shim reads RHINO pickle outputs and converts each simulation run to an
openPMD-compliant ADIOS2 BP5 series. It records run inputs, derived outputs,
subsystem metadata, steady-state inventories, time-series inventories, and
provenance in a consistent representation.

**Input:** Native RHINO output files (`.pkl`)  
**Output:** One BP5 series (`.bp5`) per run

See the [shim layer guide](1_shim/README.md) for conversion commands and usage.

### 2. Campaign Management

The campaign layer packages standardized BP5 runs into checksummed,
TAR-indexed archival replicas, groups the runs into campaign archives, and
builds a SQLite-backed campaign index. The index makes run metadata and
attributes searchable with reusable SQL queries, so downstream work can select
scientifically relevant runs without scanning every source file.

**Input:** Standardized RHINO BP5 series  
**Output:** TAR replicas and checksums (`.tar`, `.sha256`, `.tar.idx`), campaign
archives (`.aca`), and a campaign index (`.acx`)

See the [campaign layer guide](2_campaign/README.md) for configuration,
archive creation, indexing, and queries.

### 3. Feature Extraction

The feature-extraction layer turns selected campaign data into a tabular
dataset. It queries scalar inputs and targets from the campaign index, reads
array-valued scientific data from the archived ADIOS datasets when needed, and
joins those values by RHINO run. Feature specifications keep the choice of
model inputs and outputs explicit.

**Input:** Campaign index, campaign archives, SQL queries, and feature
specifications  
**Output:** An ML-ready feature table, such as `rhino_features.csv`

See the [feature-extraction layer guide](3_feature_extraction/README.md) for
feature definitions and dataset construction.

## Why the Layers Are Separate

- The shim provides a common scientific data format independent of later ML
  choices.
- Campaign management adds checksummed archival replicas, scalable
  organization, discovery, and selection.
- Feature extraction can evolve for different models or analyses without
  rerunning RHINO or rewriting the ingestion workflow.

Together, these layers create a traceable path from a source simulation run to
the exact values used in a machine-learning dataset. The resulting feature
table is consumed by the separate
[surrogate-training workflow](../ML/surrogate_training/README.md).
