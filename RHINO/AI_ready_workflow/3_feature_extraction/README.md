# RHINO Feature Extraction Layer

This layer builds an ML-ready table from RHINO runs registered with HPC
Campaign. Scalar attributes are selected from the SQLite-backed campaign index,
while array-backed values are read from the referenced ADIOS datasets only when
the feature specification requests them.

Feature selection is declared in JSON. Python is responsible for validating the
specification, executing its SQL and ADIOS readers, joining rows by campaign
run, and writing the resulting CSV.

## Workflow

```text
Campaign index (.acx) + completed campaign archives (.aca)
                   + feature_spec.json
                   + reusable SQL queries
    |
    v
Select run metadata, scalar inputs, and outputs from the campaign index
    |
    v
Read requested array values through the ADIOS campaign reader
    |
    v
Outer-join values by archive, dataset ID, and run ID
    |
    v
ML-ready feature table (.csv)
    |
    v
Surrogate training and downstream analysis
```

The campaign store may resolve live BP5 replicas, TAR-backed replicas, or other
locations recorded in the `.aca`; the feature layer does not encode those
storage details.

## Files

### `feature_spec.json`

Declares run metadata, model inputs, and outputs without embedding workflow
choices in Python. Each entry has a stable output `key`, a `role`, and a
supported `source`.

The current specification produces these columns:

| Role | Output column | Source |
| --- | --- | --- |
| Metadata | `simulation_date` | Series `date` attribute via Campaign SQL |
| Metadata | `simulation_time` | Series `date` attribute via Campaign SQL |
| Metadata | `simulation_datetime` | Series `date` attribute via Campaign SQL |
| Input | `tritium_burning_rate` | Campaign SQL attribute |
| Input | `burn_fraction` | Campaign SQL attribute |
| Output | `plant_doubling_time_days` | Campaign SQL attribute |
| Output | `minimum_startup_inventory_g` | Campaign SQL attribute |
| Output | `tritium_in_isotope_separation` | ADIOS array element |

### `spec_loader.py`

Loads and validates the JSON document. The optional `metadata` list remains
separate from required model `inputs` and `outputs`. The loader rejects
unsupported schema versions, duplicate keys, incorrect roles, unsupported
sources, missing source-specific fields, and invalid ADIOS indexes before
reading campaign data.

### `campaign_reader.py`

Executes SQL-backed feature queries, checks their result columns, rejects
duplicate run rows, and joins all features using campaign run identity:

```text
archive + datasetid + run_id
```

It then groups runs by `.aca` and requests any ADIOS-backed features.

### `campaign_adios.py`

Reads array-backed variables from a campaign archive with `adios2.FileReader`.
Campaign variables are addressed as:

```text
<run_id>/<variable path>
```

For the current tritium target, `index: 10` selects the isotope-separation
subsystem from the steady-state mass array. The accompanying `subsystem` value
documents the scientific meaning of that index and is included in error
messages.

### `build_features.py`

Provides the command-line workflow, loads the JSON specification, validates
input paths, builds the table, reports missing values, and writes CSV output.
It has no user-specific campaign paths embedded in the source.

### Campaign SQL queries

Reusable queries live in the sibling campaign layer:

```text
../2_campaign/queries/
```

Every SQL-backed feature query must return:

```text
archive, datasetid, run_id, <configured feature column>
```

The supplied feature queries accept the named parameter `:archive_name`, which
is populated from the command-line archive pattern.

## JSON Specification

The top-level document contains a schema version and separate metadata, input,
and output lists:

```json
{
  "schema_version": 1,
  "metadata": [
    {
      "key": "simulation_datetime",
      "role": "metadata",
      "source": "campaign_sql",
      "query": "list_simulation_datetime.sql",
      "column": "simulation_datetime"
    }
  ],
  "inputs": [
    {
      "key": "burn_fraction",
      "role": "input",
      "source": "campaign_sql",
      "query": "list_burn_fraction.sql",
      "column": "burn_fraction"
    }
  ],
  "outputs": [
    {
      "key": "tritium_in_isotope_separation",
      "role": "output",
      "source": "campaign_adios",
      "variable": "/data/inventory/Tritium/mass_steady",
      "subsystem": "Isotope_Seperation",
      "index": 10
    }
  ]
}
```

Supported sources are:

- `campaign_sql`: requires `query` and `column`.
- `campaign_adios`: requires `variable` and optionally accepts a
  zero-based `index`.

The keys in this JSON become the final table column names. Metadata columns
describe each run but are not model inputs. Input and output keys should remain
aligned with the surrogate-training interface.

## Simulation Datetime Provenance

The shim writes a series-level `/date` attribute by combining the source
scenario directory (`YYYY-MM-DD`) with the run filename prefix (`HH-MM-SS`). For
example:

```text
scenario directory: .../2026-04-29/
run filename:       14-24-04_IFE_...pkl
stored attribute:   date = 2026-04-29T14:24:04
```

During campaign ingestion this attribute is indexed in the `.acx`. The
`list_simulation_datetime.sql` query reads it and produces three CSV columns:
`simulation_date`, `simulation_time`, and `simulation_datetime`. Existing BP5
products created with the old hardcoded date must be regenerated and
re-ingested before these columns contain the corrected per-run values.

## Usage

Run from this directory:

```bash
cd RHINO/AI_ready_workflow/3_feature_extraction
```

Build all configured features:

```bash
python build_features.py \
  --acx /path/to/campaign-store/IFE/rhino.acx \
  --campaign-store /path/to/campaign-store/IFE \
  --output outputs/rhino_features.csv
```

`--campaign-store` is the base directory used to resolve archive names stored
inside the index. For example, if the index contains `rhino1.aca`, the command
above reads `/path/to/campaign-store/IFE/rhino1.aca`.

Select only matching archives with a SQL `LIKE` pattern:

```bash
python build_features.py \
  --acx /path/to/rhino.acx \
  --campaign-store /path/to/campaign-store/IFE \
  --archive-name 'rhino1.aca' \
  --output outputs/rhino1_features.csv
```

Use a different specification or query directory:

```bash
python build_features.py \
  --acx /path/to/rhino.acx \
  --campaign-store /path/to/campaign-store/IFE \
  --spec /path/to/feature_spec.json \
  --query-dir /path/to/queries \
  --output outputs/custom_features.csv
```

The defaults for `--spec` and `--query-dir` point to `feature_spec.json` and
`../2_campaign/queries`, respectively. `--archive-name` defaults to `%`, which
selects every archive.

## Inputs and Outputs

Inputs:

- An HPC Campaign index (`.acx`)
- Campaign archives (`.aca`) and their registered data replicas
- A validated feature JSON specification
- SQL files referenced by SQL-backed feature entries

Output:

- One CSV row per selected RHINO campaign run
- Identity columns: `archive`, `datasetid`, and `run_id`
- Metadata, input, and output columns in the same order as their JSON lists

Missing ADIOS variables are reported and represented as missing values in the
CSV. Missing query files, incorrect SQL result columns, duplicate feature keys,
duplicate campaign rows, and out-of-range ADIOS indexes stop the build with a
descriptive error.
