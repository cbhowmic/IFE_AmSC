# RHINO Campaign Management Layer

This layer organizes the openPMD/ADIOS BP5 datasets produced by the shim layer
into HPC Campaign archives and a searchable campaign index. Both Python and
Bash entry points are provided. They perform the same two-stage workflow but
use different configuration files.

## Workflow

```text
RHINO outputs from the shim layer (ADIOS BP5 directories)
    |
    v
Uncompressed TAR files, checksums, and TAR indexes
    |
    v
Campaign archives with live and TAR-backed replicas (.aca)
    |
    v
Campaign index (.acx)
    |
    v
SQL queries and feature extraction
```

The archive stage packages each configured BP5 directory in an uncompressed
TAR, verifies it with SHA-256, assigns each dataset a RHINO run ID, groups the
runs into campaign archives, and registers the TAR copies as archived replicas.
The index stage registers the completed campaign archives in a SQLite-backed
`.acx` file for inspection and querying.

## Requirements

- Python 3.10 or newer for the Python entry points and `hpc-campaign` 0.7
- Bash for the shell entry points
- `tar` and `sha256sum` for the Bash archive entry point
- `hpc_campaign` available on `PATH`
- RHINO BP5 datasets in the configured input directories
- `sqlite3` to run the supplied SQL queries directly

Run the examples below from this directory:

```bash
cd RHINO/AI_ready_workflow/2_campaign
```

## Python Workflow

The Python workflow uses `campaign_spec.json`. Edit that file before running
the scripts for a new campaign or filesystem location.

### `campaign_spec.json`

| Setting | Purpose |
| --- | --- |
| `RHINO_DATA_ROOT` | Root directory containing the RHINO BP5 output directories |
| `CAMPAIGN_STORE` | HPC Campaign storage directory |
| `CAMPAIGN_NAMESPACE` | Logical name retained across the Python and Bash campaign configurations |
| `ARCHIVE_PREFIX` | Prefix for generated archives, such as `rhino` |
| `DATASETS_PER_ARCHIVE` | Maximum number of BP5 datasets placed in each archive |
| `CAMPAIGN_INDEX` | Index path; a relative path is resolved from `CAMPAIGN_STORE` |
| `TAR_OUTPUT_DIR` | TAR storage directory; a relative path is resolved from `RHINO_DATA_ROOT` |
| `TAR_PREFIX` | Prefix for date-based TAR names, such as `rhino` |
| `TAR_STORAGE_SYSTEM` | HPC Campaign archival-storage type, such as `fs` |
| `TAR_STORAGE_HOST` | Short, unique host name recorded for the TAR location |
| `INPUT_DIRS` | Directories to scan under `RHINO_DATA_ROOT` |

Only top-level `*.bp5` entries in each `INPUT_DIRS` directory are discovered;
the search is not recursive.

### `create_archives.py`

`create_archives.py` performs the complete archive stage:

1. Loads and validates the archive settings in the JSON specification.
2. Discovers and sorts BP5 datasets to make archive grouping deterministic.
3. Creates one uncompressed TAR per `INPUT_DIRS` entry, writes its `.sha256`
   checksum, and creates its `.tar.idx` with `hpc_campaign taridx`.
4. Derives the run ID from the timestamp-like portion of each dataset name,
   falling back to the complete filename stem when no timestamp is present.
5. Groups runs according to `DATASETS_PER_ARCHIVE`.
6. Creates `ARCHIVE_PREFIX1.aca`, `ARCHIVE_PREFIX2.aca`, and so on, and adds
   each dataset with its run ID as the campaign dataset name.
7. Registers each relevant TAR and its index with each `.aca`. HPC Campaign
   automatically creates archived replicas for matching BP5 member paths.

Preview dataset discovery, grouping, and generated `hpc_campaign` commands
without creating archives:

```bash
python create_archives.py --dry-run
```

Create the archives:

```bash
python create_archives.py
```

The first dataset in each group is added with `--truncate`, so an archive with
the same generated name is replaced when that group is rebuilt. Existing TARs
are reused only after their SHA-256 checksums pass. To deliberately replace the
TARs, checksums, and TAR indexes, run:

```bash
python create_archives.py --rebuild-tars
```

### `create_index.py`

`create_index.py` performs the index stage:

1. Loads and validates the index settings in the JSON specification.
2. Finds and sorts archives matching `ARCHIVE_PREFIX*.aca` directly inside
   `CAMPAIGN_STORE`.
3. Replaces the configured campaign index.
4. Registers all discovered archives with `hpc_campaign index ... add`.
5. Runs `hpc_campaign index ... ls` to inspect the resulting index.

Preview the index location, discovered archives, and generated commands without
deleting or creating an index:

```bash
python create_index.py --dry-run
```

Build the index after the archives have been created:

```bash
python create_index.py
```

Both Python scripts accept a different JSON specification:

```bash
python create_archives.py --spec /path/to/campaign_spec.json
python create_index.py --spec /path/to/campaign_spec.json
```

The scripts stop immediately if configuration validation or an
`hpc_campaign` command fails.

## Bash Workflow

The original Bash workflow remains available. It reads the equivalent settings
from `config_campaign.sh` rather than from `campaign_spec.json`.

### `config_campaign.sh`

Edit `config_campaign.sh` to set the RHINO data root, campaign store,
namespace, archive prefix, archive size, index path, BP5 input directories, and
TAR storage settings.

### `create_archives.sh`

`create_archives.sh` discovers the configured BP5 datasets, sorts them, derives
their run IDs, creates or verifies the TAR/checksum/index products, divides the
datasets into archive groups, creates each `.aca`, and registers the relevant
TAR-backed replicas.

Run it with:

```bash
bash create_archives.sh
```

The Bash entry point supports the same controls as the Python entry point:

```bash
bash create_archives.sh --dry-run
bash create_archives.sh --rebuild-tars
```

### `create_index.sh`

`create_index.sh` finds archives matching the configured prefix, removes the
existing index, registers the archives in a new index, and lists the index
contents.

Run it after archive creation:

```bash
bash create_index.sh
```

Review `config_campaign.sh` before running the Bash workflow. Use either the
Python workflow or the Bash workflow for a given archive/index build; running
both will recreate the same configured outputs.

## Outputs

With an `ARCHIVE_PREFIX` of `rhino`, archive creation produces files such as:

```text
rhino-2026-04-29.tar
rhino-2026-04-29.tar.sha256
rhino-2026-04-29.tar.idx
rhino1.aca
rhino2.aca
rhino3.aca
```

Index creation produces the `.acx` file configured by `CAMPAIGN_INDEX`. The
index contains queryable RHINO run metadata and attributes and can be inspected
with `hpc_campaign index` or queried directly with `sqlite3`.

## SQL Queries

The `queries/` directory contains reusable queries for campaign inspection,
metadata discovery, scientific analysis, feature extraction, and surrogate
model dataset generation. These queries are used by the feature extraction layer
to compute the features needed for training the surrogate. They include both
reading metadata directly from `.aca` files and reading data from remote
locations.

For the default index located at `<CAMPAIGN_STORE>/rhino.acx`, run a query with:

```bash
sqlite3 /path/to/campaign-store/IFE/rhino.acx < queries/list_archives.sql
```

See `queries/README.md` for the available queries and their expected outputs.
