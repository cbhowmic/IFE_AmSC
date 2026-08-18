# RHINO Campaign Queries

This directory contains reusable SQL queries for RHINO campaign indexes.

The queries are designed to support:

- campaign inspection and validation
- metadata exploration
- scientific analysis
- feature extraction
- surrogate-model dataset generation

## Query Categories

Current queries include:

- archive and dataset inspection
- attribute discovery
- simulation date and time extraction
- Tritium inventory analysis
- steady-state analysis
- scenario filtering
- run selection

## Usage

Run queries directly against the RHINO campaign index. With the default
configuration, the index is `rhino.acx` inside `CAMPAIGN_STORE`:

```bash
sqlite3 /path/to/campaign-store/IFE/rhino.acx < queries/list_archives.sql
```

Example:

```bash
sqlite3 /path/to/campaign-store/IFE/rhino.acx \
  < queries/steady_state_time_all_runs.sql
```
