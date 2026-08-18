"""Build tabular RHINO features from campaign SQL and ADIOS data."""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from campaign_adios import load_adios_feature_table


IDENTITY_COLUMNS = ["archive", "datasetid", "run_id"]


def read_sql_file(query_dir: str | Path, query_name: str) -> str:
    """Read one named SQL query from the configured query directory."""
    query_path = Path(query_dir) / query_name
    if not query_path.is_file():
        raise FileNotFoundError(f"Feature query does not exist: {query_path}")
    return query_path.read_text(encoding="utf-8")


def load_sql_feature(
    acx_path: str | Path,
    query_dir: str | Path,
    spec: dict[str, Any],
    archive_name: str,
) -> pd.DataFrame:
    """Execute one SQL feature specification and normalize its output column."""
    query = read_sql_file(query_dir, spec["query"])
    with sqlite3.connect(acx_path) as connection:
        frame = pd.read_sql_query(
            query,
            connection,
            params={"archive_name": archive_name},
        )

    required_columns = [*IDENTITY_COLUMNS, spec["column"]]
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"Query {spec['query']!r} did not return required column(s): "
            f"{', '.join(missing)}"
        )

    frame = frame[required_columns].copy()
    if frame.duplicated(IDENTITY_COLUMNS).any():
        raise ValueError(
            f"Query {spec['query']!r} returned duplicate campaign run rows"
        )
    return frame.rename(columns={spec["column"]: spec["key"]})


def build_sql_feature_table(
    acx_path: str | Path,
    query_dir: str | Path,
    archive_name: str,
    feature_specs: list[dict[str, Any]],
) -> pd.DataFrame:
    """Outer-join all SQL-backed features by campaign run identity."""
    sql_specs = [
        spec for spec in feature_specs if spec["source"] == "campaign_sql"
    ]
    if not sql_specs:
        raise ValueError("At least one campaign_sql feature is required")

    feature_frame: pd.DataFrame | None = None
    for spec in sql_specs:
        frame = load_sql_feature(
            acx_path=acx_path,
            query_dir=query_dir,
            spec=spec,
            archive_name=archive_name,
        )
        feature_frame = (
            frame
            if feature_frame is None
            else feature_frame.merge(frame, on=IDENTITY_COLUMNS, how="outer")
        )

    if feature_frame is None:
        raise RuntimeError("SQL feature table construction produced no table")
    return feature_frame


def build_feature_table(
    acx_path: str | Path,
    query_dir: str | Path,
    campaign_store: str | Path,
    feature_specs: list[dict[str, Any]],
    archive_name: str = "%",
) -> pd.DataFrame:
    """Build the complete SQL- and ADIOS-backed feature table."""
    feature_frame = build_sql_feature_table(
        acx_path=acx_path,
        query_dir=query_dir,
        archive_name=archive_name,
        feature_specs=feature_specs,
    )

    adios_specs = [
        spec for spec in feature_specs if spec["source"] == "campaign_adios"
    ]
    for spec in adios_specs:
        all_rows: list[dict[str, Any]] = []
        for archive, group in feature_frame.groupby("archive", sort=True):
            archive_path = Path(campaign_store) / archive
            rows = load_adios_feature_table(
                archive_path=archive_path,
                run_ids=group["run_id"].tolist(),
                spec=spec,
            )
            for row in rows:
                row["archive"] = archive
            all_rows.extend(rows)

        adios_frame = pd.DataFrame(
            all_rows,
            columns=["archive", "run_id", spec["key"]],
        )
        feature_frame = feature_frame.merge(
            adios_frame,
            on=["archive", "run_id"],
            how="left",
            validate="one_to_one",
        )

    ordered_columns = [
        *IDENTITY_COLUMNS,
        *(spec["key"] for spec in feature_specs),
    ]
    return feature_frame.reindex(columns=ordered_columns).sort_values(
        IDENTITY_COLUMNS
    ).reset_index(drop=True)
