from pathlib import Path
import sqlite3
import pandas as pd

from feature_specs import INPUT_SPECS, OUTPUT_SPECS, FEATURE_SPECS
from campaign_adios import load_adios_feature_table



def read_sql_file(query_dir, query_name):
    query_path = Path(query_dir) / query_name
    with open(query_path, "r") as f:
        return f.read()


def load_sql_feature(acx_path, query_dir, spec, archive_name):
    query = read_sql_file(query_dir, spec["query"])

    with sqlite3.connect(acx_path) as conn:
        df = pd.read_sql_query(
            query,
            conn,
            params={"archive_name": archive_name},
        )

    keep_cols = ["archive", "datasetid", "run_id", spec["column"]]
    df = df[keep_cols].copy()

    df = df.rename(columns={spec["column"]: spec["key"]})

    return df


def build_sql_feature_table(acx_path, query_dir, archive_name):
    specs = INPUT_SPECS + OUTPUT_SPECS

    sql_specs = [
        spec for spec in specs
        if spec["source"] == "campaign_sql"
    ]

    feature_df = None

    for spec in sql_specs:
        df = load_sql_feature(
            acx_path=acx_path,
            query_dir=query_dir,
            spec=spec,
            archive_name=archive_name,
        )

        if feature_df is None:
            feature_df = df
        else:
            feature_df = feature_df.merge(
                df,
                on=["archive", "datasetid", "run_id"],
                how="outer",
            )

    return feature_df


def build_feature_table(acx_path, query_dir, campaign_store, archive_name=None):
    feature_df = build_sql_feature_table(
        acx_path=acx_path,
        query_dir=query_dir,
        archive_name=archive_name,
    )

    adios_specs = [
        spec for spec in FEATURE_SPECS
        if spec["source"] == "campaign_adios"
    ]

    for spec in adios_specs:
        all_rows = []

        for archive, group in feature_df.groupby("archive"):
            archive_path = Path(campaign_store) / archive

            rows = load_adios_feature_table(
                archive_path=str(archive_path),
                run_ids=group["run_id"].tolist(),
                spec=spec,
            )

            for row in rows:
                row["archive"] = archive

            all_rows.extend(rows)

        adios_df = pd.DataFrame(all_rows)

        feature_df = feature_df.merge(
            adios_df,
            on=["archive", "run_id"],
            how="left",
        )

    return feature_df