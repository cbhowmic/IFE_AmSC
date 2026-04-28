#!/usr/bin/env python3

import argparse
import re
import sqlite3
from pathlib import Path

import pandas as pd


INPUT_SPECS = [
    {
        "key": "Ndotminus",
        "display_name": "Tritium burning rate [g/d]",
        "source": "series_input_attribute",
        "attribute_key": "/input:Ndotminus:Tritium burned per day",
    },
    {
        "key": "beta",
        "display_name": "Burn fraction [-]",
        "source": "series_input_attribute",
        "attribute_key": "/input:beta:Burn fraction",
    },
]


OUTPUT_SPECS = [
    {
        "key": "tritium_in_isotope_separation_g",
        "display_name": "Tritium in isotope separation [g]",
        "source": "particle_record_steady_state",
        "species": "Tritium",
        "subsystem": "Isotope_Seperation",
        "attribute_key": "/data/inventory/Tritium/subsystems/Isotope_Seperation",
    },
    {
        "key": "plant_doubling_time_days",
        "display_name": "Plant doubling time [d]",
        "source": "series_output_attribute",
        "attribute_key": "/output:plant_doubling_time (days)",
    },
    {
        "key": "minimum_startup_inventory_g",
        "display_name": "Minimum startup inventory [g]",
        "source": "series_output_attribute",
        "attribute_key": "/output:I_startup (g)",
    },
]


def clean_numeric(value):
    if value is None:
        return None

    s = str(value).strip()

    # Handles values like:
    # [26.557]
    # {272.634}
    # 7.8488
    # "7.8488"
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not match:
        return None

    return float(match.group(0))


def build_training_dataset(db_path, archive_like):
    specs = INPUT_SPECS + OUTPUT_SPECS

    attr_to_key = {
        spec["attribute_key"]: spec["key"]
        for spec in specs
        if "attribute_key" in spec
    }

    attr_names = list(attr_to_key.keys())
    placeholders = ",".join(["?"] * len(attr_names))

    query = f"""
    SELECT
      d.datasetid,
      d.dsid,
      d.name AS dataset,
      a.name AS archive,
      at.name AS attribute_name,
      at.value AS raw_value
    FROM datasets d
    JOIN archives a
      ON d.archiveid = a.archiveid
    JOIN attributes at
      ON at.archiveid = d.archiveid
     AND at.datasetid = d.datasetid
    WHERE a.name LIKE ?
      AND at.name IN ({placeholders})
    ORDER BY a.name, d.datasetid, at.name;
    """

    params = [archive_like] + attr_names

    with sqlite3.connect(db_path) as conn:
        long_df = pd.read_sql_query(query, conn, params=params)

    long_df["variable"] = long_df["attribute_name"].map(attr_to_key)
    long_df["value"] = long_df["raw_value"].apply(clean_numeric)

    wide_df = (
        long_df.pivot_table(
            index=["archive", "datasetid", "dsid", "dataset"],
            columns="variable",
            values="value",
            aggfunc="first",
        )
        .reset_index()
    )

    wide_df.columns.name = None

    ordered_columns = ["archive", "datasetid", "dsid", "dataset"]
    ordered_columns += [spec["key"] for spec in INPUT_SPECS]
    ordered_columns += [spec["key"] for spec in OUTPUT_SPECS]

    ordered_columns = [c for c in ordered_columns if c in wide_df.columns]
    wide_df = wide_df[ordered_columns]

    return wide_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db",
        default="IFE/rhino.acx",
        help="Path to RHINO SQLite campaign index",
    )
    parser.add_argument(
        "--archive-like",
        default="IFE/rhino%.aca",
        help="Archive filter",
    )
    parser.add_argument(
        "--output",
        default="rhino_surrogate_dataset.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--drop-missing",
        action="store_true",
        help="Drop rows with missing input/output values",
    )

    args = parser.parse_args()

    df = build_training_dataset(args.db, args.archive_like)

    if args.drop_missing:
        model_columns = [spec["key"] for spec in INPUT_SPECS + OUTPUT_SPECS]
        model_columns = [c for c in model_columns if c in df.columns]
        df = df.dropna(subset=model_columns)

    output_path = Path(args.output)
    df.to_csv(output_path, index=False)

    print(f"Wrote {len(df)} rows to {output_path}")
    print(df.head())


if __name__ == "__main__":
    main()