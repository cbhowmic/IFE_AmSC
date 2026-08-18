"""Command-line entry point for building the RHINO ML feature table."""

from __future__ import annotations

import argparse
from pathlib import Path

from spec_loader import load_feature_spec


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_SPEC = BASE_DIR / "feature_spec.json"
DEFAULT_QUERY_DIR = BASE_DIR.parent / "2_campaign" / "queries"


def parse_args() -> argparse.Namespace:
    """Parse feature-table build options."""
    parser = argparse.ArgumentParser(
        description="Build a RHINO feature table from an HPC Campaign index."
    )
    parser.add_argument(
        "--acx",
        type=Path,
        required=True,
        help="Path to the SQLite-backed campaign index (.acx).",
    )
    parser.add_argument(
        "--campaign-store",
        type=Path,
        required=True,
        help=(
            "Base directory used to resolve archive names stored in the index."
        ),
    )
    parser.add_argument(
        "--spec",
        type=Path,
        default=DEFAULT_SPEC,
        help=f"Feature JSON specification (default: {DEFAULT_SPEC}).",
    )
    parser.add_argument(
        "--query-dir",
        type=Path,
        default=DEFAULT_QUERY_DIR,
        help=f"Directory containing feature SQL queries (default: {DEFAULT_QUERY_DIR}).",
    )
    parser.add_argument(
        "--archive-name",
        default="%",
        help="SQL LIKE pattern used to select archive names (default: %%).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("rhino_features.csv"),
        help="Output CSV file (default: rhino_features.csv).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.acx.is_file():
        raise FileNotFoundError(f"Campaign index does not exist: {args.acx}")
    if not args.campaign_store.is_dir():
        raise FileNotFoundError(
            f"Campaign store does not exist: {args.campaign_store}"
        )
    if not args.query_dir.is_dir():
        raise FileNotFoundError(
            f"Campaign query directory does not exist: {args.query_dir}"
        )

    from campaign_reader import build_feature_table

    specification = load_feature_spec(args.spec)
    frame = build_feature_table(
        acx_path=args.acx,
        query_dir=args.query_dir,
        campaign_store=args.campaign_store,
        archive_name=args.archive_name,
        feature_specs=specification["features"],
    )
    if frame.empty:
        raise ValueError(
            f"No campaign runs matched archive pattern {args.archive_name!r}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)

    print(frame.head())
    print(f"\nRows: {len(frame)}, columns: {len(frame.columns)}")
    missing = frame.isna().sum()
    missing = missing[missing > 0]
    if missing.empty:
        print("Missing values: none")
    else:
        print("Missing values:")
        print(missing.to_string())
    print(f"\nSaved feature table to: {args.output}")


if __name__ == "__main__":
    main()
