"""Count LASER BP5 acquisitions and valid shots recorded by the shim."""

from __future__ import annotations

import argparse
from pathlib import Path

import openpmd_api as io


def count_outputs(output_dir: Path) -> tuple[int, int]:
    datasets = sorted(output_dir.glob("*.bp5"))
    total_shots = 0

    for path in datasets:
        series = io.Series(str(path), io.Access_Type.read_only)
        try:
            shots = int(series.get_attribute("input:num_valid_samples"))
        finally:
            series.close()
        total_shots += shots
        print(f"{path.name}: {shots} valid shot(s)")

    print(f"\nBP5 acquisitions: {len(datasets)}")
    print(f"Total valid shots: {total_shots}")
    return len(datasets), total_shots


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count LASER BP5 acquisitions and their valid shots."
    )
    parser.add_argument("output_dir", type=Path, help="Directory containing *.bp5 outputs")
    args = parser.parse_args()

    if not args.output_dir.is_dir():
        parser.error(f"output directory does not exist: {args.output_dir}")
    count_outputs(args.output_dir)


if __name__ == "__main__":
    main()
