"""Batch converter for multiple Laser HDF5 files.

This script discovers input .h5 files and calls laser_to_adios(...) once per file.
All schema and metadata logic stays in laserWrite.py, metadata.py, and utils.py.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from laserWrite import laser_to_adios
from utils import default_bp_filename, get_default_out_dir


def discover_h5_files(
    input_path: str | os.PathLike[str],
    recursive: bool = False,
    exclude_dirs: set[str] | None = None,
) -> list[Path]:
    """Return sorted HDF5 files from a file or directory path.

    When recursive=True, skip output/documentation directories by name so that
    running from the NERSC Laser root does not accidentally walk bp_output or
    documentation.
    """
    path = Path(input_path)
    if path.is_file():
        return [path]

    exclude_dirs = exclude_dirs or {"bp_output", "documentation", "__pycache__"}
    pattern = "**/*.h5" if recursive else "*.h5"
    files = []
    for candidate in path.glob(pattern):
        if any(part in exclude_dirs for part in candidate.parts):
            continue
        files.append(candidate)
    return sorted(files)


def convert_many(
    input_path: str | os.PathLike[str],
    out_dir: str | os.PathLike[str] | None = None,
    *,
    recursive: bool = False,
    excel_candidates: list[str | os.PathLike[str]] | None = None,
    documentation_dir: str | os.PathLike[str] | None = None,
    domain_notes_source: str = "",
    overwrite: bool = True,
    skip_existing: bool = False,
    strict: bool = False,
    continue_on_error: bool = True,
    verbose: bool = True,
) -> list[str]:
    """Convert every discovered .h5 file into an output BP5 series."""
    h5_files = discover_h5_files(input_path, recursive=recursive)
    out_dir = Path(out_dir) if out_dir is not None else get_default_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[str] = []
    failures: list[tuple[Path, Exception]] = []
    skipped: list[Path] = []

    for h5_path in h5_files:
        out_bp_path = out_dir / default_bp_filename(h5_path)
        if skip_existing and out_bp_path.exists():
            skipped.append(out_bp_path)
            if verbose:
                print(f"[skip] Output already exists: {out_bp_path}")
            continue
        try:
            written_path = laser_to_adios(
                h5_path=h5_path,
                out_bp_path=out_bp_path,
                excel_candidates=excel_candidates,
                documentation_dir=documentation_dir,
                domain_notes_source=domain_notes_source,
                overwrite=overwrite,
                strict=strict,
                verbose=verbose,
            )
            written.append(written_path)
        except Exception as exc:
            failures.append((h5_path, exc))
            if not continue_on_error:
                raise
            print(f"[error] Failed to convert {h5_path}: {exc}")

    if verbose:
        print(
            f"[summary] Converted {len(written)} file(s); "
            f"skipped {len(skipped)} existing file(s); "
            f"failed {len(failures)} file(s)."
        )
    return written


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert multiple Laser HDF5 files to openPMD/ADIOS2 BP5.")
    parser.add_argument("input_path", help="Input .h5 file or directory containing .h5 files.")
    parser.add_argument("--out-dir", default=None, help="Output directory for .bp5 files. Defaults to $LASER_BP_OUT_DIR or <LASER_ROOT>/bp_output.")
    parser.add_argument("--recursive", action="store_true", help="Recursively search input directory for .h5 files.")
    parser.add_argument("--excel", action="append", dest="excel_candidates", help="Optional PV metadata Excel path. May be repeated.")
    parser.add_argument("--documentation-dir", help="Optional documentation directory.")
    parser.add_argument("--domain-notes", default="", help="Optional path to laser-pv-details.docx for provenance.")
    parser.add_argument("--no-overwrite", action="store_true", help="Fail if output already exists.")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Leave existing BP5 outputs untouched and continue with remaining inputs.",
    )
    parser.add_argument("--strict", action="store_true", help="Raise on optional shape/classification problems.")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop at the first failed conversion.")
    parser.add_argument("--quiet", action="store_true", help="Suppress routine conversion messages.")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    convert_many(
        input_path=args.input_path,
        out_dir=args.out_dir,
        recursive=args.recursive,
        excel_candidates=args.excel_candidates,
        documentation_dir=args.documentation_dir,
        domain_notes_source=args.domain_notes,
        overwrite=not args.no_overwrite,
        skip_existing=args.skip_existing,
        strict=args.strict,
        continue_on_error=not args.stop_on_error,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
