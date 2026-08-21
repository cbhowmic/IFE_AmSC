"""Single-file Laser HDF5 -> openPMD/ADIOS2 BP5 converter.

This module has one conversion function: laser_to_adios(...). 
Metadata parsing, HDF5 classification, and openPMD helper
writing are delegated to metadata.py and utils.py.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from metadata import load_conversion_metadata, read_laser_h5
from utils import default_bp_filename, get_default_out_dir, remove_path, write_laser_model_to_openpmd


def laser_to_adios(
    h5_path: str | os.PathLike[str],
    out_bp_path: str | os.PathLike[str] | None = None,
    *,
    excel_candidates: list[str | os.PathLike[str]] | None = None,
    documentation_dir: str | os.PathLike[str] | None = None,
    domain_notes_source: str = "",
    overwrite: bool = True,
    strict: bool = False,
    verbose: bool = True,
) -> str:
    """Convert one Laser HDF5 file into one openPMD/ADIOS2 BP5 series.

    Parameters
    ----------
    h5_path:
        Input Laser .h5 file.
    out_bp_path:
        Output .bp5 path. Defaults to bp_output/<input-derived-name>.bp5.
    excel_candidates:
        Optional PV metadata spreadsheet candidates. The first existing file is
        used for Excel enrichment.
    documentation_dir:
        Optional directory containing README/schematic documentation references.
    domain_notes_source:
        Optional path to the laser domain-notes document. The notes are encoded
        in metadata.py as a machine-readable domain manifest; this records the
        source path for provenance.
    overwrite:
        Remove existing output path before writing.
    strict:
        Raise on shape/classification issues instead of skipping problematic
        optional products.
    verbose:
        Print a short conversion summary.

    Returns
    -------
    str
        Path to the written BP5 output.
    """
    h5_path = Path(h5_path)
    if out_bp_path is None:
        out_bp_path = get_default_out_dir() / default_bp_filename(h5_path)
    out_bp_path = Path(out_bp_path)

    if out_bp_path.exists():
        if not overwrite:
            raise FileExistsError(f"Output already exists: {out_bp_path}")
        remove_path(out_bp_path)
    out_bp_path.parent.mkdir(parents=True, exist_ok=True)

    metadata_kwargs: dict[str, Any] = {
        "excel_candidates": excel_candidates,
        "domain_notes_source": domain_notes_source,
        "verbose": verbose,
    }
    if documentation_dir is not None:
        metadata_kwargs["documentation_dir"] = documentation_dir

    excel_meta, documentation_meta = load_conversion_metadata(**metadata_kwargs)
    model = read_laser_h5(h5_path, excel_meta=excel_meta, strict=strict)
    written = write_laser_model_to_openpmd(model, out_bp_path, documentation_meta)

    if verbose:
        print(
            "[done] Wrote Laser BP5: "
            f"{written} "
            f"({model.num_valid_samples}/{model.n_requested} valid samples, "
            f"{len(model.scalars)} scalars, {len(model.traces)} traces, "
            f"{len(model.coordinates)} coordinates, {len(model.images)} images, "
            f"{len(model.timestamps)} timestamps)"
        )
    return written
