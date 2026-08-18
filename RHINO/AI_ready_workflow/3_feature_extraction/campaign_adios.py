"""Read array-backed RHINO features through the ADIOS campaign reader."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable


def build_adios_variable_name(run_id: str, variable: str) -> str:
    """Build the campaign-reader variable name for one represented RHINO run."""
    return f"{run_id}/{variable.lstrip('/')}"


def read_adios_feature(
    reader: Any,
    available_vars: dict[str, Any],
    run_id: str,
    spec: dict[str, Any],
) -> Any:
    """Read and optionally index one ADIOS-backed feature value."""
    variable_name = build_adios_variable_name(run_id, spec["variable"])
    if variable_name not in available_vars:
        print(f"Missing ADIOS variable: {variable_name}")
        return None

    value = reader.read(variable_name)
    if "index" not in spec:
        return value

    try:
        return value[spec["index"]]
    except IndexError as error:
        subsystem = spec.get("subsystem", "requested subsystem")
        raise IndexError(
            f"Index {spec['index']} for {subsystem!r} is outside ADIOS variable "
            f"{variable_name!r}"
        ) from error


def load_adios_feature_table(
    archive_path: str | Path,
    run_ids: Iterable[str],
    spec: dict[str, Any],
) -> list[dict[str, Any]]:
    """Read one ADIOS-backed feature for every requested campaign run."""
    import adios2

    reader = adios2.FileReader(str(archive_path))
    try:
        available_vars = reader.available_variables()
        return [
            {
                "run_id": run_id,
                spec["key"]: read_adios_feature(
                    reader=reader,
                    available_vars=available_vars,
                    run_id=run_id,
                    spec=spec,
                ),
            }
            for run_id in run_ids
        ]
    finally:
        reader.close()
