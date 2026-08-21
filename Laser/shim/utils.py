"""Shared utilities for the Laser openPMD/ADIOS2 shim.

This module intentionally contains low-level helpers only: path handling,
HDF5/string/JSON conversion, valid-sample inference, and small openPMD write
helpers. Laser-specific metadata interpretation lives in metadata.py.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np


SCHEMA_VERSION = "0.0.1"

LASER_FILENAME_RE = re.compile(
    r"^(?P<file_kind>png|alignment)-(?P<run_id>[^-]+)-"
    r"(?P<run_date>\d{4}-\d{2}-\d{2})(?P<run_suffix>.*)$",
    re.IGNORECASE,
)

LASER_FILE_ROLES = {
    "png": "full_system_peening",
    "alignment": "alignment_acquisition",
}

# Default NERSC location for the Laser data campaign.
DEFAULT_NERSC_LASER_ROOT = Path("/global/cfs/cdirs/m3239/2026_FES-AmSC/data/laser")


def get_laser_root() -> Path:
    """Return the Laser data root directory.
    On NERSC this defaults to /global/cfs/cdirs/m3239/2026_FES-AmSC/data/laser.
    Override with LASER_ROOT for local testing or alternate layouts.
    """
    return Path(os.environ.get("LASER_ROOT", str(DEFAULT_NERSC_LASER_ROOT))).expanduser()


def get_laser_data_dir() -> Path:
    """Return the directory containing Laser .h5 files and documentation.
    Older local layouts used LASER_ROOT/Data. The NERSC layout places dated folders and documentation directly under LASER_ROOT. 
    Prefer LASER_ROOT/Data only when it exists; otherwise use LASER_ROOT itself.
    """
    root = get_laser_root()
    legacy_data = root / "Data"
    return legacy_data if legacy_data.exists() else root


def get_default_documentation_dir() -> Path:
    """Return the default documentation directory for PV spreadsheets/notes."""
    return Path(os.environ.get("LASER_DOCUMENTATION_DIR", str(get_laser_data_dir() / "documentation"))).expanduser()


def get_default_out_dir() -> Path:
    """Return the default BP5 output directory."""
    return Path(os.environ.get("LASER_BP_OUT_DIR", str(get_laser_data_dir() / "bp_output"))).expanduser()



# -----------------------------------------------------------------------------
# General path/string helpers
# -----------------------------------------------------------------------------


def find_first_existing(paths: Iterable[str | os.PathLike[str] | None]) -> str | None:
    """Return the first path that exists, or None if no candidate exists."""
    for path in paths:
        if path and Path(path).exists():
            return str(path)
    return None


def sanitize_pv_name(pv_name: str) -> str:
    """Convert an EPICS PV name into a safe openPMD mesh path component."""
    name = (pv_name or "").strip()
    name = name.replace("/", "_").replace(":", "__").replace(" ", "_")
    name = name.replace(".", "_")
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    if not name:
        name = "pv"
    if name[0].isdigit():
        name = f"_{name}"
    return name


def is_dataset(obj: Any) -> bool:
    """Return true when obj is an HDF5 dataset."""
    return isinstance(obj, h5py.Dataset)


def is_group(obj: Any) -> bool:
    """Return true when obj is an HDF5 group."""
    return isinstance(obj, h5py.Group)


def decode_hdf5_value(value: Any) -> str:
    """Decode bytes-like HDF5 values without raising on invalid bytes."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8", errors="replace")
    return str(value)


def deduplicate(items: Iterable[str]) -> list[str]:
    """Deduplicate strings while preserving first-seen order."""
    seen: set[str] = set()
    result: list[str] = []
    for item in items:
        normalized = (item or "").strip()
        if normalized and normalized not in seen:
            seen.add(normalized)
            result.append(normalized)
    return result


def parse_laser_filename(path: str | os.PathLike[str]) -> dict[str, str]:
    """Parse Laser source filenames such as png-<run>-YYYY-MM-DD.h5."""
    stem = Path(path).stem
    info = {
        "source_stem": stem,
        "file_kind": "",
        "file_role": "unknown",
        "run_id": "",
        "run_date": "",
        "run_suffix": "",
    }

    match = LASER_FILENAME_RE.match(stem)
    if match:
        file_kind = match.group("file_kind").lower()
        info.update(
            {
                "file_kind": file_kind,
                "file_role": LASER_FILE_ROLES.get(file_kind, "unknown"),
                "run_id": match.group("run_id"),
                "run_date": match.group("run_date"),
                "run_suffix": match.group("run_suffix") or "",
            }
        )
        return info

    prefix = stem.split("-", 1)[0].lower()
    if prefix in LASER_FILE_ROLES:
        info["file_kind"] = prefix
        info["file_role"] = LASER_FILE_ROLES[prefix]

    date_match = re.search(r"\d{4}-\d{2}-\d{2}", stem)
    if date_match:
        info["run_date"] = date_match.group(0)

    return info


def default_bp_filename(h5_path: str | os.PathLike[str]) -> str:
    """Build the default BP5 filename for one Laser HDF5 source file."""
    info = parse_laser_filename(h5_path)
    if info["run_id"] and info["run_date"]:
        return f"{info['file_kind']}-{info['run_id']}-{info['run_date']}{info['run_suffix']}.bp5"
    return f"{Path(h5_path).stem}.bp5"


def remove_path(path: str | os.PathLike[str]) -> None:
    """Remove a file or directory tree if it exists."""
    path_obj = Path(path)
    if path_obj.is_dir():
        shutil.rmtree(path_obj)
    elif path_obj.exists():
        path_obj.unlink()


def current_timestamp_strings() -> tuple[str, str]:
    """Return local and UTC timestamp strings for provenance attributes."""
    now_local = datetime.now().astimezone()
    now_utc = datetime.now(timezone.utc)
    return (
        now_local.strftime("%Y-%m-%d %H:%M:%S %Z"),
        now_utc.strftime("%Y-%m-%d %H:%M:%S UTC"),
    )


# -----------------------------------------------------------------------------
# JSON/HDF5 conversion helpers
# -----------------------------------------------------------------------------


def to_jsonable(value: Any) -> Any:
    """Convert NumPy/HDF5-ish values to JSON-serializable Python values."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        as_float = float(value)
        if np.isnan(as_float):
            return None
        return as_float
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (bytes, np.bytes_)):
        return decode_hdf5_value(value)
    if isinstance(value, np.ndarray):
        return [to_jsonable(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): to_jsonable(val) for key, val in value.items()}
    if value is None:
        return None
    return value


def encode_text_blob(text: str) -> np.ndarray:
    """Encode text as a uint8 array suitable for an openPMD record component."""
    return np.frombuffer((text or "").encode("utf-8"), dtype=np.uint8).copy()


def encode_json_blob(obj: Any) -> np.ndarray:
    """Encode JSON-serializable data as a uint8 array."""
    text = json.dumps(to_jsonable(obj), ensure_ascii=False, indent=2)
    return encode_text_blob(text)


def parse_json_list(value: Any) -> list[str]:
    """Parse a JSON-encoded HDF5 attribute list, returning [] on failure."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, np.ndarray)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = decode_hdf5_value(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return [text]
    if isinstance(parsed, list):
        return [str(item).strip() for item in parsed if str(item).strip()]
    return [str(parsed).strip()] if parsed else []


def h5_attrs_to_dict(attrs: h5py.AttributeManager) -> dict[str, Any]:
    """Return JSON-safe HDF5 attributes."""
    return {str(key): to_jsonable(value) for key, value in attrs.items()}


def get_root_datasets(root: h5py.File) -> dict[str, h5py.Dataset]:
    """Return datasets stored directly at the HDF5 root."""
    return {name: value for name, value in root.items() if is_dataset(value)}


def get_image_groups(root: h5py.File) -> dict[str, h5py.Group]:
    """Return groups containing datasets named image 0, image 1, ..."""
    image_groups: dict[str, h5py.Group] = {}
    for name, value in root.items():
        if not is_group(value):
            continue
        for dataset_name, dataset in value.items():
            if is_dataset(dataset) and re.match(r"^image\s+\d+$", dataset_name):
                image_groups[name] = value
                break
    return image_groups


def infer_requested_sample_count(
    root_attrs: Mapping[str, Any],
    root_datasets: Mapping[str, h5py.Dataset],
) -> int:
    """Infer the requested/source-buffer sample count."""
    if root_attrs.get("N_requested") not in (None, ""):
        try:
            return int(root_attrs["N_requested"])
        except (TypeError, ValueError):
            pass

    first_dims: list[int] = []
    for dataset in root_datasets.values():
        if dataset.ndim >= 1:
            first_dims.append(int(dataset.shape[0]))
    if not first_dims:
        return 0
    return max(set(first_dims), key=first_dims.count)


def infer_valid_sample_mask(
    root_datasets: Mapping[str, h5py.Dataset],
    image_groups: Mapping[str, h5py.Group],
    n_requested: int,
) -> np.ndarray:
    """Infer which rows in the HDF5 sample buffer contain acquired data.

    Preference order:
    1. phoeniX:epoch, because it is the acquisition/sample clock.
    2. PNG:ShotNumber and PNG:RunNumber, because they are expected to be zero
       in padded rows in the current files.
    3. image frame count, when images are present.
    4. all rows valid.
    """
    if n_requested <= 0:
        return np.zeros(0, dtype=bool)

    candidate_names = ["phoeniX:epoch", "epoch", "PNG:ShotNumber", "PNG:RunNumber"]
    for name in candidate_names:
        dataset = root_datasets.get(name)
        if dataset is None or dataset.ndim != 1 or dataset.shape[0] != n_requested:
            continue
        data = np.asarray(dataset[()])
        if np.issubdtype(data.dtype, np.number):
            mask = np.isfinite(data) & (data != 0)
            if 0 < int(mask.sum()) <= n_requested:
                return mask.astype(bool)

    frame_counts = [count_image_frames(group) for group in image_groups.values()]
    frame_counts = [count for count in frame_counts if count > 0]
    if frame_counts:
        n_valid = min(max(frame_counts), n_requested)
        mask = np.zeros(n_requested, dtype=bool)
        mask[:n_valid] = True
        return mask

    return np.ones(n_requested, dtype=bool)


def count_image_frames(group: h5py.Group) -> int:
    """Return the number of image <index> datasets in an image group."""
    return len(image_dataset_indices(group))


def image_dataset_indices(group: h5py.Group) -> list[int]:
    """Return sorted integer frame indices for datasets named image <index>."""
    indices: list[int] = []
    for name, dataset in group.items():
        if not is_dataset(dataset):
            continue
        match = re.match(r"^image\s+(\d+)$", name)
        if match:
            indices.append(int(match.group(1)))
    return sorted(indices)


def trim_to_valid_samples(data: np.ndarray, valid_mask: np.ndarray) -> np.ndarray:
    """Trim an array's first axis to valid sample rows when possible."""
    if data.ndim >= 1 and data.shape[0] == valid_mask.shape[0]:
        return np.asarray(data[valid_mask])
    return np.asarray(data)


# -----------------------------------------------------------------------------
# openPMD helpers
# -----------------------------------------------------------------------------


def _openpmd_api():
    """Import openpmd_api lazily so metadata-only tests do not need it."""
    try:
        import openpmd_api as io  # type: ignore
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "openpmd_api is required to write BP5 output. Install openPMD-api "
            "in the conversion environment before running laser_to_adios."
        ) from exc
    return io


def make_series(out_bp_path: str | os.PathLike[str]):
    """Create an openPMD Series using ADIOS2 BP5."""
    io = _openpmd_api()
    adios2_config = r'''
    {
      "iteration_encoding": "variable_based",
      "adios2": {
        "engine": {
          "type": "bp5",
          "parameters": {
            "StatsLevel": "1"
          }
        }
      }
    }
    '''
    return io.Series(str(out_bp_path), io.Access_Type.create_linear, adios2_config)


def setup_mesh(mesh: Any, axis_labels: Sequence[str]) -> None:
    """Set common mesh geometry metadata."""
    io = _openpmd_api()
    mesh.geometry = io.Geometry.cartesian
    mesh.axis_labels = list(axis_labels)
    mesh.grid_spacing = [1.0] * len(axis_labels)
    mesh.grid_global_offset = [0.0] * len(axis_labels)


def write_record(mesh: Any, record_name: str, data: np.ndarray) -> None:
    """Write one NumPy array into one openPMD record component."""
    io = _openpmd_api()
    array = np.ascontiguousarray(data)
    record = mesh[record_name]
    record.reset_dataset(io.Dataset(array.dtype, array.shape))
    record.store_chunk(array)


def set_attr_safe(target: Any, key: str, value: Any) -> None:
    """Set an openPMD attribute after converting to supported scalar/string values."""
    if value is None:
        return
    value = to_jsonable(value)
    if isinstance(value, float) and np.isnan(value):
        return
    if isinstance(value, (dict, list, tuple)):
        target.set_attribute(key, json.dumps(value, ensure_ascii=False))
    elif isinstance(value, (str, int, float, bool)):
        target.set_attribute(key, value)
    else:
        target.set_attribute(key, str(value))


def set_attrs(target: Any, attrs: Mapping[str, Any]) -> None:
    """Set many openPMD attributes safely."""
    for key, value in attrs.items():
        set_attr_safe(target, str(key), value)


def write_array_mesh(
    iteration: Any,
    mesh_path: str,
    record_name: str,
    data: np.ndarray,
    axis_labels: Sequence[str],
    attrs: Mapping[str, Any] | None = None,
) -> Any:
    """Create a mesh, write one record component, and attach attributes."""
    mesh = iteration.meshes[mesh_path]
    setup_mesh(mesh, axis_labels)
    write_record(mesh, record_name, data)
    if attrs:
        set_attrs(mesh, attrs)
    return mesh


def write_json_mesh(
    iteration: Any,
    mesh_path: str,
    obj: Any,
    label: str,
    attrs: Mapping[str, Any] | None = None,
) -> Any:
    """Write a JSON object as a UTF-8 uint8 mesh."""
    merged_attrs = {
        "label": label,
        "encoding": "utf-8",
        "contentType": "application/json",
    }
    if attrs:
        merged_attrs.update(dict(attrs))
    return write_array_mesh(
        iteration=iteration,
        mesh_path=mesh_path,
        record_name="json",
        data=encode_json_blob(obj),
        axis_labels=["byte"],
        attrs=merged_attrs,
    )

# -----------------------------------------------------------------------------
# Laser model -> openPMD writer helpers
# -----------------------------------------------------------------------------


EPICS_ATTR_MAP = {
    "pvname": "pv:name",
    "description": "hdf5:description",
    "units": "unitLabel",
    "precision": "epics:precision",
    "lower_ctrl_limit": "epics:lower_ctrl_limit",
    "upper_ctrl_limit": "epics:upper_ctrl_limit",
    "lower_alarm_limit": "epics:lower_alarm_limit",
    "upper_alarm_limit": "epics:upper_alarm_limit",
    "lower_warning_limit": "epics:lower_warning_limit",
    "upper_warning_limit": "epics:upper_warning_limit",
    "lower_disp_limit": "epics:lower_disp_limit",
    "upper_disp_limit": "epics:upper_disp_limit",
    "enum_strs": "epics:enum_strs",
    "type": "epics:type",
    "count": "epics:count",
    "host": "epics:host",
    "access": "epics:access",
}


def pv_mesh_attrs(pv: Any, data_role: str, axis_name: str = "sample") -> dict[str, Any]:
    """Build common openPMD mesh attributes for a PVData-like object."""
    h5_attrs = getattr(pv, "attrs", {}) or {}
    excel = getattr(pv, "excel", {}) or {}
    description = h5_attrs.get("description") or excel.get("excel:description") or f"Laser {data_role} PV"
    attrs: dict[str, Any] = {
        "description": description,
        "pv:name": getattr(pv, "name", ""),
        "dataRole": data_role,
        "sampleAxis": 0 if axis_name == "sample" else None,
        "metadata:sources": [source for source, present in [("hdf5", bool(h5_attrs)), ("excel", bool(excel))] if present],
        "semantic:roles": getattr(pv, "semantic_roles", []),
    }
    for h5_key, openpmd_key in EPICS_ATTR_MAP.items():
        if h5_key in h5_attrs:
            attrs[openpmd_key] = h5_attrs[h5_key]
    for excel_key, value in excel.items():
        attrs[excel_key] = value
    if "unitLabel" not in attrs or attrs.get("unitLabel") in (None, ""):
        attrs["unitLabel"] = "unknown"
        attrs["unitStatus"] = "unmapped"
    else:
        attrs["unitStatus"] = "from_hdf5"
    return attrs


def set_series_attributes(series: Any, model: Any, documentation_meta: Mapping[str, Any]) -> None:
    """Set series-level openPMD/input/output/provenance attributes."""
    creation_date, creation_time_utc = current_timestamp_strings()
    h5_path = Path(model.h5_path)
    file_info = dict(model.file_info)
    root_attrs = dict(model.root_attrs)

    set_attrs(
        series,
        {
            "software": "Laser",
            "softwareVersion": "unknown",
            "softwareDescription": "Laser experimental acquisition conversion to openPMD/ADIOS2 BP5",
            "author": "IFE_AmSC Laser data pipeline",
            "authorAffiliation": "LANL",
            "date": creation_date.split()[0],
            "facility": "LANL",
            "instrument": "IFE_AmSC Laser",
            "experiment": "Laser Peening",
            "comment": f"Converted from Laser HDF5 source {h5_path.name}",
            "laser:schema_version": SCHEMA_VERSION,
            "input:source_h5": h5_path.name,
            "input:source_path": str(h5_path.resolve()),
            "input:file_kind": file_info.get("file_kind", ""),
            "input:file_role": file_info.get("file_role", "unknown"),
            "input:run_id": file_info.get("run_id", ""),
            "input:run_date": file_info.get("run_date", ""),
            "input:run_suffix": file_info.get("run_suffix", ""),
            "input:created_unix_time": root_attrs.get("created_unix_time"),
            "input:created_iso": root_attrs.get("created_iso"),
            "input:source_filename_recorded": root_attrs.get("filename"),
            "input:trigger_pv": root_attrs.get("trigger"),
            "input:hostname": root_attrs.get("hostname"),
            "input:user": root_attrs.get("user"),
            "input:script_name": root_attrs.get("script_name"),
            "input:script_path": root_attrs.get("script_path"),
            "input:num_requested_samples": int(model.n_requested),
            "input:num_valid_samples": int(model.num_valid_samples),
            "input:num_declared_scalar_pvs": len(model.scalar_pvs_declared),
            "input:num_declared_array_pvs": len(model.array_pvs_declared),
            "input:num_declared_image_pvs": len(model.image_pvs_declared),
            "input:pv_metadata_source": documentation_meta.get("pv_metadata_source", ""),
            "output:numScalarPVs": len(model.scalars),
            "output:numTracePVs": len(model.traces),
            "output:numCoordinatePVs": len(model.coordinates),
            "output:numImagePVsPresent": len(model.images),
            "output:numImagePVsDeclared": len(model.image_pvs_declared),
            "output:numTimestampPVs": len(model.timestamps),
            "output:dataProducts": [
                "samples",
                "shots",
                "scalars",
                "traces",
                "coordinates",
                "images",
                "timestamps",
                "metadata",
            ],
            "output:schemaVersion": SCHEMA_VERSION,
            "provenance:converter": "laserWrite.py",
            "provenance:converterSchemaVersion": SCHEMA_VERSION,
            "provenance:creationDate": creation_date,
            "provenance:creationTimeUTC": creation_time_utc,
            "provenance:originalDataDirectory": str(h5_path.parent.resolve()),
            "provenance:originalDataFiles": [h5_path.name],
        },
    )


def write_metadata_meshes(iteration: Any, model: Any, documentation_meta: Mapping[str, Any]) -> None:
    """Write JSON metadata manifests."""
    write_json_mesh(
        iteration,
        "metadata/source/root_attrs",
        model.root_attrs,
        "Original HDF5 root attributes",
    )
    write_json_mesh(
        iteration,
        "metadata/source/file_manifest",
        {
            "source_h5": str(model.h5_path),
            "file_info": model.file_info,
            "documentation": dict(documentation_meta),
        },
        "Source file and documentation manifest",
    )
    write_json_mesh(iteration, "metadata/pv_manifest", model.pv_manifest, "PV classification manifest")
    write_json_mesh(iteration, "metadata/image_manifest", model.image_manifest, "Image/camera manifest")
    write_json_mesh(
        iteration,
        "metadata/control_response_manifest",
        model.conversion_metadata.get("domain_manifest", {}).get("control_response", []),
        "Laser control-response relationship manifest",
    )
    write_json_mesh(
        iteration,
        "metadata/full_system_criteria",
        model.conversion_metadata.get("domain_manifest", {}).get("full_system_criteria", {}),
        "Full-system shot criteria",
    )
    write_json_mesh(
        iteration,
        "metadata/domain_notes",
        model.conversion_metadata.get("domain_manifest", {}),
        "Domain metadata derived from laser notes",
    )
    write_json_mesh(iteration, "metadata/conversion", model.conversion_metadata, "Conversion metadata")


def write_sample_meshes(iteration: Any, model: Any) -> None:
    """Write universal acquisition/sample coordinate meshes."""
    write_array_mesh(
        iteration,
        "samples/index",
        "value",
        np.arange(model.num_valid_samples, dtype=np.int64),
        ["sample"],
        {"description": "Valid acquisition sample index", "unitLabel": "index"},
    )
    write_array_mesh(
        iteration,
        "samples/source_index",
        "value",
        np.asarray(model.source_indices, dtype=np.int64),
        ["sample"],
        {
            "description": "Original source row index in the HDF5 acquisition buffer",
            "unitLabel": "index",
        },
    )
    write_array_mesh(
        iteration,
        "samples/valid",
        "value",
        np.asarray(model.valid_sample_mask, dtype=np.uint8),
        ["source_sample"],
        {
            "description": "Boolean validity mask over the original HDF5 requested sample buffer",
            "unitLabel": "bool",
            "trueValue": 1,
            "falseValue": 0,
        },
    )
    epoch = model.scalars.get("phoeniX:epoch") or model.scalars.get("epoch")
    if epoch is not None:
        write_array_mesh(
            iteration,
            "samples/epoch",
            "value",
            np.asarray(epoch.data),
            ["sample"],
            {
                "description": "Acquisition/sample Unix epoch time",
                "unitLabel": "s",
                "timeReference": "Unix epoch",
                "pv:name": epoch.name,
            },
        )


def write_shot_meshes(iteration: Any, model: Any) -> None:
    """Write optional shot-level convenience meshes when shot PVs exist."""
    if "PNG:RunNumber" in model.scalars:
        write_array_mesh(
            iteration,
            "shots/run_number",
            "value",
            np.asarray(model.scalars["PNG:RunNumber"].data),
            ["sample"],
            {"description": "Run number per valid sample", "pv:name": "PNG:RunNumber"},
        )
    if "PNG:ShotNumber" in model.scalars:
        write_array_mesh(
            iteration,
            "shots/shot_number",
            "value",
            np.asarray(model.scalars["PNG:ShotNumber"].data),
            ["sample"],
            {"description": "Shot number per valid sample", "pv:name": "PNG:ShotNumber"},
        )
    if model.full_system_mask is not None:
        write_array_mesh(
            iteration,
            "shots/full_system_mask",
            "value",
            np.asarray(model.full_system_mask, dtype=np.uint8),
            ["sample"],
            {
                "description": "Mask of samples satisfying main-amplifier/full-system criteria",
                "unitLabel": "bool",
                "trueValue": 1,
                "falseValue": 0,
            },
        )


def write_scalar_meshes(iteration: Any, model: Any) -> None:
    """Write scalar PV arrays."""
    for pv_name, pv in model.scalars.items():
        mesh_path = f"scalars/{sanitize_pv_name(pv_name)}"
        write_array_mesh(
            iteration,
            mesh_path,
            "value",
            np.asarray(pv.data),
            ["sample"],
            pv_mesh_attrs(pv, "scalar"),
        )


def write_coordinate_meshes(iteration: Any, model: Any) -> None:
    """Write coordinate arrays independently so they are not lost."""
    for pv_name, pv in model.coordinates.items():
        role = getattr(pv, "data_role", "coordinate") or "coordinate"
        axis_labels = ["sample", role.split(":", 1)[-1]] if np.asarray(pv.data).ndim == 2 else ["sample"]
        write_array_mesh(
            iteration,
            f"coordinates/{sanitize_pv_name(pv_name)}",
            "value",
            np.asarray(pv.data),
            axis_labels,
            pv_mesh_attrs(pv, role),
        )


def write_trace_meshes(iteration: Any, model: Any) -> None:
    """Write waveform/trace arrays and attach matched coordinates when available."""
    for pv_name, pv in model.traces.items():
        mesh_path = f"traces/{sanitize_pv_name(pv_name)}"
        mesh = iteration.meshes[mesh_path]
        setup_mesh(mesh, ["sample", "point"])
        write_record(mesh, "signal", np.asarray(pv.data))
        attrs = pv_mesh_attrs(pv, "trace")
        attrs.update({"timeAxis": 1, "pointAxis": 1})
        if pv.coordinate_data is not None and pv.coordinate_name is not None:
            record_name = pv.coordinate_role or "coordinate"
            if record_name == "histogram_bin":
                record_name = "bin"
            write_record(mesh, record_name, np.asarray(pv.coordinate_data))
            attrs["coordinate:name"] = pv.coordinate_name
            attrs["coordinate:role"] = pv.coordinate_role or "coordinate"
        set_attrs(mesh, attrs)


def write_image_meshes(iteration: Any, model: Any) -> None:
    """Write image stacks and frame-level metadata."""
    for image_name, image in model.images.items():
        safe = sanitize_pv_name(image_name)
        base_attrs = {
            "description": "Image stack for camera/image PV group",
            "pv:name": image.name,
            "dataRole": "image",
            "frameAxis": 0,
            "imageYAxis": 1,
            "imageXAxis": 2,
            "unitLabel": "counts",
            "numFrames": int(image.frame_index.size),
            "declaredInRootImagePVs": bool(image.declared),
            "semantic:roles": image.semantic_roles,
            "timestamp:meaning": "script_read_attempt_time_not_camera_exposure_time",
        }
        write_array_mesh(iteration, f"images/{safe}/image", "value", image.data, ["frame", "y", "x"], base_attrs)
        write_array_mesh(
            iteration,
            f"images/{safe}/frame_index",
            "value",
            image.frame_index,
            ["frame"],
            {"description": "Image frame index from HDF5 dataset name image <index>", "pv:name": image.name},
        )
        write_array_mesh(
            iteration,
            f"images/{safe}/timestamp",
            "value",
            image.timestamp,
            ["frame"],
            {
                "description": "Timestamp stored on each image dataset; for CCDs this is script read time",
                "unitLabel": "s",
                "timeReference": "Unix epoch",
                "pv:name": image.name,
            },
        )
        if image.array_counter is not None and image.array_counter_pv:
            write_array_mesh(
                iteration,
                f"images/{safe}/array_counter",
                "value",
                image.array_counter,
                ["frame"],
                {
                    "description": "Camera array counter used to detect new camera acquisitions",
                    "pv:name": image.array_counter_pv,
                    "validityRule": "new frame if array counter increments",
                },
            )
        if image.new_frame_mask is not None:
            write_array_mesh(
                iteration,
                f"images/{safe}/new_frame_mask",
                "value",
                np.asarray(image.new_frame_mask, dtype=np.uint8),
                ["frame"],
                {
                    "description": "Boolean frame-validity mask derived from array counter increments",
                    "unitLabel": "bool",
                    "trueValue": 1,
                    "falseValue": 0,
                },
            )


def write_timestamp_meshes(iteration: Any, model: Any) -> None:
    """Write PV timestamp arrays from the /timestamps group."""
    for source_name, pv in model.timestamps.items():
        mesh_path = f"timestamps/{sanitize_pv_name(source_name)}"
        attrs = pv_mesh_attrs(pv, "timestamp")
        attrs.update(
            {
                "description": "PV last-update timestamp per valid sample",
                "meaning": "last_update_time_of_pv",
                "sourceGroup": "timestamps",
                "sourceDataset": source_name,
                "unitLabel": "s",
                "timeReference": "Unix epoch",
            }
        )
        write_array_mesh(iteration, mesh_path, "timestamp", np.asarray(pv.data), ["sample"], attrs)


def write_laser_model_to_openpmd(
    model: Any,
    out_bp_path: str | os.PathLike[str],
    documentation_meta: Mapping[str, Any],
) -> str:
    """Serialize a LaserH5Model into one openPMD/ADIOS2 BP5 series."""
    series = make_series(out_bp_path)
    set_series_attributes(series, model, documentation_meta)

    iteration = series.iterations[0]
    iteration.time = 0.0
    iteration.dt = 1.0
    iteration.time_unit_SI = 1.0
    set_attrs(
        iteration,
        {
            "timeUnitLabel": "sample",
            "timeUnitSI": 1.0,
            "description": "One converted Laser HDF5 acquisition/run stored as iteration 0",
        },
    )

    write_metadata_meshes(iteration, model, documentation_meta)
    write_sample_meshes(iteration, model)
    write_shot_meshes(iteration, model)
    write_scalar_meshes(iteration, model)
    write_coordinate_meshes(iteration, model)
    write_trace_meshes(iteration, model)
    write_image_meshes(iteration, model)
    write_timestamp_meshes(iteration, model)

    iteration.close()
    series.close()
    return str(out_bp_path)
