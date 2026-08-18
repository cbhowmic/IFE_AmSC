"""Load and validate declarative RHINO feature specifications."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


SUPPORTED_SCHEMA_VERSION = 1
SUPPORTED_SOURCES = {"campaign_sql", "campaign_adios"}


def _require(spec: dict[str, Any], fields: set[str], *, context: str) -> None:
    missing = sorted(fields - spec.keys())
    if missing:
        raise ValueError(f"{context} is missing required field(s): {', '.join(missing)}")


def _validate_feature(spec: Any, *, expected_role: str, position: int) -> None:
    context = f"{expected_role} feature #{position}"
    if not isinstance(spec, dict):
        raise TypeError(f"{context} must be a JSON object")

    _require(spec, {"key", "role", "source"}, context=context)
    if not isinstance(spec["key"], str) or not spec["key"]:
        raise TypeError(f"{context} key must be a non-empty string")
    if spec["role"] != expected_role:
        raise ValueError(
            f"{context} role must be {expected_role!r}, not {spec['role']!r}"
        )
    if spec["source"] not in SUPPORTED_SOURCES:
        allowed = ", ".join(sorted(SUPPORTED_SOURCES))
        raise ValueError(f"{context} source must be one of: {allowed}")

    if spec["source"] == "campaign_sql":
        _require(spec, {"query", "column"}, context=context)
        for field in ("query", "column"):
            if not isinstance(spec[field], str) or not spec[field]:
                raise TypeError(f"{context} {field} must be a non-empty string")
    else:
        _require(spec, {"variable"}, context=context)
        if not isinstance(spec["variable"], str) or not spec["variable"]:
            raise TypeError(f"{context} variable must be a non-empty string")
        if "index" in spec and (
            isinstance(spec["index"], bool)
            or not isinstance(spec["index"], int)
            or spec["index"] < 0
        ):
            raise TypeError(f"{context} index must be a non-negative integer")


def load_feature_spec(path: str | Path) -> dict[str, Any]:
    """Load a feature JSON file and return validated metadata/input/output lists."""
    spec_path = Path(path)
    with spec_path.open(encoding="utf-8") as stream:
        document = json.load(stream)

    if not isinstance(document, dict):
        raise TypeError("Feature specification must be a JSON object")
    _require(document, {"schema_version", "inputs", "outputs"}, context="document")
    if document["schema_version"] != SUPPORTED_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported feature specification schema_version "
            f"{document['schema_version']!r}; expected {SUPPORTED_SCHEMA_VERSION}"
        )
    if not isinstance(document["inputs"], list):
        raise TypeError("Feature specification inputs must be a list")
    if not isinstance(document["outputs"], list):
        raise TypeError("Feature specification outputs must be a list")
    metadata = document.get("metadata", [])
    if not isinstance(metadata, list):
        raise TypeError("Feature specification metadata must be a list")
    if not document["inputs"]:
        raise ValueError("Feature specification must define at least one input")
    if not document["outputs"]:
        raise ValueError("Feature specification must define at least one output")

    for position, feature in enumerate(metadata, start=1):
        _validate_feature(feature, expected_role="metadata", position=position)
    for position, feature in enumerate(document["inputs"], start=1):
        _validate_feature(feature, expected_role="input", position=position)
    for position, feature in enumerate(document["outputs"], start=1):
        _validate_feature(feature, expected_role="output", position=position)

    features = [*metadata, *document["inputs"], *document["outputs"]]
    keys = [feature["key"] for feature in features]
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    if duplicates:
        raise ValueError(f"Duplicate feature key(s): {', '.join(duplicates)}")
    if not any(feature["source"] == "campaign_sql" for feature in features):
        raise ValueError(
            "At least one campaign_sql feature is required to establish campaign runs"
        )

    return {
        "schema_version": document["schema_version"],
        "metadata": metadata,
        "inputs": document["inputs"],
        "outputs": document["outputs"],
        "features": features,
    }
