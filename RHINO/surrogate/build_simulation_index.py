#!/usr/bin/env python
# coding: utf-8

# In[1]:


from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import openpmd_api as io


# In[2]:


INPUT_SPECS = [
    {
        "key": "I0_SD",
        "display_name": "Startup inventory [g]",
        "source": "series_input_attribute",
        "attribute_key": "input:I0_SD:Starting inventory"
    },
    {
        "key": "Ndotminus",
        "display_name": "Tritium burning rate [g/d]",
        "source": "series_input_attribute",
        "attribute_key": "input:Ndotminus:Tritium burned per day"
    },
    {
        "key": "beta",
        "display_name": "Burn fraction [-]",
        "source": "series_input_attribute",
        "attribute_key": "input:beta:Burn fraction"
    },
]


OUTPUT_SPECS = [
    {
        "key": "tritium_in_isotope_separation",
        "display_name": "Tritium in isotope separation [g]",
        "source": "particle_record_steady_state",
        "species": "Tritium",
        "subsystem": "Isotope_Seperation",
    },
    {
        "key": "plant_doubling_time_days",
        "display_name": "Plant doubling time [d]",
        "source": "series_output_attribute",
        "attribute_key": "output:plant_doubling_time (days)",
    },
    {
        "key": "minimum_startup_inventory_g",
        "display_name": "Minimum startup inventory [g]",
        "source": "series_output_attribute",
        "attribute_key": "output:I_startup (g)",
    },
]


# In[3]:


def _get_first_snapshot(series):
    snapshots = series.snapshots()
    if not snapshots:
        raise ValueError("The series does not contain any snapshots.")
    return snapshots[0]


def _validate_species(it, species):
    avail_species = [sp for sp in it.particles if sp != "Times"]
    if species == "Times" or species not in avail_species:
        raise ValueError(
            f"Invalid species {species!r}. Available species are: {avail_species}"
        )


def _validate_subsystem(it, species, subsystem):
    avail_subsystems = [sub for sub in it.particles[species]["subsystems"]]
    if subsystem not in avail_subsystems:
        raise ValueError(
            f"Invalid subsystem {subsystem!r}. "
            f"Available subsystems for species {species!r} are: {avail_subsystems}"
        )


def get_series_input_attributes(series):
    values = {}
    for attr in series.attributes:
        if attr.startswith("input:"):
            values[attr] = series.get_attribute(attr)
    return values


def get_series_output_attributes(series):
    values = {}
    for attr in series.attributes:
        if attr.startswith("output:"):
            values[attr] = series.get_attribute(attr)
    return values


def get_steady_state(series, species="Tritium", subsystem=None):
    it = _get_first_snapshot(series)
    _validate_species(it, species)

    mass_steady = it.particles[species]["mass_steady"][io.Record_Component.SCALAR]
    mass_steady_data = mass_steady.load_chunk()
    series.flush()

    if subsystem is None:
        return mass_steady_data

    _validate_subsystem(it, species, subsystem)
    subsystem_id = it.particles[species]["subsystems"][subsystem].get_attribute("id")
    return mass_steady_data[subsystem_id]



# In[4]:


def extract_inputs_from_series(series, input_specs):
    series_inputs = get_series_input_attributes(series)
    values = {}
    for spec in input_specs:
        if spec["source"] != "series_input_attribute":
            raise ValueError(f"Unsupported input source: {spec['source']}")
        values[spec["key"]] = float(series_inputs[spec["attribute_key"]])

    return values


def extract_outputs_from_series(series, output_specs):
    series_outputs = get_series_output_attributes(series)
    values = {}

    for spec in output_specs:
        source = spec["source"]

        if source == "series_output_attribute":
            values[spec["key"]] = float(series_outputs[spec["attribute_key"]])

        elif source == "particle_record_steady_state":
            values[spec["key"]] = float(
                np.asarray(
                    get_steady_state(
                        series,
                        species=spec["species"],
                        subsystem=spec["subsystem"],
                    )
                )
            )

        else:
            raise ValueError(f"Unsupported output source: {source}")

    return values


def build_index_entry(bp5_file: Path, scenario: str) -> dict:
    series = io.Series(
        str(bp5_file),
        io.Access.read_only,
        '{"verify_homogeneous_extents": false}'
    )

    try:
        inputs = extract_inputs_from_series(series, INPUT_SPECS)
        outputs = extract_outputs_from_series(series, OUTPUT_SPECS)
        return {
            "simulation_id": f"{scenario}/{bp5_file.name}",
            "path": str(bp5_file),
            "scenario": scenario,
            "inputs": inputs,
            "outputs": outputs,
        }
    finally:
        series.close()


def build_simulation_index(root_path: Path, scenarios: list[str]) -> list[dict]:
    records = []

    for scenario in scenarios:
        data_path = root_path / scenario
        for bp5_file in sorted(data_path.glob("*.bp5")):
            try:
                record = build_index_entry(bp5_file, scenario)
                records.append(record)
            except Exception as exc:
                print(f"Skipping {bp5_file}: {exc}")

    return records


def main():
    root_path = Path("/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/surrogate_bp_output")
    scenarios = ["2026-04-29", "2026-04-30", "2026-05-01"]

    records = build_simulation_index(root_path, scenarios)

    payload = {
        "input_specs": INPUT_SPECS,
        "output_specs": OUTPUT_SPECS,
        "records": records,
    }

    out_path = Path("simulation_index.json")
    with out_path.open("w") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved {len(records)} records to {out_path}")


if __name__ == "__main__":
    main()


# In[ ]:




