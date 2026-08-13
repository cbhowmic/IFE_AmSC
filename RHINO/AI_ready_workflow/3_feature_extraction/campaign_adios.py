import adios2
import pandas as pd


def build_adios_variable_name(run_id, variable):
    variable = variable.lstrip("/")
    return f"{run_id}/{variable}"


def read_adios_feature(reader, available_vars, run_id, spec):
    variable_name = build_adios_variable_name(
        run_id=run_id,
        variable=spec["variable"],
    )

    if variable_name not in available_vars:
        print(f"Missing ADIOS variable: {variable_name}")
        return None

    arr = reader.read(variable_name)

    if "index" in spec:
        return arr[spec["index"]]

    return arr


def load_adios_feature_table(archive_path, run_ids, spec):
    rows = []

    reader = adios2.FileReader(archive_path)
    available_vars = reader.available_variables()

    for run_id in run_ids:
        value = read_adios_feature(
            reader=reader,
            available_vars=available_vars,
            run_id=run_id,
            spec=spec,
        )

        rows.append(
            {
                "run_id": run_id,
                spec["key"]: value,
            }
        )

    return rows