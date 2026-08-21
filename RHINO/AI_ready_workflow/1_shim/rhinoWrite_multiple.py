"""Batch conversion utilities for the RHINO openPMD/ADIOS2 shim layer."""

from pathlib import Path

from .rhinoWrite import rhino_to_adios


def convert_scenarios(
    root_path,
    scenarios,
    output_root,
    skip_runs=None,
) -> None:
    """Convert all matching RHINO runs for one or more scenario directories.

    Parameters
    ----------
    root_path
        Directory containing scenario subdirectories.
    scenarios
        Iterable of scenario directory names.
    output_root
        Directory where BP5 output directories/files will be written.
    skip_runs
        Optional set of ``(scenario, prefix)`` pairs to skip.
    """
    root_path = Path(root_path)
    output_root = Path(output_root)
    skip_runs = set(skip_runs or [])

    for scenario in scenarios:
        print(scenario)

        try:
            scenario_path = root_path / scenario
            data_root = scenario_path

            if not data_root.exists():
                print(f"Scenario folder missing or no Data/: {data_root} -- skipping")
                continue

            for tfile in sorted(data_root.glob("*_T_reduced.pkl")):
                try:
                    run_prefix_data = tfile.name.split("_T_reduced.pkl")[0]
                    run_time_prefix = run_prefix_data.split("_")[0]
                    prefix = run_time_prefix
                    infix = "_".join(run_prefix_data.split("_")[1:])

                    if (scenario, prefix) in skip_runs:
                        print(f"Skipping run {prefix} in scenario {scenario}")
                        continue

                    safe_param = scenario.replace(" ", "_").replace("&", "And")
                    output_path = output_root / safe_param / f"{run_time_prefix}.bp5"

                    rhino_to_adios(
                        DATA_PATH=scenario_path,
                        PREFIX=prefix,
                        INFIX=infix,
                        OUTPUT_PATH=output_path,
                    )

                except Exception as exc:
                    print(
                        f"ERROR processing run '{tfile.name}' "
                        f"in scenario '{scenario}': {exc}"
                    )
                    continue

        except Exception as exc:
            print(f"ERROR processing scenario '{scenario}': {exc}")
            continue


def main() -> None:
    """Command-line entry point for batch RHINO shim conversion."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert multiple RHINO runs into openPMD/ADIOS2 BP5 format."
    )
    parser.add_argument(
        "--root-path",
        required=True,
        help="Root directory containing RHINO scenario directories.",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        required=True,
        help="Scenario directory names to process.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        help="Directory where BP5 outputs will be written.",
    )
    parser.add_argument(
        "--skip-run",
        action="append",
        default=[],
        help=(
            "Run to skip, formatted as SCENARIO:PREFIX. "
            "May be supplied multiple times."
        ),
    )

    args = parser.parse_args()

    skip_runs = set()
    for item in args.skip_run:
        try:
            scenario, prefix = item.split(":", 1)
        except ValueError as exc:
            raise ValueError(
                f"Invalid --skip-run value '{item}'. Expected SCENARIO:PREFIX."
            ) from exc

        skip_runs.add((scenario, prefix))

    convert_scenarios(
        root_path=args.root_path,
        scenarios=args.scenarios,
        output_root=args.output_root,
        skip_runs=skip_runs,
    )


if __name__ == "__main__":
    main()
