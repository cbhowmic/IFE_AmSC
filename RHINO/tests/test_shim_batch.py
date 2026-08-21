from rhino.shim import rhinoWrite_multiple


def test_failed_run_does_not_stop_later_runs(tmp_path, monkeypatch, capsys):
    scenario = "2026-04-30"
    scenario_path = tmp_path / "raw" / scenario
    scenario_path.mkdir(parents=True)

    prefixes = ["11-00-38", "11-00-39", "11-00-40"]
    for prefix in prefixes:
        (scenario_path / f"{prefix}_IFE_AmSC_T_reduced.pkl").touch()

    converted = []

    def fake_rhino_to_adios(*, DATA_PATH, PREFIX, INFIX, OUTPUT_PATH):
        converted.append(PREFIX)
        if PREFIX == "11-00-39":
            raise FileNotFoundError("missing steady-state pickle")

    monkeypatch.setattr(
        rhinoWrite_multiple,
        "rhino_to_adios",
        fake_rhino_to_adios,
    )

    rhinoWrite_multiple.convert_scenarios(
        root_path=tmp_path / "raw",
        scenarios=[scenario],
        output_root=tmp_path / "bp5",
    )

    assert converted == prefixes
    output = capsys.readouterr().out
    assert "ERROR processing run '11-00-39_IFE_AmSC_T_reduced.pkl'" in output
    assert f"in scenario '{scenario}': missing steady-state pickle" in output
    assert f"ERROR processing scenario '{scenario}'" not in output
