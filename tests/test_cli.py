"""The CLI surface, including the prerequisites that used to raise tracebacks."""

from __future__ import annotations

import pytest

from chembench.cli import build_parser, main


def test_report_only_does_not_rerun_the_grid(tmp_path):
    """The flag was declared and ignored, so asking for a report re-ran a three-hour job."""
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "evaluate", "--report-only"])
    assert "no findings" in str(excinfo.value)


def test_evaluate_without_curated_data_names_the_fix(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--data-dir", str(tmp_path), "--results-dir", str(tmp_path), "evaluate"])
    message = str(excinfo.value)
    assert "no curated targets" in message
    assert "chembench curate" in message


def test_both_subcommands_parse():
    parser = build_parser()
    assert parser.parse_args(["curate"]).command == "curate"
    assert parser.parse_args(["evaluate", "--report-only"]).report_only is True


def test_curate_activity_type_is_settable_and_defaults_to_ki():
    """--targets used to hardcode Ki, which silently fetched the wrong assay type."""
    parser = build_parser()
    assert parser.parse_args(["curate", "--targets", "CHEMBL228"]).activity_type == "Ki"
    args = parser.parse_args(["curate", "--targets", "CHEMBL228", "--activity-type", "IC50"])
    assert args.activity_type == "IC50"


def test_curate_no_longer_accepts_a_flag_it_ignores():
    """--generic-scaffolds affects splitting, not curation, and curate silently ignored it."""
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["curate", "--generic-scaffolds"])


def test_missing_optional_dependency_names_the_install_command(monkeypatch):
    """A missing backend is an install problem and should read as one, not as a traceback."""
    import chembench.cli as cli

    def missing(_args):
        raise ModuleNotFoundError("No module named 'rdkit'", name="rdkit")

    monkeypatch.setattr(cli, "_dispatch", missing)
    with pytest.raises(SystemExit) as excinfo:
        main(["curate"])
    message = str(excinfo.value)
    assert "rdkit" in message
    assert 'pip install -e ".[chem,models]"' in message


def test_curation_summary_totals_the_targets():
    from chembench.curate import CurationReport, curation_summary

    first = CurationReport(target_id="CHEMBL1", fetched=10, kept=6, duplicates_collapsed=3)
    first.reject("censored measurement (relation is not '=')")
    second = CurationReport(target_id="CHEMBL2", fetched=5, kept=4, duplicates_collapsed=0)
    second.reject("censored measurement (relation is not '=')")
    summary = curation_summary([first, second], {"chembl_db_version": "ChEMBL_37"})
    assert summary["chembl_db_version"] == "ChEMBL_37"
    assert summary["totals"] == {
        "fetched": 15,
        "kept": 10,
        "duplicates_collapsed": 3,
        "rejected": {"censored measurement (relation is not '=')": 2},
    }
    assert [t["target_id"] for t in summary["targets"]] == ["CHEMBL1", "CHEMBL2"]
