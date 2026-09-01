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
