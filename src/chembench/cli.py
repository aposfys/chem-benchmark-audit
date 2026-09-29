"""Command line entry point: ``python -m chembench.cli`` or ``chembench``."""

from __future__ import annotations

import argparse
from pathlib import Path

from chembench import __version__
from chembench.models import MODEL_NAMES

SPLIT_NAMES = ("random", "scaffold", "activity_cliff")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chembench",
        description="Leakage-aware evaluation of molecular property prediction",
    )
    parser.add_argument("--version", action="version", version=f"chembench {__version__}")
    parser.add_argument(
        "--data-dir", type=Path, default=Path("data"), help="cache for curated ChEMBL sets"
    )
    parser.add_argument(
        "--results-dir", type=Path, default=Path("results"), help="where findings are written"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    curate = sub.add_parser(
        "curate",
        help="fetch and curate ChEMBL targets",
        description="Fetch and curate ChEMBL Ki data. Curating the default panel also "
        "writes RESULTS_DIR/curation_report.json with the counts and the ChEMBL release.",
    )
    curate.add_argument("--targets", nargs="*", help="ChEMBL target ids; default is the panel")
    curate.add_argument(
        "--activity-type",
        default="Ki",
        help="ChEMBL standard_type to fetch for --targets (default Ki)",
    )

    evaluate = sub.add_parser("evaluate", help="score every model under every split regime")
    evaluate.add_argument(
        "--models", nargs="*", choices=MODEL_NAMES, default=list(MODEL_NAMES)
    )
    evaluate.add_argument(
        "--splits", nargs="*", choices=SPLIT_NAMES, default=list(SPLIT_NAMES)
    )
    evaluate.add_argument("--targets", type=int, default=0, help="limit to N targets; 0 = all")
    evaluate.add_argument("--seed", type=int, default=0)
    evaluate.add_argument(
        "--generic-scaffolds",
        action="store_true",
        help="erase atom types when computing Murcko scaffolds (a harder split)",
    )
    evaluate.add_argument(
        "--report-only",
        action="store_true",
        help="re-render RESULTS.md from an existing findings.json",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        return _dispatch(args)
    except RuntimeError as exc:
        # Expected failures -- ChEMBL unreachable, nothing curated yet -- are reported with
        # the command that fixes them rather than as a traceback.
        raise SystemExit(str(exc)) from exc
    except ImportError as exc:
        # A missing optional backend (RDKit for curation, the model libraries for the grid)
        # is an install problem, so name the install command instead of a traceback.
        raise SystemExit(
            f"missing optional dependency: {exc.name or exc}.\n"
            'Install the full toolchain with:  pip install -e ".[chem,models]"'
        ) from exc
    except OSError as exc:
        raise SystemExit(f"could not reach ChEMBL: {exc}") from exc


def _dispatch(args: argparse.Namespace) -> int:
    if args.command == "curate":
        import json

        # Curation needs RDKit. Fail here, before a nine-minute download, not after it.
        import rdkit  # noqa: F401

        from chembench.curate import (
            DEFAULT_TARGETS,
            chembl_release,
            curate_target,
            curation_summary,
            fetch_target,
            write_curated,
        )

        panel = (
            [(tid, tid, args.activity_type) for tid in args.targets]
            if args.targets
            else list(DEFAULT_TARGETS)
        )
        release = chembl_release()
        reports = []
        for target_id, name, activity_type in panel:
            raw = fetch_target(
                target_id,
                args.data_dir / "raw" / f"{target_id}.json",
                activity_type=activity_type,
            )
            records, report = curate_target(target_id, raw, activity_type=activity_type)
            write_curated(records, report, args.data_dir / "curated")
            reports.append(report)
            print(
                f"{target_id} {name}: fetched {report.fetched} -> kept {report.kept} "
                f"(collapsed {report.duplicates_collapsed}, "
                f"rejected {sum(report.rejected.values())})"
            )
        if args.targets:
            # The committed summary describes the whole panel. A partial run must not
            # overwrite it with a subset.
            return 0
        summary_path = args.results_dir / "curation_report.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(curation_summary(reports, release), indent=1) + "\n"
        )
        print(f"wrote {summary_path} ({release['chembl_db_version']})")
        return 0

    if args.command == "evaluate":
        from chembench.experiment import run
        from chembench.report import write

        curated = args.data_dir / "curated"
        if not args.report_only and not sorted(curated.glob("CHEMBL*.json")):
            raise SystemExit(
                f"no curated targets in {curated}.\n"
                "Fetch and curate them first:  chembench curate\n"
                "Or re-render an existing run:  chembench evaluate --report-only"
            )
        if args.report_only:
            # Re-render from an existing findings.json. The flag existed and was ignored,
            # so asking for a report silently re-ran a three-hour grid.
            findings_path = args.results_dir / "findings.json"
            if not findings_path.exists():
                raise SystemExit(
                    f"no findings at {findings_path}. Run 'chembench evaluate' first."
                )
            out = write(findings_path, args.results_dir / "RESULTS.md")
            print(f"wrote {out}")
            return 0

        findings = run(
            args.data_dir / "curated",
            args.results_dir,
            models=args.models,
            splits=args.splits,
            generic_scaffolds=getattr(args, "generic_scaffolds", False),
            seed=args.seed,
            max_targets=args.targets,
        )
        out = write(args.results_dir / "findings.json", args.results_dir / "RESULTS.md")
        print(f"wrote {out} ({len(findings['cells'])} cells)")
        return 0

    raise SystemExit(f"unknown command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
