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

    curate = sub.add_parser("curate", help="fetch and curate ChEMBL targets")
    curate.add_argument("--targets", nargs="*", help="ChEMBL target ids; default is the panel")
    curate.add_argument(
        "--generic-scaffolds",
        action="store_true",
        help="erase atom types when computing Murcko scaffolds (a harder split)",
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

    if args.command == "curate":
        from chembench.curate import (
            DEFAULT_TARGETS,
            curate_target,
            fetch_target,
            write_curated,
        )

        panel = (
            [(tid, tid, "Ki") for tid in args.targets]
            if args.targets
            else list(DEFAULT_TARGETS)
        )
        for target_id, name, activity_type in panel:
            raw = fetch_target(
                target_id,
                args.data_dir / "raw" / f"{target_id}.json",
                activity_type=activity_type,
            )
            records, report = curate_target(target_id, raw, activity_type=activity_type)
            write_curated(records, report, args.data_dir / "curated")
            print(
                f"{target_id} {name}: fetched {report.fetched} -> kept {report.kept} "
                f"(collapsed {report.duplicates_collapsed}, "
                f"rejected {sum(report.rejected.values())})"
            )
        return 0

    if args.command == "evaluate":
        from chembench.experiment import run
        from chembench.report import write

        if args.report_only:
            # Re-render from an existing findings.json. The flag existed and was ignored,
            # so asking for a report silently re-ran a three-hour grid.
            out = write(args.results_dir / "findings.json", args.results_dir / "RESULTS.md")
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
