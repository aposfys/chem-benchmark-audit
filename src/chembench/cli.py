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

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    raise SystemExit(f"'{args.command}' is not implemented yet; see README milestones")


if __name__ == "__main__":
    raise SystemExit(main())
