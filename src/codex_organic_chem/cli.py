from __future__ import annotations

import argparse
import json
import sys

from .external import tool_statuses
from .service import (
    chem_compute,
    chem_draw,
    chem_figure_tool_status,
    chem_input_review,
    chem_literature_search,
    chem_mechanism_draft,
    chem_mechanism_render,
    chem_mechanism_spec_example,
    chem_normalize_structure,
    chem_parse_image,
    chem_reaction_analyze,
    chem_synthesis_suggest,
    chem_tool_doctor,
)
from .utils import json_dumps


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="codex-chem", description="Codex organic chemistry assistant CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("tool-status", help="Show availability of local chemistry tools")
    sub.add_parser("doctor", help="Explain installed/missing chemistry tools and install commands")
    sub.add_parser("figure-tools", help="Show ChemDraw/ChemDoodle/Marvin/Ketcher/Inkscape integration status")

    p = sub.add_parser("parse-image", help="Parse molecule or reaction image into machine-readable candidates")
    p.add_argument("path")
    p.add_argument("--kind", choices=["auto", "molecule", "reaction"], default="auto")

    p = sub.add_parser("input-review", help="Render recognized input and block downstream work until user confirmation")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--reaction-smiles")
    group.add_argument("--molfile")
    group.add_argument("--image-path")
    p.add_argument("--kind", choices=["auto", "molecule", "reaction"], default="auto")

    p = sub.add_parser("normalize", help="Normalize SMILES or Molfile into a MoleculeRecord")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--molfile")

    p = sub.add_parser("draw", help="Draw molecule or reaction as SVG/Molfile")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--reaction-smiles")
    p.add_argument("--output", choices=["svg", "png", "molfile"], default="svg")
    p.add_argument("--output-file")

    p = sub.add_parser("compute", help="Run lightweight calculations")
    p.add_argument("--smiles", required=True)
    p.add_argument("--task", action="append", choices=["descriptors", "conformers", "charges", "xtb_opt", "crest"])
    p.add_argument("--num-confs", type=int, default=8)
    p.add_argument("--max-iters", type=int, default=200)

    p = sub.add_parser("reaction-analyze", help="Analyze reaction plausibility, conditions, or retrosynthesis hints")
    p.add_argument("--reaction", required=True)
    p.add_argument("--mode", choices=["forward", "retro", "conditions", "sanity_check"], default="sanity_check")

    p = sub.add_parser("mechanism-draft", help="Draft a rule-based organic mechanism")
    p.add_argument("--reaction", required=True)
    p.add_argument("--style", choices=["stepwise", "teaching", "research_note"], default="stepwise")
    p.add_argument("--quality", choices=["draft", "publication"], default="draft")
    p.add_argument("--confirmed", action="store_true", help="Assert the rendered input has been explicitly user-confirmed")

    p = sub.add_parser("mechanism-render", help="Render an explicit atom-mapped mechanism canvas with intermediates, lone pairs, charges, and arrows")
    p.add_argument("--spec", help="Path to mechanism spec JSON. Omit with --example to print an example spec.")
    p.add_argument("--output-dir")
    p.add_argument("--example", action="store_true")

    p = sub.add_parser("literature-search", help="Search literature metadata for chemistry evidence")
    p.add_argument("--query", required=True)
    p.add_argument("--rows", type=int, default=5)

    p = sub.add_parser("synthesis-suggest", help="Suggest exactly one layer of retrosynthetic options with tool/literature evidence")
    p.add_argument("--target-smiles", required=True)
    p.add_argument("--stage", default="first_disconnection")
    p.add_argument("--selected-option")
    p.add_argument("--confirmed", action="store_true", help="Assert the rendered target has been explicitly user-confirmed")
    p.add_argument("--no-literature", action="store_true")
    p.add_argument("--literature-rows", type=int, default=4)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "tool-status":
        payload = tool_statuses()
    elif args.command == "doctor":
        payload = chem_tool_doctor()
    elif args.command == "figure-tools":
        payload = chem_figure_tool_status()
    elif args.command == "parse-image":
        payload = chem_parse_image(args.path, kind=args.kind)
    elif args.command == "input-review":
        payload = chem_input_review(
            smiles=args.smiles,
            reaction_smiles=args.reaction_smiles,
            molfile=args.molfile,
            image_path=args.image_path,
            kind=args.kind,
        )
    elif args.command == "normalize":
        payload = chem_normalize_structure(smiles=args.smiles, molfile=args.molfile)
    elif args.command == "draw":
        payload = chem_draw(
            smiles=args.smiles,
            reaction_smiles=args.reaction_smiles,
            output=args.output,
            output_file=args.output_file,
        )
    elif args.command == "compute":
        payload = chem_compute(args.smiles, tasks=args.task, num_confs=args.num_confs, max_iters=args.max_iters)
    elif args.command == "reaction-analyze":
        payload = chem_reaction_analyze(args.reaction, mode=args.mode)
    elif args.command == "mechanism-draft":
        payload = chem_mechanism_draft(
            args.reaction,
            style=args.style,
            quality=args.quality,
            structure_confirmed=args.confirmed,
        )
    elif args.command == "mechanism-render":
        if args.example:
            payload = chem_mechanism_spec_example()
        elif not args.spec:
            payload = {"status": "error", "warnings": ["Provide --spec path or use --example."]}
        else:
            with open(args.spec, encoding="utf-8") as handle:
                spec = json.load(handle)
            payload = chem_mechanism_render(spec=spec, output_dir=args.output_dir)
    elif args.command == "literature-search":
        payload = chem_literature_search(args.query, rows=args.rows)
    elif args.command == "synthesis-suggest":
        payload = chem_synthesis_suggest(
            args.target_smiles,
            stage=args.stage,
            selected_option=args.selected_option,
            confirmed=args.confirmed,
            literature=not args.no_literature,
            literature_rows=args.literature_rows,
        )
    else:  # pragma: no cover
        parser.error(f"Unknown command: {args.command}")
        return
    sys.stdout.write(json_dumps(payload))
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
