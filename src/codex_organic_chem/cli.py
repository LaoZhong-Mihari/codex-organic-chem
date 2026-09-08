from __future__ import annotations

import argparse
import json
import sys

from .external import tool_statuses
from .service import (
    COMPUTE_TASKS,
    chem_compute,
    chem_draw,
    chem_figure_tool_status,
    chem_input_review,
    chem_literature_search,
    chem_mechanism_draft,
    chem_mechanism_render,
    chem_mechanism_spec_example,
    chem_mechanism_validate,
    chem_normalize_structure,
    chem_ocsr_benchmark,
    chem_parse_image,
    chem_parse_scheme,
    chem_reaction_analyze,
    chem_route_figure,
    chem_route_figure_spec_example,
    chem_synthesis_suggest,
    chem_structure_review_batch,
    chem_structure_review_result,
    chem_tool_doctor,
)
from .scheme_ocsr import benchmark_to_tsv
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

    p = sub.add_parser("ocsr-benchmark", help="Benchmark OCSR adapters against golden crop SMILES")
    p.add_argument("--gold-smiles", required=True, help="Python/TSV/CSV file containing gold SMILES")
    p.add_argument("--image-dir", required=True, help="Directory containing crop images")
    p.add_argument("--out", choices=["json", "tsv"], default="json", help="Output format")

    p = sub.add_parser("parse-scheme", help="Parse a reaction scheme image and molecule crops with context resolution")
    p.add_argument("--image", required=True, help="Original full scheme image")
    p.add_argument("--crops", required=True, help="Directory containing molecule crop images")
    p.add_argument("--gold-map", help="Optional JSON/TSV compound-to-placeholder definitions")
    p.add_argument("--out", choices=["json"], default="json", help="Output format")

    p = sub.add_parser("input-review", help="Render recognized input and optionally wait for editor confirmation")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--reaction-smiles")
    group.add_argument("--molfile")
    group.add_argument("--image-path")
    p.add_argument("--kind", choices=["auto", "molecule", "reaction"], default="auto")
    p.add_argument("--wait", action="store_true", help="Keep the local review editor live and wait for completion")
    p.add_argument("--timeout-s", type=int, default=1800, help="Review session timeout in seconds")

    p = sub.add_parser("review-batch", help="Start a local Ketcher batch review session for molecule SMILES")
    p.add_argument("--input", required=True, help="JSON file containing a list of review items or an object with an items list")
    p.add_argument("--wait", action="store_true", help="Wait until the browser review session completes or times out")
    p.add_argument("--timeout-s", type=int, default=1800, help="Review session timeout in seconds")

    p = sub.add_parser("review-result", help="Fetch or wait for a local Ketcher review session result")
    p.add_argument("--session-id", required=True)
    p.add_argument("--review-token")
    p.add_argument("--wait", action="store_true")
    p.add_argument("--timeout-s", type=int, default=1800)

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

    p = sub.add_parser("compute", help="Run RDKit descriptors or GFN2-xTB/CREST calculations")
    p.add_argument("--smiles", required=True)
    p.add_argument("--task", action="append", choices=list(COMPUTE_TASKS))
    p.add_argument("--num-confs", type=int, default=8)
    p.add_argument("--max-iters", type=int, default=200)
    p.add_argument("--solvent", help="Implicit solvent for xTB/CREST (ALPB name, e.g. water, thf, dmso)")
    p.add_argument("--timeout", type=int, dest="timeout_s", help="Per-calculation timeout in seconds")

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

    p = sub.add_parser("mechanism-validate", help="Check mapped states, electron moves, and graph edits without rendering")
    p.add_argument("--spec", required=True)

    p = sub.add_parser("route-figure", help="Render a publication-style route/scheme figure from a JSON spec")
    p.add_argument("--spec", help="Path to route figure spec JSON. Omit with --example to print an example spec.")
    p.add_argument("--output-dir")
    p.add_argument("--basename", default="route_figure")
    p.add_argument("--format", action="append", choices=["svg", "png"], dest="formats")
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
    elif args.command == "ocsr-benchmark":
        payload = chem_ocsr_benchmark(args.gold_smiles, args.image_dir)
        if args.out == "tsv":
            sys.stdout.write(benchmark_to_tsv(payload))
            return
    elif args.command == "parse-scheme":
        payload = chem_parse_scheme(args.image, args.crops, gold_map=args.gold_map)
    elif args.command == "input-review":
        payload = chem_input_review(
            smiles=args.smiles,
            reaction_smiles=args.reaction_smiles,
            molfile=args.molfile,
            image_path=args.image_path,
            kind=args.kind,
            interactive=args.wait,
            timeout_s=args.timeout_s,
        )
        if args.wait and payload.get("session_id"):
            sys.stderr.write(f"Review URL: {payload['review_url']}\n")
            sys.stderr.flush()
            payload = chem_structure_review_result(
                session_id=payload["session_id"],
                review_token=payload.get("review_token"),
                wait=True,
                timeout_s=args.timeout_s,
            )
        elif not args.wait and payload.get("status") != "error":
            payload.setdefault("warnings", []).append(
                "Static preview only: use --wait to keep the local review editor live until confirmation."
            )
    elif args.command == "review-batch":
        with open(args.input, encoding="utf-8") as handle:
            review_payload = json.load(handle)
        items = review_payload.get("items") if isinstance(review_payload, dict) else review_payload
        payload = chem_structure_review_batch(items=items, wait=False, timeout_s=args.timeout_s)
        if args.wait:
            sys.stderr.write(f"Review URL: {payload['review_url']}\n")
            sys.stderr.flush()
            payload = chem_structure_review_result(
                session_id=payload["session_id"],
                review_token=payload.get("review_token"),
                wait=True,
                timeout_s=args.timeout_s,
            )
        else:
            payload.setdefault("warnings", []).append(
                "CLI review server exits with this command; use --wait or call through the long-running MCP server."
            )
    elif args.command == "review-result":
        payload = chem_structure_review_result(
            session_id=args.session_id,
            review_token=args.review_token,
            wait=args.wait,
            timeout_s=args.timeout_s,
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
        payload = chem_compute(
            args.smiles,
            tasks=args.task,
            num_confs=args.num_confs,
            max_iters=args.max_iters,
            solvent=args.solvent,
            timeout_s=args.timeout_s,
        )
    elif args.command == "reaction-analyze":
        payload = chem_reaction_analyze(args.reaction, mode=args.mode)
    elif args.command == "mechanism-draft":
        payload = chem_mechanism_draft(
            args.reaction,
            style=args.style,
            quality=args.quality,
            structure_confirmed=args.confirmed,
        )
    elif args.command == "mechanism-validate":
        with open(args.spec, encoding="utf-8") as handle:
            payload = chem_mechanism_validate(json.load(handle))
    elif args.command == "mechanism-render":
        if args.example:
            payload = chem_mechanism_spec_example()
        elif not args.spec:
            payload = {"status": "error", "warnings": ["Provide --spec path or use --example."]}
        else:
            with open(args.spec, encoding="utf-8") as handle:
                spec = json.load(handle)
            payload = chem_mechanism_render(spec=spec, output_dir=args.output_dir)
    elif args.command == "route-figure":
        if args.example:
            payload = chem_route_figure_spec_example()
        elif not args.spec:
            payload = {"status": "error", "warnings": ["Provide --spec path or use --example."]}
        else:
            with open(args.spec, encoding="utf-8") as handle:
                spec = json.load(handle)
            payload = chem_route_figure(
                spec=spec,
                output_dir=args.output_dir,
                basename=args.basename,
                formats=args.formats,
            )
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
