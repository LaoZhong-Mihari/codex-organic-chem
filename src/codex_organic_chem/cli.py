from __future__ import annotations

import argparse
import json
import sys

from .external import tool_statuses
from .service import (
    COMPUTE_TASKS,
    ENUMERATE_TASKS,
    chem_compare_structures,
    chem_compute,
    chem_crop_image,
    chem_draw,
    chem_enumerate,
    chem_figure_tool_status,
    chem_input_review,
    chem_review_demo,
    chem_review_mode,
    chem_route_review,
    chem_literature_search,
    chem_mechanism_draft,
    chem_mechanism_render,
    chem_mechanism_spec_example,
    chem_mechanism_steps_check,
    chem_mechanism_validate,
    chem_model_kit,
    chem_normalize_structure,
    chem_ocsr_benchmark,
    chem_parse_image,
    chem_parse_scheme,
    chem_reaction_analyze,
    chem_route_check,
    chem_route_check_example,
    chem_route_figure,
    chem_route_figure_spec_example,
    chem_synthesis_suggest,
    chem_structure_review_batch,
    chem_structure_review_result,
    chem_tool_doctor,
)
from .model_kit import TASKS as MODEL_KIT_TASKS
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
    p.add_argument("--crops", required=True, help="Directory of crop images or JSON list of {id, box: [x1, y1, x2, y2]}")
    p.add_argument("--gold-map", help="Optional JSON/TSV compound-to-placeholder definitions")
    p.add_argument("--out", choices=["json"], default="json", help="Output format")

    p = sub.add_parser("crop-image", help="Crop and magnify subregions of a page/scheme image or overlay a coordinate grid")
    p.add_argument("image_path", help="Source page or scheme image")
    p.add_argument("--boxes", help="JSON list of {id, box: [x1, y1, x2, y2]} (0..1 normalized or pixels)")
    p.add_argument("--output-dir", help="Directory to write cropped PNGs / grid image")
    p.add_argument("--grid", action="store_true", help="Also write grid_overlay.png with 0.1-step coordinate rulers")
    p.add_argument("--scale", type=float, default=2.0, help="Upscale factor for crops (default: 2.0)")
    p.add_argument("--padding-px", type=int, default=12, help="Padding pixels around each crop box")

    def review_options(parser: argparse.ArgumentParser, lifetime: int = 1800) -> None:
        parser.add_argument("--wait", action="store_true", help="Block until the user finishes the review, for at most --wait-s")
        parser.add_argument("--wait-s", type=int, default=540, help="Longest time --wait blocks; the session stays open after it")
        parser.add_argument("--timeout-s", type=int, default=lifetime, help="How long the review session stays open, in seconds")
        parser.add_argument("--yes", action="store_true", help="The user asked for the page or agreed to it (needed when review mode is ask or off)")
        parser.add_argument("--lang", choices=["en", "zh"], help="Page language (default: the browser's)")

    p = sub.add_parser("input-review", help="Open the review page on one structure, reaction or image and return a session")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles")
    group.add_argument("--reaction-smiles")
    group.add_argument("--molfile")
    group.add_argument("--image-path")
    p.add_argument("--kind", choices=["auto", "molecule", "reaction"], default="auto")
    p.add_argument("--source-image", help="With --smiles: the image the structure was read from, shown beside it")
    p.add_argument("--note", help="How the structure was read and what was assumed; shown to the reviewer")
    p.add_argument("--static", action="store_true", help="Only return the rendered preview; do not open a review page")
    review_options(p)

    p = sub.add_parser("review-batch", help="Open the review page on several structures (source crops, notes, alternatives)")
    p.add_argument("--input", required=True, help="JSON: a list of items, or {items, title, intro, images, language}")
    p.add_argument("--title")
    p.add_argument("--intro", help="One or two sentences telling the reviewer what to look at")
    p.add_argument("--image", action="append", help="Whole source page shown above the queue (repeatable)")
    review_options(p)

    p = sub.add_parser("route-review", help="Open the review page on a route: audit per step, editable structures and conditions")
    p.add_argument("--route", required=True, help="Route JSON (same format as route-check)")
    p.add_argument("--intro", help="One or two sentences telling the reviewer what to look at")
    p.add_argument("--image", action="append", help="Source image the target was read from, shown beside it (repeatable)")
    review_options(p, lifetime=3600)

    p = sub.add_parser("review-result", help="Fetch, or with --wait block for, the result of a review session")
    p.add_argument("--session-id", required=True)
    p.add_argument("--review-token")
    p.add_argument("--wait", action="store_true")
    p.add_argument("--timeout-s", type=int, default=540, help="How long to wait; a timeout leaves the session open")

    p = sub.add_parser("review-mode", help="Show or set when the review page opens: auto, ask or off")
    p.add_argument("mode", nargs="?", choices=["auto", "ask", "off"])

    p = sub.add_parser("review-demo", help="Open two practice review pages with planted mistakes")
    p.add_argument("--lang", choices=["en", "zh"], default="en")
    p.add_argument("--which", choices=["both", "structures", "route"], default="both")

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
    p.add_argument(
        "--mode",
        choices=["forward", "retro", "conditions", "sanity_check", "mass_charge_balance"],
        default="sanity_check",
    )

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

    p = sub.add_parser("mechanism-steps-check", help="Audit an intermediate SMILES sequence step by step (valence, skeleton, delta, product closure)")
    p.add_argument("--steps", required=True, help="Path to JSON file or inline JSON list of intermediate SMILES / step dicts")
    p.add_argument("--reactants", help="Dot-separated reactant SMILES")
    p.add_argument("--products", help="Dot-separated target product SMILES")
    p.add_argument("--conditions", help="Free-text reagents and conditions")
    p.add_argument("--reaction", help="Reaction SMILES shorthand (reactants>reagents>products)")
    p.add_argument("--issues-only", action="store_true", help="Omit per-step details and report status + issues only")

    p = sub.add_parser("route-figure", help="Render a publication-style route/scheme figure from a JSON spec")
    p.add_argument("--spec", help="Path to route figure spec JSON. Omit with --example to print an example spec.")
    p.add_argument("--output-dir")
    p.add_argument("--basename", default="route_figure")
    p.add_argument("--format", action="append", choices=["svg", "png"], dest="formats")
    p.add_argument("--example", action="store_true")

    p = sub.add_parser("route-check", help="Audit a route step by step: atoms, stereocentres, groups, target identity")
    p.add_argument("--route", help="Path to route JSON. Omit with --example to print the format.")
    p.add_argument("--figure-dir", help="Also draw every numbered compound as a scheme into this directory")
    p.add_argument("--mcs-timeout", type=int, default=6, help="Seconds allowed per atom-mapping search")
    p.add_argument("--example", action="store_true")
    p.add_argument("--issues-only", action="store_true", help="Print status and issues without per-step detail")

    p = sub.add_parser("compare", help="Relationship between two structures (or stereo audit of one)")
    p.add_argument("--a", required=True, help="SMILES")
    p.add_argument("--b", help="SMILES; omit for a stereo audit of --a")

    p = sub.add_parser("enumerate", help="Count stereoisomers, substitution products, NMR environments or groups")
    p.add_argument("--smiles", required=True)
    p.add_argument("--task", choices=list(ENUMERATE_TASKS), default="stereoisomers")
    p.add_argument("--substituent", default="Cl")
    p.add_argument("--n", type=int, default=1, dest="n_substitutions", choices=[1, 2])

    p = sub.add_parser("model-kit", help="3D questions: cis/trans on a ring, face accessibility, reach, dihedral, planarity, feasibility, isomer energies")
    p.add_argument("--task", choices=list(MODEL_KIT_TASKS), required=True)
    p.add_argument("--smiles", required=True, action="append", help="Repeat for isomer_energies")
    p.add_argument("--atoms", help="Comma-separated atom indices (SMILES order) or map numbers with --by-map")
    p.add_argument("--by-map", action="store_true")
    p.add_argument("--method", choices=["mmff", "xtb"], default="mmff")
    p.add_argument("--solvent")
    p.add_argument("--threshold", type=float, default=4.0)
    p.add_argument("--refine", choices=["auto", "xtb", "none"], default="auto")

    p = sub.add_parser("literature-search", help="Search literature metadata for chemistry evidence")
    p.add_argument("--query", required=True)
    p.add_argument("--rows", type=int, default=5)
    p.add_argument("--source", choices=["auto", "openalex", "crossref"], default="auto")
    p.add_argument("--year-from", type=int)

    p = sub.add_parser("synthesis-suggest", help="Suggest exactly one layer of retrosynthetic options with tool/literature evidence")
    p.add_argument("--target-smiles", required=True)
    p.add_argument("--stage", default="first_disconnection")
    p.add_argument("--selected-option")
    p.add_argument("--confirmed", action="store_true", help="Assert the rendered target has been explicitly user-confirmed")
    p.add_argument("--no-literature", action="store_true")
    p.add_argument("--literature-rows", type=int, default=4)

    return parser


def _follow_review(payload: dict, args: argparse.Namespace) -> dict:
    """Announce an opened review page and, with --wait, block for its result.

    The page is served by the detached review server, so it stays up after
    this command returns; `review-result --session-id ... --wait` collects it.
    """
    if not payload.get("session_id"):
        return payload
    sys.stderr.write(f"Review URL: {payload['review_url']}\n")
    sys.stderr.flush()
    collect = f"codex-chem review-result --session-id {payload['session_id']} --wait"
    if args.wait:
        result = chem_structure_review_result(
            session_id=payload["session_id"], wait=True, timeout_s=min(args.wait_s, args.timeout_s)
        )
        result["browser_opened"] = payload.get("browser_opened", False)
        if result.get("wait_timed_out"):
            result["collect_with"] = collect
        return result
    payload["collect_with"] = f"codex-chem review-result --session-id {payload['session_id']} --wait"
    return payload


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
    elif args.command == "crop-image":
        payload = chem_crop_image(
            image_path=args.image_path,
            boxes=args.boxes,
            output_dir=args.output_dir,
            grid=args.grid,
            scale=args.scale,
            padding_px=args.padding_px,
        )
    elif args.command == "input-review":
        payload = chem_input_review(
            smiles=args.smiles,
            reaction_smiles=args.reaction_smiles,
            molfile=args.molfile,
            image_path=args.image_path,
            kind=args.kind,
            interactive=not args.static,
            timeout_s=args.timeout_s,
            confirmed=args.yes,
            detached=True,
            language=args.lang,
            note=args.note,
            source_image=args.source_image,
        )
        payload = _follow_review(payload, args)
    elif args.command == "review-batch":
        with open(args.input, encoding="utf-8") as handle:
            review_payload = json.load(handle)
        options = review_payload if isinstance(review_payload, dict) else {}
        items = options.get("items") if isinstance(review_payload, dict) else review_payload
        payload = chem_structure_review_batch(
            items=items,
            wait=False,
            timeout_s=args.timeout_s,
            title=args.title or options.get("title"),
            intro=args.intro or options.get("intro"),
            language=args.lang or options.get("language"),
            images=args.image or options.get("images"),
            confirmed=args.yes,
            detached=True,
        )
        payload = _follow_review(payload, args)
    elif args.command == "route-review":
        with open(args.route, encoding="utf-8") as handle:
            route_payload = json.load(handle)
        payload = chem_route_review(
            route=route_payload,
            wait=False,
            timeout_s=args.timeout_s,
            intro=args.intro,
            language=args.lang,
            confirmed=args.yes,
            detached=True,
            images=args.image,
        )
        payload = _follow_review(payload, args)
    elif args.command == "review-result":
        try:
            payload = chem_structure_review_result(
                session_id=args.session_id,
                review_token=args.review_token,
                wait=args.wait,
                timeout_s=args.timeout_s,
            )
        except KeyError:
            payload = {"status": "error", "warnings": [f"Unknown review session: {args.session_id}"]}
    elif args.command == "review-mode":
        payload = chem_review_mode(args.mode)
    elif args.command == "review-demo":
        payload = chem_review_demo(language=args.lang, which=args.which)
        for name, session in payload.get("sessions", {}).items():
            sys.stderr.write(f"Review URL ({name}): {session['review_url']}\n")
        sys.stderr.flush()
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
    elif args.command == "mechanism-steps-check":
        raw_steps = args.steps.strip()
        if raw_steps.startswith(("[", "{")):
            parsed = json.loads(raw_steps)
            if isinstance(parsed, dict):
                payload = chem_mechanism_steps_check(
                    steps=parsed.get("steps") or parsed.get("mechanism"),
                    reactants=args.reactants or parsed.get("reactants") or parsed.get("reactants_smiles"),
                    products=args.products or parsed.get("products") or parsed.get("products_smiles"),
                    conditions=args.conditions or parsed.get("conditions"),
                    reaction=args.reaction or parsed.get("reaction"),
                    issues_only=args.issues_only,
                )
            else:
                payload = chem_mechanism_steps_check(
                    steps=parsed,
                    reactants=args.reactants,
                    products=args.products,
                    conditions=args.conditions,
                    reaction=args.reaction,
                    issues_only=args.issues_only,
                )
        else:
            payload = chem_mechanism_steps_check(
                steps_path=raw_steps,
                reactants=args.reactants,
                products=args.products,
                conditions=args.conditions,
                reaction=args.reaction,
                issues_only=args.issues_only,
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
    elif args.command == "route-check":
        if args.example:
            payload = chem_route_check_example()
        elif not args.route:
            payload = {"status": "error", "warnings": ["Provide --route path or use --example."]}
        else:
            with open(args.route, encoding="utf-8") as handle:
                payload = chem_route_check(
                    json.load(handle), figure_output_dir=args.figure_dir, mcs_timeout_s=args.mcs_timeout
                )
            if args.issues_only:
                payload = {k: v for k, v in payload.items() if k != "steps"}
    elif args.command == "compare":
        payload = chem_compare_structures(args.a, args.b)
    elif args.command == "enumerate":
        payload = chem_enumerate(
            args.smiles, task=args.task, substituent=args.substituent, n_substitutions=args.n_substitutions
        )
    elif args.command == "model-kit":
        atoms = [int(x) for x in args.atoms.split(",")] if args.atoms else None
        payload = chem_model_kit(
            args.task,
            args.smiles if args.task == "isomer_energies" else args.smiles[0],
            atoms=atoms,
            by_map=args.by_map,
            method=args.method,
            solvent=args.solvent,
            threshold=args.threshold,
            refine=args.refine,
        )
    elif args.command == "literature-search":
        payload = chem_literature_search(args.query, rows=args.rows, source=args.source, year_from=args.year_from)
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
