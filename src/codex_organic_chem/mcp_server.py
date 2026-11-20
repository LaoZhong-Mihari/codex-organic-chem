from __future__ import annotations

import os
import sys
from typing import Any

from .service import (
    chem_compare_structures as _chem_compare_structures,
    chem_compute as _chem_compute,
    chem_enumerate as _chem_enumerate,
    chem_model_kit as _chem_model_kit,
    chem_python as _chem_python,
    chem_route_check as _chem_route_check,
    python_tool_enabled as _python_tool_enabled,
    chem_draw as _chem_draw,
    chem_figure_tool_status as _chem_figure_tool_status,
    chem_input_review as _chem_input_review,
    chem_review_mode as _chem_review_mode,
    chem_route_review as _chem_route_review,
    chem_literature_search as _chem_literature_search,
    chem_mechanism_draft as _chem_mechanism_draft,
    chem_mechanism_render as _chem_mechanism_render,
    chem_mechanism_spec_example as _chem_mechanism_spec_example,
    chem_mechanism_validate as _chem_mechanism_validate,
    chem_normalize_structure as _chem_normalize_structure,
    chem_ocsr_benchmark as _chem_ocsr_benchmark,
    chem_parse_image as _chem_parse_image,
    chem_parse_scheme as _chem_parse_scheme,
    chem_reaction_analyze as _chem_reaction_analyze,
    chem_route_figure as _chem_route_figure,
    chem_route_figure_spec_example as _chem_route_figure_spec_example,
    chem_synthesis_suggest as _chem_synthesis_suggest,
    chem_structure_review_batch as _chem_structure_review_batch,
    chem_structure_review_result as _chem_structure_review_result,
    chem_tool_doctor as _chem_tool_doctor,
)


SERVER_INSTRUCTIONS = (
    "Structure-level checks for organic chemistry: normalize and compare structures, count isomers, "
    "stereoisomers and NMR environments, audit a synthetic route step by step, answer 3D questions from "
    "conformer ensembles, run xTB/CREST, validate mechanisms, read structures from images and draw figures. "
    "If the organic-chemistry-assistant skill is available, load it before your first call to these tools: "
    "it says which tool can refute which kind of claim and how to read each result."
)


def server_instructions() -> str | None:
    """Instructions sent to the client at initialization; CODEX_CHEM_MCP_INSTRUCTIONS=off leaves them out."""
    if os.environ.get("CODEX_CHEM_MCP_INSTRUCTIONS", "on").strip().lower() in {"0", "off", "false", "no"}:
        return None
    return SERVER_INSTRUCTIONS


def main() -> None:
    try:
        from mcp.server.fastmcp import FastMCP
    except Exception as exc:  # pragma: no cover - optional dependency
        sys.stderr.write(
            "The MCP dependency is not installed. Run `uv sync --extra mcp` "
            f"or install `mcp`. Import error: {exc}\n"
        )
        raise SystemExit(1)

    mcp = FastMCP("codex-organic-chem", instructions=server_instructions())

    @mcp.tool()
    def chem_parse_image(
        path: str | None = None,
        kind: str = "auto",
        image_path: str | None = None,
        include_drawings: bool = False,
    ) -> dict[str, Any]:
        """OCSR readings of a structure image, ranked; a second reader for your own transcription.

        Give the image as `path` (alias `image_path`). Each candidate carries
        SMILES, confidence and consensus metadata; SVG, Molfile and per-atom
        geometry are omitted unless include_drawings=true.
        """
        return _chem_parse_image(path=path, kind=kind, image_path=image_path, include_drawings=include_drawings)

    @mcp.tool()
    def chem_ocsr_benchmark(gold_smiles: str, image_dir: str) -> dict[str, Any]:
        return _chem_ocsr_benchmark(gold_smiles=gold_smiles, image_dir=image_dir)

    @mcp.tool()
    def chem_parse_scheme(image: str, crops: str, gold_map: str | None = None) -> dict[str, Any]:
        return _chem_parse_scheme(image=image, crops=crops, gold_map=gold_map)

    @mcp.tool()
    def chem_input_review(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        molfile: str | None = None,
        image_path: str | None = None,
        kind: str = "auto",
        interactive: bool = True,
        confirmed: bool = False,
        language: str | None = None,
        note: str | None = None,
        source_image: str | None = None,
    ) -> dict[str, Any]:
        """Human check of one structure, reaction or image, in the browser.

        After reading a drawing yourself: smiles=<your reading> with
        source_image=<the image>, so the user sees both side by side.
        image_path alone runs OCSR and reviews its top reading instead.

        Opens the review page (source image beside the structure, click a
        stereocentre to flip it, OCSR alternatives, SMILES or drawing editor)
        and returns a session; collect the answer with
        chem_structure_review_result. The user's review mode decides whether
        the page opens: on 'ask_user_first' put user_prompt to the user and call
        again with confirmed=true if they agree. `note` tells the reviewer how
        you read the structure; `language` is 'en' or 'zh'.
        """
        return _chem_input_review(
            smiles=smiles,
            reaction_smiles=reaction_smiles,
            molfile=molfile,
            image_path=image_path,
            kind=kind,
            interactive=interactive,
            confirmed=confirmed,
            language=language,
            note=note,
            source_image=source_image,
        )

    @mcp.tool()
    def chem_structure_review_batch(
        items: list[dict[str, Any]],
        wait: bool = False,
        timeout_s: int = 1800,
        title: str | None = None,
        intro: str | None = None,
        language: str | None = None,
        images: list[str] | None = None,
        confirmed: bool = False,
    ) -> dict[str, Any]:
        """Human check of several structures in one review page.

        items: [{smiles, id?, label?, image_path? (crop it was read from),
        note? (how you read it, what you assumed), question? (one thing to ask),
        source?, candidates? (alternative readings)}]. `intro` says what to look
        at; `images` are whole source pages. Returns a session at once
        (wait=false); collect with chem_structure_review_result. Statuses
        'ask_user_first' and 'review_disabled' come from the user's review mode.
        """
        return _chem_structure_review_batch(
            items=items,
            wait=wait,
            timeout_s=timeout_s,
            title=title,
            intro=intro,
            language=language,
            images=images,
            confirmed=confirmed,
        )

    @mcp.tool()
    def chem_route_review(
        route: dict[str, Any] | None = None,
        wait: bool = False,
        timeout_s: int = 3600,
        intro: str | None = None,
        language: str | None = None,
        confirmed: bool = False,
        route_path: str | None = None,
        images: list[str] | None = None,
    ) -> dict[str, Any]:
        """Human check of a route, step by step, in the browser.

        Same route format as chem_route_check. The page shows the audit for
        every step; the reviewer can correct any structure or the conditions
        (the audit re-runs), answer audit prompts, comment on steps, and approve
        the route or send it back. Collect with chem_structure_review_result:
        route_result.route is the reviewer's version of the route. Pass the
        route inline or as `route_path`; `images` are source images (the
        drawing the target was read from) shown beside the target.
        """
        return _chem_route_review(
            route=route,
            wait=wait,
            timeout_s=timeout_s,
            intro=intro,
            language=language,
            confirmed=confirmed,
            route_path=route_path,
            images=images,
        )

    @mcp.tool()
    def chem_review_mode(mode: str | None = None) -> dict[str, Any]:
        """Show, or with mode='auto'|'ask'|'off' save, when review pages open."""
        return _chem_review_mode(mode)

    @mcp.tool()
    def chem_structure_review_result(
        session_id: str,
        review_token: str | None = None,
        wait: bool = False,
        timeout_s: int = 1800,
    ) -> dict[str, Any]:
        """Result of a structure or route review session.

        wait=true blocks until the user finishes or timeout_s passes; a wait
        that times out returns status 'pending' with wait_timed_out and leaves
        the session open. Use the returned structures or route, never your
        pre-review copies.
        """
        return _chem_structure_review_result(
            session_id=session_id,
            review_token=review_token,
            wait=wait,
            timeout_s=timeout_s,
        )

    @mcp.tool()
    def chem_normalize_structure(smiles: str | None = None, molfile: str | None = None) -> dict[str, Any]:
        return _chem_normalize_structure(smiles=smiles, molfile=molfile)

    @mcp.tool()
    def chem_draw(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        output: str = "svg",
        output_file: str | None = None,
    ) -> dict[str, Any]:
        return _chem_draw(smiles=smiles, reaction_smiles=reaction_smiles, output=output, output_file=output_file)

    @mcp.tool()
    def chem_compute(
        smiles: str,
        tasks: list[str] | None = None,
        num_confs: int = 8,
        max_iters: int = 200,
        solvent: str | None = None,
        timeout_s: int | None = None,
    ) -> dict[str, Any]:
        """Run RDKit or GFN2-xTB/CREST calculations and return structure-mapped numbers.

        tasks: descriptors, conformers, charges, xtb_opt, xtb_reactivity,
        xtb_thermo, crest. xtb_reactivity returns per-atom Fukui indices and
        partial charges keyed to RDKit atom indices; solvent is an ALPB name.
        """
        return _chem_compute(
            smiles=smiles,
            tasks=tasks,
            num_confs=num_confs,
            max_iters=max_iters,
            solvent=solvent,
            timeout_s=timeout_s,
        )

    @mcp.tool()
    def chem_tool_doctor() -> dict[str, Any]:
        return _chem_tool_doctor()

    @mcp.tool()
    def chem_figure_tool_status() -> dict[str, Any]:
        return _chem_figure_tool_status()

    @mcp.tool()
    def chem_literature_search(
        query: str, rows: int = 5, source: str = "auto", year_from: int | None = None
    ) -> dict[str, Any]:
        """Search OpenAlex (with abstracts, citations, open-access links) plus Crossref.

        Use before committing to a selectivity- or catalyst-dependent prediction and
        for every route step that rests on precedent: search the exact catalyst/ligand
        or reagent with the substrate class, then the reaction class. Abstract-level
        evidence only; cite DOI and say whether it is exact precedent or an analogue.
        """
        return _chem_literature_search(query=query, rows=rows, source=source, year_from=year_from)

    @mcp.tool()
    def chem_reaction_analyze(input: str, mode: str = "sanity_check") -> dict[str, Any]:
        return _chem_reaction_analyze(input=input, mode=mode)

    @mcp.tool()
    def chem_mechanism_draft(
        reaction: str,
        style: str = "stepwise",
        quality: str = "draft",
        structure_confirmed: bool = False,
    ) -> dict[str, Any]:
        return _chem_mechanism_draft(
            reaction=reaction,
            style=style,
            quality=quality,
            structure_confirmed=structure_confirmed,
        )

    @mcp.tool()
    def chem_mechanism_render(spec: dict[str, Any], output_dir: str | None = None) -> dict[str, Any]:
        return _chem_mechanism_render(spec=spec, output_dir=output_dir)

    @mcp.tool()
    def chem_mechanism_validate(spec: dict[str, Any]) -> dict[str, Any]:
        """Validate state identity, mapped H, electron bookkeeping and graph edits.

        A valid result is internal consistency, not evidence for kinetics or
        selectivity. Inspect semantic_checks separately from figure audits.
        """
        return _chem_mechanism_validate(spec)

    @mcp.tool()
    def chem_mechanism_spec_example() -> dict[str, Any]:
        return _chem_mechanism_spec_example()

    @mcp.tool()
    def chem_route_figure(
        spec: dict[str, Any],
        output_dir: str | None = None,
        basename: str = "route_figure",
        formats: list[str] | None = None,
    ) -> dict[str, Any]:
        return _chem_route_figure(spec=spec, output_dir=output_dir, basename=basename, formats=formats)

    @mcp.tool()
    def chem_route_figure_spec_example() -> dict[str, Any]:
        return _chem_route_figure_spec_example()

    @mcp.tool()
    def chem_synthesis_suggest(
        target_smiles: str,
        stage: str = "first_disconnection",
        selected_option: str | None = None,
        confirmed: bool = False,
        literature: bool = True,
        literature_rows: int = 4,
    ) -> dict[str, Any]:
        return _chem_synthesis_suggest(
            target_smiles=target_smiles,
            stage=stage,
            selected_option=selected_option,
            confirmed=confirmed,
            literature=literature,
            literature_rows=literature_rows,
        )

    @mcp.tool()
    def chem_route_check(
        route: dict[str, Any] | None = None,
        figure_output_dir: str | None = None,
        route_path: str | None = None,
    ) -> dict[str, Any]:
        """Audit a synthetic route step by step before presenting it.

        route = {"target": SMILES, "compounds": {"1": SMILES | {"smiles", "name"}, ...},
        "steps": [{"step": 1, "reactants": ["1", "2"], "product": "3",
        "conditions": "NaH, THF, 0 C", "reagents": [SMILES]?, "byproducts": [SMILES]?}]}.
        Checks per step: formula delta read against the conditions (exact closure
        when reagents/byproducts are given), every stereocentre retained / inverted /
        destroyed / created / undefined under an atom mapping, bystander groups the
        reagents attack, and finally last product vs target centre by centre.
        'blocking' issues are facts; 'review' items are heuristic prompts.
        figure_output_dir also draws every numbered compound as a scheme.
        Pass the route inline (`route`) or as the path of its JSON file (`route_path`).
        """
        return _chem_route_check(route=route, figure_output_dir=figure_output_dir, route_path=route_path)

    @mcp.tool()
    def chem_compare_structures(
        smiles_a: str | None = None,
        smiles_b: str | None = None,
        a: str | None = None,
        b: str | None = None,
    ) -> dict[str, Any]:
        """Identical, enantiomers, diastereomers (which centres), tautomers, isomers or different formula.

        With one SMILES: stereo audit (every centre, defined or not, CIP label, meso,
        how many stereoisomers fit the drawing). Use after transcribing a structure
        from an image and to compare an answer with a reference. `a` / `b` are
        accepted as aliases of `smiles_a` / `smiles_b`.
        """
        return _chem_compare_structures(smiles_a=smiles_a, smiles_b=smiles_b, a=a, b=b)

    @mcp.tool()
    def chem_enumerate(
        smiles: str,
        task: str = "stereoisomers",
        substituent: str = "Cl",
        n_substitutions: int = 1,
    ) -> dict[str, Any]:
        """Count by enumeration instead of estimating.

        task: stereoisomers (pairs enantiomers, counts meso forms) | substitution
        (distinct mono/di-substitution products, with and without stereoisomers) |
        environments (symmetry-distinct C and H: expected NMR signal counts, with
        enantiotopic/diastereotopic H) | groups (functional and protecting groups,
        epimerisable or elimination-prone centres).
        """
        return _chem_enumerate(
            smiles=smiles, task=task, substituent=substituent, n_substitutions=n_substitutions
        )

    @mcp.tool()
    def chem_model_kit(
        task: str,
        smiles: str | list[str],
        atoms: list[int] | None = None,
        by_map: bool = False,
        method: str = "mmff",
        solvent: str | None = None,
        threshold: float = 4.0,
    ) -> dict[str, Any]:
        """3D model-kit questions from a conformer ensemble (sterics/ground state only).

        task: atoms (index table) | cis_trans (are the substituents on two ring
        atoms on the same face; atoms=[i, j]) | face (which face of a trigonal atom
        is more open; atoms=[i]) | distance (can two atoms reach; atoms=[i, j]) | dihedral
        (anti/syn-periplanar populations; atoms=[a, b, c, d]) | planarity (flat or
        not, xTB-refined when available) | embed (can the drawn ring fusion exist) |
        isomer_energies (smiles = list of isomers; method mmff or xtb).
        Atom indices follow the SMILES reading order; by_map uses atom-map numbers.
        """
        return _chem_model_kit(
            task=task,
            smiles=smiles,
            atoms=atoms,
            by_map=by_map,
            method=method,
            solvent=solvent,
            threshold=threshold,
        )

    if _python_tool_enabled():

        @mcp.tool()
        def chem_python(code: str, timeout_s: int = 120) -> dict[str, Any]:
            """Run a Python snippet with RDKit and `from codex_organic_chem import kit` available.

            Use it to write whatever check the question needs (enumerate, count,
            compare, measure) when no fixed tool fits. Print what you want to read.
            """
            return _chem_python(code=code, timeout_s=timeout_s)

    mcp.run()


if __name__ == "__main__":
    main()
