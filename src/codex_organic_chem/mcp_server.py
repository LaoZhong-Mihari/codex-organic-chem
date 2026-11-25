from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Literal

from .service import (
    chem_compare_structures as _chem_compare_structures,
    chem_compute as _chem_compute,
    chem_crop_image as _chem_crop_image,
    chem_draw as _chem_draw,
    chem_enumerate as _chem_enumerate,
    chem_figure_tool_status as _chem_figure_tool_status,
    chem_input_review as _chem_input_review,
    chem_literature_search as _chem_literature_search,
    chem_mechanism_draft as _chem_mechanism_draft,
    chem_mechanism_render as _chem_mechanism_render,
    chem_mechanism_spec_example as _chem_mechanism_spec_example,
    chem_mechanism_steps_check as _chem_mechanism_steps_check,
    chem_mechanism_validate as _chem_mechanism_validate,
    chem_model_kit as _chem_model_kit,
    chem_normalize_structure as _chem_normalize_structure,
    chem_ocsr_benchmark as _chem_ocsr_benchmark,
    chem_parse_image as _chem_parse_image,
    chem_parse_scheme as _chem_parse_scheme,
    chem_python as _chem_python,
    chem_reaction_analyze as _chem_reaction_analyze,
    chem_review_mode as _chem_review_mode,
    chem_route_check as _chem_route_check,
    chem_route_figure as _chem_route_figure,
    chem_route_figure_spec_example as _chem_route_figure_spec_example,
    chem_route_review as _chem_route_review,
    chem_structure_review_batch as _chem_structure_review_batch,
    chem_structure_review_result as _chem_structure_review_result,
    chem_synthesis_suggest as _chem_synthesis_suggest,
    chem_tool_doctor as _chem_tool_doctor,
    python_tool_enabled as _python_tool_enabled,
)


SERVER_INSTRUCTIONS = (
    "Structure-level checks for organic chemistry: compare structures and stereocentres, count isomers, "
    "stereoisomers and NMR environments, audit a synthetic route (`chem_route_check`) or a stepwise mechanism "
    "intermediate list (`chem_mechanism_steps_check`), answer 3D questions from conformer ensembles (`chem_model_kit`), "
    "run xTB/CREST (`chem_compute`), crop/zoom dense page images (`chem_crop_image`), read structures from images "
    "and draw journal-style figures. If the organic-chemistry-assistant skill is available, load it before your "
    "first call to these tools: it says which tool can refute which kind of claim and how to read each result."
)


def server_instructions() -> str | None:
    """Instructions sent to the client at initialization; CODEX_CHEM_MCP_INSTRUCTIONS=off leaves them out."""
    if os.environ.get("CODEX_CHEM_MCP_INSTRUCTIONS", "on").strip().lower() in {"0", "off", "false", "no"}:
        return None
    return SERVER_INSTRUCTIONS


def _compact_mode() -> bool:
    return os.environ.get("CODEX_CHEM_MCP_COMPACT", "").strip().lower() in {"1", "true", "yes", "on"}


def _skill_root() -> Path | None:
    here = Path(__file__).resolve()
    candidates = [
        here.parents[2] / ".agents" / "skills" / "organic-chemistry-assistant",
        Path(sys.prefix) / ".agents" / "skills" / "organic-chemistry-assistant",
    ]
    for candidate in candidates:
        if (candidate / "SKILL.md").is_file():
            return candidate
    return None


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
    compact = _compact_mode()

    skill_dir = _skill_root()
    if skill_dir is not None:

        @mcp.resource("chem://skill/SKILL.md")
        def skill_main_doc() -> str:
            """Main instructions for the organic-chemistry-assistant skill."""
            return (skill_dir / "SKILL.md").read_text(encoding="utf-8")

        refs_dir = skill_dir / "references"
        if refs_dir.is_dir():
            for ref_file in sorted(refs_dir.glob("*.md")):
                ref_name = ref_file.name

                def _make_reader(path: Path) -> Any:
                    def _read_ref() -> str:
                        return path.read_text(encoding="utf-8")

                    _read_ref.__name__ = f"skill_ref_{path.stem.replace('-', '_')}"
                    _read_ref.__doc__ = f"Reference protocol: {path.name}"
                    return _read_ref

                mcp.resource(f"chem://skill/references/{ref_name}")(_make_reader(ref_file))

    @mcp.tool()
    def chem_crop_image(
        image_path: str | None = None,
        boxes: list[dict[str, Any] | list[float]] | str | None = None,
        output_dir: str | None = None,
        grid: bool = False,
        scale: float = 2.0,
        padding_px: int = 12,
        normalized: bool | None = None,
        path: str | None = None,
    ) -> dict[str, Any]:
        """Crop and magnify subregions of a dense exam/journal page image, or draw a coordinate ruler grid.

        Call with `grid=true` (or omit `boxes`) first on multi-structure pages to write `grid_overlay.png`
        with 0.1-step normalized and pixel rulers; then pass `boxes=[{"id": "A", "box": [x1, y1, x2, y2]}]`
        (0..1 fractions or pixels) to get 2x upscaled PNG crops for accurate wedge/hash and ring reading.
        """
        return _chem_crop_image(
            image_path=image_path,
            boxes=boxes,
            output_dir=output_dir,
            grid=grid,
            scale=scale,
            padding_px=padding_px,
            normalized=normalized,
            path=path,
        )

    @mcp.tool()
    def chem_parse_image(
        path: str | None = None,
        kind: Literal["auto", "molecule", "reaction"] = "auto",
        image_path: str | None = None,
        include_drawings: bool = False,
    ) -> dict[str, Any]:
        """OCSR readings of a structure image, ranked; a second reader for your own transcription.

        Give the image as `path` (alias `image_path`). For multi-structure pages, crop individual
        molecules first with `chem_crop_image`. SVG, Molfile and per-atom geometry are omitted unless
        `include_drawings=true`.
        """
        return _chem_parse_image(path=path, kind=kind, image_path=image_path, include_drawings=include_drawings)

    if not compact:

        @mcp.tool()
        def chem_ocsr_benchmark(gold_smiles: str, image_dir: str) -> dict[str, Any]:
            """Offline OCSR evaluation helper comparing crop images against a TSV/CSV/literal gold SMILES table."""
            return _chem_ocsr_benchmark(gold_smiles=gold_smiles, image_dir=image_dir)

    @mcp.tool()
    def chem_parse_scheme(
        image: str,
        crops: str | list[dict[str, Any]],
        gold_map: str | None = None,
        include_drawings: bool = False,
    ) -> dict[str, Any]:
        """Run OCSR across multiple molecule crops from a reaction scheme and resolve R-group legends.

        `crops` may be either a directory path containing crop images OR a list/JSON string of
        `{"id": "A", "box": [x1, y1, x2, y2]}` bounding boxes to crop directly from `image`.
        """
        return _chem_parse_scheme(
            image=image,
            crops=crops,
            gold_map=gold_map,
            include_drawings=include_drawings,
        )

    @mcp.tool()
    def chem_input_review(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        molfile: str | None = None,
        image_path: str | None = None,
        kind: Literal["auto", "molecule", "reaction"] = "auto",
        interactive: bool = True,
        confirmed: bool = False,
        language: Literal["en", "zh"] | None = None,
        note: str | None = None,
        source_image: str | None = None,
    ) -> dict[str, Any]:
        """Human check of one structure, reaction or image, in the browser.

        After reading a drawing yourself: `smiles=<your reading>` with `source_image=<the image>`,
        so the user sees both side by side. `image_path` alone runs OCSR and reviews its top reading.
        Collect the answer with `chem_structure_review_result`. On `'ask_user_first'`, put `user_prompt`
        to the user and call again with `confirmed=true` if they agree.
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
        language: Literal["en", "zh"] | None = None,
        images: list[str] | None = None,
        confirmed: bool = False,
    ) -> dict[str, Any]:
        """Human check of several structures in one review page.

        items: `[{smiles, id?, label?, image_path?, note?, question?, source?, candidates?}]`.
        Returns a session immediately (`wait=false`); collect with `chem_structure_review_result`.
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
        language: Literal["en", "zh"] | None = None,
        confirmed: bool = False,
        route_path: str | None = None,
        images: list[str] | None = None,
    ) -> dict[str, Any]:
        """Human check of a route, step by step, in the browser.

        Same route format as `chem_route_check`. Collect with `chem_structure_review_result`:
        `route_result.route` is the reviewer's version of the route.
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
    def chem_review_mode(mode: Literal["auto", "ask", "off"] | None = None) -> dict[str, Any]:
        """Show, or with `mode='auto'|'ask'|'off'` save, when review pages open."""
        return _chem_review_mode(mode)

    @mcp.tool()
    def chem_structure_review_result(
        session_id: str,
        review_token: str | None = None,
        wait: bool = False,
        timeout_s: int = 1800,
    ) -> dict[str, Any]:
        """Result of a structure or route review session.

        `wait=true` blocks until the user finishes or `timeout_s` passes; a wait that times out
        returns status `'pending'` with `wait_timed_out` and leaves the session open.
        """
        return _chem_structure_review_result(
            session_id=session_id,
            review_token=review_token,
            wait=wait,
            timeout_s=timeout_s,
        )

    if not compact:

        @mcp.tool()
        def chem_normalize_structure(
            smiles: str | None = None,
            molfile: str | None = None,
            include_drawings: bool = False,
        ) -> dict[str, Any]:
            """Parse and canonicalize a single SMILES or Molfile.

            Note: parsing a SMILES you just wrote proves only syntax. Prefer `chem_compare_structures`
            for stereocenter audits / pairwise comparison, or `chem_mechanism_steps_check` to check
            an entire list of mechanism intermediates in one call.
            """
            return _chem_normalize_structure(smiles=smiles, molfile=molfile, include_drawings=include_drawings)

    @mcp.tool()
    def chem_draw(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        output: Literal["svg", "png", "molfile"] = "svg",
        output_file: str | None = None,
    ) -> dict[str, Any]:
        """Render a single molecule or reaction SMILES to fixed-geometry ACS-style SVG, PNG, or Molfile."""
        return _chem_draw(smiles=smiles, reaction_smiles=reaction_smiles, output=output, output_file=output_file)

    @mcp.tool()
    def chem_compute(
        smiles: str,
        tasks: list[
            Literal["descriptors", "conformers", "charges", "xtb_opt", "xtb_reactivity", "xtb_thermo", "crest"]
        ]
        | None = None,
        num_confs: int = 8,
        max_iters: int = 200,
        solvent: str | None = None,
        timeout_s: int | None = None,
    ) -> dict[str, Any]:
        """Run RDKit or GFN2-xTB/CREST calculations and return structure-mapped numbers.

        `tasks`: `descriptors`, `conformers`, `charges`, `xtb_opt`, `xtb_reactivity` (per-atom Fukui
        indices and partial charges keyed to RDKit atom indices), `xtb_thermo`, `crest`.
        """
        return _chem_compute(
            smiles=smiles,
            tasks=list(tasks) if tasks is not None else None,
            num_confs=num_confs,
            max_iters=max_iters,
            solvent=solvent,
            timeout_s=timeout_s,
        )

    @mcp.tool()
    def chem_tool_doctor() -> dict[str, Any]:
        """Report installed/missing Python dependencies and external engines (xTB, CREST, Open Babel, OCSR)."""
        return _chem_tool_doctor()

    if not compact:

        @mcp.tool()
        def chem_figure_tool_status() -> dict[str, Any]:
            """Report local availability of ChemDraw, ChemDoodle, Marvin, Ketcher, and Inkscape integrations."""
            return _chem_figure_tool_status()

    @mcp.tool()
    def chem_literature_search(
        query: str,
        rows: int = 5,
        source: Literal["auto", "openalex", "crossref"] = "auto",
        year_from: int | None = None,
    ) -> dict[str, Any]:
        """Search OpenAlex (with abstracts, citations, open-access links) plus Crossref.

        Use before committing to a selectivity- or catalyst-dependent prediction and for every route
        step that rests on precedent. Abstract-level evidence only; cite DOI and precedent grade.
        """
        return _chem_literature_search(query=query, rows=rows, source=source, year_from=year_from)

    @mcp.tool()
    def chem_reaction_analyze(
        input: str,
        mode: Literal["sanity_check", "mass_charge_balance", "forward", "retro", "conditions"] = "sanity_check",
        include_drawings: bool = False,
    ) -> dict[str, Any]:
        """Check equation atom/charge balance (`sanity_check` / `mass_charge_balance`) and heuristic reaction classes.

        `input` is a reaction SMILES `reactants>reagents>products` or `reactants>>products`.
        For multi-step synthetic routes use `chem_route_check`; for stepwise mechanism intermediates
        use `chem_mechanism_steps_check`.
        """
        return _chem_reaction_analyze(input=input, mode=mode, include_drawings=include_drawings)

    if not compact:

        @mcp.tool()
        def chem_mechanism_draft(
            reaction: str,
            style: Literal["stepwise", "teaching", "research_note"] = "stepwise",
            quality: Literal["draft", "publication"] = "draft",
            structure_confirmed: bool = False,
        ) -> dict[str, Any]:
            """Rule-template hypothesis sketch for 5 textbook reaction classes (not a general mechanism solver).

            Prefer deriving explicit intermediate SMILES yourself and checking them with
            `chem_mechanism_steps_check` (or `chem_mechanism_validate` for curved-arrow figures).
            """
            return _chem_mechanism_draft(
                reaction=reaction,
                style=style,
                quality=quality,
                structure_confirmed=structure_confirmed,
            )

    @mcp.tool()
    def chem_mechanism_steps_check(
        steps: list[str | dict[str, Any]] | None = None,
        reactants: str | list[str] | None = None,
        products: str | list[str] | None = None,
        conditions: str | None = None,
        reaction: str | None = None,
        steps_path: str | None = None,
        issues_only: bool = False,
    ) -> dict[str, Any]:
        """Audit a proposed list of mechanism intermediates step by step (no curved-arrow spec needed).

        Pass `steps` as SMILES strings or `[{"step": 1, "type": "...", "subtype": "...", "intermediate_smiles": "..."}]`
        together with `reactants`, `products`, and `conditions` (or `reaction="reactants>reagents>products"`).
        Checks RDKit valence, odd-electron radical parity, proton-transfer skeleton invariance, element provenance,
        step-to-step atom/charge continuity, and target product closure (flagging constitutional isomers when a
        cascade/rearrangement step is still missing).
        """
        return _chem_mechanism_steps_check(
            steps=steps,
            reactants=reactants,
            products=products,
            conditions=conditions,
            reaction=reaction,
            steps_path=steps_path,
            issues_only=issues_only,
        )

    @mcp.tool()
    def chem_mechanism_render(spec: dict[str, Any], output_dir: str | None = None) -> dict[str, Any]:
        """Render a validated atom-mapped MechanismSpec 2.1 into curved-arrow SVG, PNG, CDXML, and trace JSON."""
        return _chem_mechanism_render(spec=spec, output_dir=output_dir)

    @mcp.tool()
    def chem_mechanism_validate(spec: dict[str, Any]) -> dict[str, Any]:
        """Validate a curved-arrow MechanismSpec 2.1 (atom-mapped states, electron pushes, and graph edits).

        For plain intermediate SMILES lists without curved-arrow artwork, use `chem_mechanism_steps_check` instead.
        """
        return _chem_mechanism_validate(spec)

    if not compact:

        @mcp.tool()
        def chem_mechanism_spec_example() -> dict[str, Any]:
            """Return a complete MechanismSpec 2.1 JSON example (SN2 curved-arrow template)."""
            return _chem_mechanism_spec_example()

    @mcp.tool()
    def chem_route_figure(
        spec: dict[str, Any],
        output_dir: str | None = None,
        basename: str = "route_figure",
        formats: list[Literal["svg", "png"]] | None = None,
    ) -> dict[str, Any]:
        """Render a multi-step synthesis scheme from a route figure JSON spec (ACS, RSC, or presentation style)."""
        return _chem_route_figure(
            spec=spec,
            output_dir=output_dir,
            basename=basename,
            formats=list(formats) if formats is not None else None,
        )

    if not compact:

        @mcp.tool()
        def chem_route_figure_spec_example() -> dict[str, Any]:
            """Return a complete JSON spec example for `chem_route_figure`."""
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
            """Single-step retrosynthetic disconnection brainstormer for interactive exploration.

            For full route planning, commit to a complete route and audit it with `chem_route_check`.
            """
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
        mcs_timeout_s: int = 6,
        issues_only: bool = False,
    ) -> dict[str, Any]:
        """Audit a synthetic route step by step before presenting it.

        route = `{"target": SMILES, "compounds": {"1": SMILES | {"smiles", "name"}, ...},
        "steps": [{"step": 1, "reactants": ["1", "2"], "product": "3",
        "conditions": "NaH, THF, 0 C", "reagents": [SMILES]?, "byproducts": [SMILES]?}]}`.
        Checks formula delta against conditions, mapped stereocentre retention/inversion, bystander
        group compatibility, and final product vs target. Pass inline (`route`) or via `route_path`.
        """
        return _chem_route_check(
            route=route,
            figure_output_dir=figure_output_dir,
            mcs_timeout_s=mcs_timeout_s,
            route_path=route_path,
            issues_only=issues_only,
        )

    @mcp.tool()
    def chem_compare_structures(
        smiles_a: str | None = None,
        smiles_b: str | None = None,
        a: str | None = None,
        b: str | None = None,
    ) -> dict[str, Any]:
        """Identical, enantiomers, diastereomers (which centres), tautomers, isotopologues, isomers, or formula delta.

        With one SMILES: stereo audit (every stereo element, defined or not, CIP label, meso,
        stereoisomer count). `a` / `b` are accepted as aliases of `smiles_a` / `smiles_b`.
        """
        return _chem_compare_structures(smiles_a=smiles_a, smiles_b=smiles_b, a=a, b=b)

    @mcp.tool()
    def chem_enumerate(
        smiles: str,
        task: Literal["stereoisomers", "substitution", "environments", "groups"] = "stereoisomers",
        substituent: str = "Cl",
        n_substitutions: Literal[1, 2] = 1,
    ) -> dict[str, Any]:
        """Count by enumeration instead of estimating in your head.

        `task`: `stereoisomers` (enantiomer pairs, meso forms) | `substitution` (distinct mono/di-substitution
        products, with and without stereoisomers) | `environments` (symmetry-distinct 13C and 1H NMR counts,
        including diastereotopic groups) | `groups` (functional/protecting groups and liabilities).
        """
        return _chem_enumerate(
            smiles=smiles, task=task, substituent=substituent, n_substitutions=n_substitutions
        )

    @mcp.tool()
    def chem_model_kit(
        task: Literal["atoms", "cis_trans", "face", "distance", "dihedral", "planarity", "embed", "isomer_energies"],
        smiles: str | list[str],
        atoms: list[int] | None = None,
        by_map: bool = False,
        method: Literal["mmff", "xtb"] = "mmff",
        solvent: str | None = None,
        threshold: float = 4.0,
        refine: Literal["auto", "xtb", "none"] = "auto",
    ) -> dict[str, Any]:
        """3D model-kit questions from a conformer ensemble (sterics/ground state only).

        `task`: `atoms` (index table) | `cis_trans` (`atoms=[i, j]`) | `face` (`atoms=[i]`, Re/Si shielding) |
        `distance` (`atoms=[i, j]`) | `dihedral` (`atoms=[a, b, c, d]`) | `planarity` | `embed` |
        `isomer_energies` (`smiles` = list of isomers; `method='mmff'|'xtb'`).
        """
        return _chem_model_kit(
            task=task,
            smiles=smiles,
            atoms=atoms,
            by_map=by_map,
            method=method,
            solvent=solvent,
            threshold=threshold,
            refine=refine,
        )

    if _python_tool_enabled():

        @mcp.tool()
        def chem_python(code: str, timeout_s: int = 120) -> dict[str, Any]:
            """Run a Python snippet with `rdkit`, `numpy`, `PIL`, and `from codex_organic_chem import kit` available.

            Use it to write custom molecular graph enumerations, image inspections, or checks when no
            predefined tool fits (`scipy` and `networkx` are not installed; use `rdkit.Chem` + `numpy`).
            Print what you want to read.
            """
            return _chem_python(code=code, timeout_s=timeout_s)

    mcp.run()


if __name__ == "__main__":
    main()
