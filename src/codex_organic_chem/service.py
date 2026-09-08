from __future__ import annotations

from pathlib import Path
from typing import Any

from .external import doctor_report, tool_statuses
from .figure_tools import figure_tool_statuses
from .input_review import prepare_input_review
from .literature import literature_search
from .mechanism import draft_mechanism
from .mechanism_canvas import mechanism_spec_example, render_mechanism_canvas, validate_mechanism_spec
from .models import CalculationRecord
from .ocsr import parse_image
from .rdkit_tools import (
    charges_record,
    conformer_record,
    descriptors_record,
    normalize_structure,
    reaction_to_svg,
)
from .reaction import analyze_reaction
from .route_figure import render_route_figure, route_figure_spec_example
from .scheme_ocsr import benchmark_ocsr, parse_scheme
from .semiempirical import (
    crest_conformer_record,
    xtb_opt_record,
    xtb_reactivity_record,
    xtb_thermo_record,
)
from .synthesis import suggest_synthesis_path
from .structure_review import review_session_result, start_structure_review_batch


def chem_parse_image(path: str, kind: str = "auto") -> dict[str, Any]:
    return parse_image(path=path, kind=kind)


def chem_ocsr_benchmark(gold_smiles: str, image_dir: str) -> dict[str, Any]:
    return benchmark_ocsr(gold_smiles=gold_smiles, image_dir=image_dir)


def chem_parse_scheme(image: str, crops: str, gold_map: str | None = None) -> dict[str, Any]:
    return parse_scheme(image=image, crops=crops, gold_map=gold_map)


def chem_input_review(
    smiles: str | None = None,
    reaction_smiles: str | None = None,
    molfile: str | None = None,
    image_path: str | None = None,
    kind: str = "auto",
    interactive: bool = True,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    return prepare_input_review(
        smiles=smiles,
        reaction_smiles=reaction_smiles,
        molfile=molfile,
        image_path=image_path,
        kind=kind,
        interactive=interactive,
        timeout_s=timeout_s,
    )


def chem_structure_review_batch(
    items: list[dict[str, Any]],
    wait: bool = False,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    return start_structure_review_batch(items=items, wait=wait, timeout_s=timeout_s)


def chem_structure_review_result(
    session_id: str,
    review_token: str | None = None,
    wait: bool = False,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    return review_session_result(session_id=session_id, token=review_token, wait=wait, timeout_s=timeout_s)


def chem_normalize_structure(smiles: str | None = None, molfile: str | None = None) -> dict[str, Any]:
    return normalize_structure(smiles=smiles, molfile=molfile, source="user").to_dict()


def chem_draw(
    smiles: str | None = None,
    reaction_smiles: str | None = None,
    output: str = "svg",
    output_file: str | None = None,
) -> dict[str, Any]:
    warnings: list[str] = []
    if reaction_smiles is None and smiles and (">" in smiles):
        reaction_smiles = smiles
        smiles = None
    if output not in {"svg", "png", "molfile"}:
        return {"status": "error", "warnings": [f"Unsupported output format: {output}"]}
    payload: str | None = None
    kind = "molecule"
    if reaction_smiles:
        kind = "reaction"
        if output == "molfile":
            return {"status": "unavailable", "warnings": ["Reaction RXN export is not implemented in MVP."]}
        payload, draw_warnings = reaction_to_svg(reaction_smiles)
        warnings.extend(draw_warnings)
    elif smiles:
        record = normalize_structure(smiles=smiles, source="draw")
        warnings.extend(record.warnings)
        if output == "molfile":
            payload = record.molblock
        else:
            # Draw single molecules through the shared fixed-geometry engine so
            # standalone drawings match route/mechanism figures exactly.
            payload = _standalone_molecule_svg(
                record.isomeric_smiles or record.canonical_smiles or smiles, warnings
            ) or record.svg
    else:
        return {"status": "error", "warnings": ["Provide smiles or reaction_smiles."]}

    if output == "png":
        return _draw_png(kind, payload, output_file, warnings)

    result = {
        "status": "ok" if payload else "error",
        "kind": kind,
        "format": output,
        "data": payload,
        "warnings": warnings,
    }
    if output_file and payload:
        path = Path(output_file).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload, encoding="utf-8")
        result["output_file"] = str(path)
    return result


def _standalone_molecule_svg(smiles: str | None, warnings: list[str]) -> str | None:
    if not smiles:
        return None
    try:
        from .figure_style import render_molecule

        depiction = render_molecule(smiles=smiles, preset="acs")
    except Exception:
        return None
    warnings.extend(depiction.warnings)
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{depiction.width}px' "
        f"height='{depiction.height}px' viewBox='0 0 {depiction.width} {depiction.height}'>\n"
        f"<rect width='100%' height='100%' fill='white'/>\n{depiction.svg_body}\n</svg>"
    )


def _draw_png(
    kind: str, svg_payload: str | None, output_file: str | None, warnings: list[str]
) -> dict[str, Any]:
    """Rasterize the checked SVG master; PNG is always derived, never drawn."""
    import tempfile

    from .figure_audit import raster_ink_check, svg_to_png

    if not svg_payload:
        return {"status": "error", "kind": kind, "format": "png", "warnings": warnings}
    if not output_file:
        return {
            "status": "error",
            "kind": kind,
            "format": "png",
            "warnings": [*warnings, "PNG output requires output_file (binary data is not returned inline)."],
        }
    png_path = Path(output_file).expanduser().resolve()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".svg", mode="w", delete=False, encoding="utf-8") as handle:
        handle.write(svg_payload)
        svg_tmp = Path(handle.name)
    try:
        raster_warnings = svg_to_png(svg_tmp, png_path)
    finally:
        svg_tmp.unlink(missing_ok=True)
    warnings.extend(raster_warnings)
    if not png_path.exists() or raster_warnings:
        return {"status": "unavailable", "kind": kind, "format": "png", "warnings": warnings}
    ink_warnings, ink_summary = raster_ink_check(png_path)
    warnings.extend(ink_warnings)
    return {
        "status": "ok",
        "kind": kind,
        "format": "png",
        "svg_master": svg_payload,
        "output_file": str(png_path),
        "raster_check": ink_summary,
        "warnings": warnings,
    }


COMPUTE_TASKS = (
    "descriptors",
    "conformers",
    "charges",
    "xtb_opt",
    "xtb_reactivity",
    "xtb_thermo",
    "crest",
)


def chem_compute(
    smiles: str,
    tasks: list[str] | None = None,
    num_confs: int = 8,
    max_iters: int = 200,
    solvent: str | None = None,
    timeout_s: int | None = None,
) -> dict[str, Any]:
    selected = tasks or ["descriptors"]
    records: list[CalculationRecord] = []
    for task in selected:
        if task == "descriptors":
            records.append(descriptors_record(smiles))
        elif task == "conformers":
            records.append(conformer_record(smiles, num_confs=num_confs, max_iters=max_iters))
        elif task == "charges":
            records.append(charges_record(smiles))
        elif task == "xtb_opt":
            records.append(xtb_opt_record(smiles, solvent=solvent, timeout_s=timeout_s or 300))
        elif task == "xtb_reactivity":
            records.append(xtb_reactivity_record(smiles, solvent=solvent, timeout_s=timeout_s or 300))
        elif task == "xtb_thermo":
            records.append(xtb_thermo_record(smiles, solvent=solvent, timeout_s=timeout_s or 600))
        elif task == "crest":
            records.append(crest_conformer_record(smiles, solvent=solvent, timeout_s=timeout_s or 600))
        else:
            records.append(
                CalculationRecord(
                    method=task,
                    parameters={"smiles": smiles},
                    status="unavailable",
                    warnings=[f"Unknown calculation task: {task}"],
                )
            )
    return {
        "smiles": smiles,
        "records": [record.to_dict() for record in records],
        "tool_status": tool_statuses(),
    }


def chem_tool_doctor() -> dict[str, Any]:
    return doctor_report()


def chem_figure_tool_status() -> dict[str, Any]:
    return figure_tool_statuses()


def chem_literature_search(query: str, rows: int = 5) -> dict[str, Any]:
    return literature_search(query=query, rows=rows)


def chem_reaction_analyze(input: str, mode: str = "sanity_check") -> dict[str, Any]:
    if mode not in {"forward", "retro", "conditions", "sanity_check"}:
        return {"status": "error", "warnings": [f"Unsupported reaction analysis mode: {mode}"]}
    return analyze_reaction(input, mode=mode).to_dict()


def chem_mechanism_draft(
    reaction: str,
    style: str = "stepwise",
    quality: str = "draft",
    structure_confirmed: bool = False,
) -> dict[str, Any]:
    if style not in {"stepwise", "teaching", "research_note"}:
        return {"status": "error", "warnings": [f"Unsupported mechanism style: {style}"]}
    if quality not in {"draft", "publication"}:
        return {"status": "error", "warnings": [f"Unsupported mechanism quality: {quality}"]}
    return draft_mechanism(reaction, style=style, quality=quality, structure_confirmed=structure_confirmed)


def chem_mechanism_render(spec: dict[str, Any], output_dir: str | None = None) -> dict[str, Any]:
    return render_mechanism_canvas(spec=spec, output_dir=output_dir)


def chem_mechanism_validate(spec: dict[str, Any]) -> dict[str, Any]:
    """Check state/electron consistency without rendering or claiming plausibility."""
    return validate_mechanism_spec(spec)


def chem_mechanism_spec_example() -> dict[str, Any]:
    return mechanism_spec_example()


def chem_route_figure(
    spec: dict[str, Any],
    output_dir: str | None = None,
    basename: str = "route_figure",
    formats: list[str] | None = None,
) -> dict[str, Any]:
    return render_route_figure(spec=spec, output_dir=output_dir, basename=basename, formats=formats)


def chem_route_figure_spec_example() -> dict[str, Any]:
    return route_figure_spec_example()


def chem_synthesis_suggest(
    target_smiles: str,
    stage: str = "first_disconnection",
    selected_option: str | None = None,
    confirmed: bool = False,
    literature: bool = True,
    literature_rows: int = 4,
) -> dict[str, Any]:
    return suggest_synthesis_path(
        target_smiles=target_smiles,
        stage=stage,
        selected_option=selected_option,
        confirmed=confirmed,
        literature=literature,
        literature_rows=literature_rows,
    )
