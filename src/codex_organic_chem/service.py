from __future__ import annotations

from pathlib import Path
from typing import Any

from .external import crest_record, doctor_report, tool_statuses, xtb_opt_record
from .figure_tools import figure_tool_statuses
from .input_review import prepare_input_review
from .literature import literature_search
from .mechanism import draft_mechanism
from .mechanism_canvas import mechanism_spec_example, render_mechanism_canvas
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
from .scheme_ocsr import benchmark_ocsr, parse_scheme
from .synthesis import suggest_synthesis_path


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
) -> dict[str, Any]:
    return prepare_input_review(
        smiles=smiles,
        reaction_smiles=reaction_smiles,
        molfile=molfile,
        image_path=image_path,
        kind=kind,
    )


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
    if output == "png":
        return {
            "status": "unavailable",
            "warnings": ["PNG export is not implemented in MVP; request SVG or Molfile."],
        }
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
        payload = record.svg if output == "svg" else record.molblock
    else:
        return {"status": "error", "warnings": ["Provide smiles or reaction_smiles."]}
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


def chem_compute(
    smiles: str,
    tasks: list[str] | None = None,
    num_confs: int = 8,
    max_iters: int = 200,
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
            records.append(xtb_opt_record(smiles))
        elif task == "crest":
            records.append(crest_record(smiles))
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


def chem_mechanism_spec_example() -> dict[str, Any]:
    return mechanism_spec_example()


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
