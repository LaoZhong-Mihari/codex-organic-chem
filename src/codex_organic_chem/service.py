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
from .structure_review import (
    open_review,
    review_gate,
    review_mode_status,
    review_session_result,
    set_review_mode,
)
from . import image_tools as _image_tools
from . import mechanism_steps as _mechanism_steps
from . import model_kit as _model_kit
from . import route_check as _route_check
from . import structure_compare as _structure


_DRAWING_KEYS = {"svg", "molblock", "reaction_svg"}
_GEOMETRY_KEYS = {"atoms", "bonds"}


def _without_drawings(payload: Any, inside_metadata: bool = False) -> Any:
    """Drop SVG, Molfile and per-atom OCSR geometry from a tool result.

    They are for renderers. In an assistant's context they are tens of
    thousands of characters that say nothing the SMILES does not.
    """
    if isinstance(payload, dict):
        return {
            key: _without_drawings(value, inside_metadata or key == "metadata")
            for key, value in payload.items()
            if key not in _DRAWING_KEYS and not (inside_metadata and key in _GEOMETRY_KEYS)
        }
    if isinstance(payload, list):
        return [_without_drawings(value, inside_metadata) for value in payload]
    return payload


def _load_route(route: dict[str, Any] | None, route_path: str | None) -> dict[str, Any]:
    """A route given inline, or read from the JSON file it was written to."""
    if route is None and route_path:
        import json

        with open(Path(route_path).expanduser(), encoding="utf-8") as handle:
            route = json.load(handle)
    if not isinstance(route, dict):
        raise ValueError("Give the route as an object (route) or as the path of its JSON file (route_path).")
    return route


def chem_parse_image(
    path: str | None = None,
    kind: str = "auto",
    image_path: str | None = None,
    include_drawings: bool = False,
) -> dict[str, Any]:
    """OCSR candidates for an image. ``image_path`` is accepted as an alias of
    ``path``; drawings and per-atom geometry are left out unless asked for."""
    result = parse_image(path=path or image_path or "", kind=kind)
    return result if include_drawings else _without_drawings(result)


def chem_ocsr_benchmark(gold_smiles: str, image_dir: str) -> dict[str, Any]:
    return benchmark_ocsr(gold_smiles=gold_smiles, image_dir=image_dir)


def chem_parse_scheme(
    image: str,
    crops: str | list[dict[str, Any]],
    gold_map: str | None = None,
    include_drawings: bool = True,
) -> dict[str, Any]:
    result = parse_scheme(image=image, crops=crops, gold_map=gold_map)
    return result if include_drawings else _without_drawings(result)


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
    """Crop and magnify subregions of a page/scheme image, or render a coordinate grid."""
    return _image_tools.crop_image(
        image_path=image_path or path or "",
        boxes=boxes,
        output_dir=output_dir,
        grid=grid,
        scale=scale,
        padding_px=padding_px,
        normalized=normalized,
    )


def chem_input_review(
    smiles: str | None = None,
    reaction_smiles: str | None = None,
    molfile: str | None = None,
    image_path: str | None = None,
    kind: str = "auto",
    interactive: bool = True,
    timeout_s: int = 1800,
    confirmed: bool = False,
    detached: bool = False,
    language: str | None = None,
    note: str | None = None,
    source_image: str | None = None,
) -> dict[str, Any]:
    result = prepare_input_review(
        smiles=smiles,
        reaction_smiles=reaction_smiles,
        molfile=molfile,
        image_path=image_path,
        kind=kind,
        interactive=interactive,
        timeout_s=timeout_s,
        confirmed=confirmed,
        detached=detached,
        language=language,
        note=note,
        source_image=source_image,
    )
    if image_path:  # OCSR candidates: keep the readings, drop their drawings
        result = _without_drawings(result)
    return result


def chem_structure_review_batch(
    items: list[dict[str, Any]],
    wait: bool = False,
    timeout_s: int = 1800,
    title: str | None = None,
    intro: str | None = None,
    language: str | None = None,
    images: list[str] | None = None,
    confirmed: bool = False,
    detached: bool = False,
) -> dict[str, Any]:
    """Open the review page on a queue of structures.

    Each item needs ``smiles`` and may carry ``id``, ``label``, ``image_path``
    (the crop it was read from), ``note`` (how it was read, what was assumed),
    ``question`` (one thing to ask the reviewer), ``source`` and ``candidates``
    (alternative readings as SMILES or ``{label, smiles}``).
    """
    gate = review_gate(confirmed, what=f"{len(items) if isinstance(items, list) else 'these'} structure(s)")
    if gate is not None:
        return gate
    try:
        return open_review(
            items=items, wait=wait, timeout_s=timeout_s, detached=detached,
            title=title, intro=intro, language=language, images=images,
        )
    except RuntimeError as exc:  # the detached review server could not be reached
        return {"status": "error", "warnings": [str(exc)]}


def chem_route_review(
    route: dict[str, Any] | None = None,
    wait: bool = False,
    timeout_s: int = 3600,
    title: str | None = None,
    intro: str | None = None,
    language: str | None = None,
    confirmed: bool = False,
    detached: bool = False,
    route_path: str | None = None,
    images: list[str] | None = None,
) -> dict[str, Any]:
    """Open the review page on a route (same format as ``chem_route_check``).

    The reviewer sees the audit step by step, can correct any structure or
    the conditions of a step (the audit re-runs), answers audit prompts,
    comments, and approves the route or sends it back.
    """
    gate = review_gate(confirmed, what="this route step by step")
    if gate is not None:
        return gate
    try:
        return open_review(
            route=_load_route(route, route_path), wait=wait, timeout_s=timeout_s, detached=detached,
            title=title, intro=intro, language=language, images=images,
        )
    except (ValueError, RuntimeError, OSError) as exc:
        return {"status": "error", "warnings": [str(exc)]}


def chem_structure_review_result(
    session_id: str,
    review_token: str | None = None,
    wait: bool = False,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    """Result of a structure or route review, from this process or from disk."""
    return review_session_result(session_id=session_id, token=review_token, wait=wait, timeout_s=timeout_s)


def chem_review_mode(mode: str | None = None) -> dict[str, Any]:
    """Show, or with ``mode`` save, when the review page opens: auto, ask or off."""
    if mode is None:
        return review_mode_status()
    try:
        return set_review_mode(mode)
    except ValueError as exc:
        return {"status": "error", "warnings": [str(exc)]}


def chem_review_demo(language: str = "en", detached: bool = True, which: str = "both") -> dict[str, Any]:
    from .review_demo import review_demo

    return review_demo(language=language, detached=detached, which=which)


def chem_normalize_structure(
    smiles: str | None = None,
    molfile: str | None = None,
    include_drawings: bool = True,
) -> dict[str, Any]:
    result = normalize_structure(smiles=smiles, molfile=molfile, source="user").to_dict()
    return result if include_drawings else _without_drawings(result)


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


def chem_literature_search(
    query: str, rows: int = 5, source: str = "auto", year_from: int | None = None
) -> dict[str, Any]:
    return literature_search(query=query, rows=rows, source=source, year_from=year_from)


def chem_reaction_analyze(
    input: str,
    mode: str = "sanity_check",
    include_drawings: bool = True,
) -> dict[str, Any]:
    canonical_mode = "sanity_check" if mode == "mass_charge_balance" else mode
    if canonical_mode not in {"forward", "retro", "conditions", "sanity_check"}:
        return {
            "status": "error",
            "warnings": [
                f"Unsupported reaction analysis mode: {mode}. "
                "Use 'sanity_check' (alias 'mass_charge_balance'), 'forward', 'retro', or 'conditions'."
            ],
        }
    result = analyze_reaction(input, mode=canonical_mode).to_dict()
    return result if include_drawings else _without_drawings(result)


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


# --------------------------------------------------------------------------
# verification kit: checks whose answer differs between a right and a wrong
# structure, route or 3D claim


def _guard(function: Any, *args: Any, **kwargs: Any) -> dict[str, Any]:
    try:
        return function(*args, **kwargs)
    except _structure.StructureError as exc:
        return {"status": "error", "warnings": [str(exc)]}


def chem_compare_structures(
    smiles_a: str | None = None,
    smiles_b: str | None = None,
    a: str | None = None,
    b: str | None = None,
) -> dict[str, Any]:
    """Relationship between two structures, or a stereo audit of one.

    ``a`` / ``b`` are accepted as aliases, matching ``kit.compare(a, b)`` and
    the CLI flags.
    """
    smiles_a, smiles_b = smiles_a or a or "", smiles_b or b
    if not smiles_b:
        return _guard(_structure.stereo_audit, smiles_a)
    return _guard(_structure.compare_structures, smiles_a, smiles_b)


ENUMERATE_TASKS = ("stereoisomers", "substitution", "environments", "groups")


def chem_enumerate(
    smiles: str,
    task: str = "stereoisomers",
    substituent: str = "Cl",
    n_substitutions: int = 1,
) -> dict[str, Any]:
    """Count instead of guessing: stereoisomers, substitution products, NMR environments, groups."""
    if task == "stereoisomers":
        return _guard(_structure.enumerate_stereoisomers, smiles)
    if task == "substitution":
        return _guard(_structure.substitution_products, smiles, substituent, n_substitutions=n_substitutions)
    if task == "environments":
        return _guard(_structure.symmetry_environments, smiles)
    if task == "groups":
        return _guard(_route_check.functional_groups, smiles)
    return {"status": "error", "warnings": [f"Unknown task {task!r}; use one of {ENUMERATE_TASKS}."]}


def chem_model_kit(
    task: str,
    smiles: str | list[str],
    atoms: list[int] | None = None,
    by_map: bool = False,
    method: str = "mmff",
    solvent: str | None = None,
    threshold: float = 4.0,
    refine: str = "auto",
) -> dict[str, Any]:
    """3D model-kit questions: face accessibility, reach, dihedrals, planarity, feasibility, isomer energies."""
    options: dict[str, Any] = {}
    if task in {"distance", "dihedral", "face", "cis_trans"}:
        options.update({"atoms": atoms or [], "by_map": by_map})
    if task == "distance":
        options["threshold"] = threshold
    if task == "planarity":
        options.update({"atoms": atoms, "by_map": by_map, "refine": refine, "solvent": solvent})
    if task == "isomer_energies":
        options.update({"method": method, "solvent": solvent})
    return _model_kit.model_kit(task, smiles, **options)


def chem_mechanism_steps_check(
    steps: list[str | dict[str, Any]] | None = None,
    reactants: str | list[str] | None = None,
    products: str | list[str] | None = None,
    conditions: str | None = None,
    reaction: str | None = None,
    steps_path: str | None = None,
    issues_only: bool = False,
) -> dict[str, Any]:
    """Audit a proposed sequence of reaction intermediates step by step."""
    if steps is None and steps_path:
        import json

        try:
            loaded = json.loads(Path(steps_path).expanduser().read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            return {"status": "error", "warnings": [f"Could not read steps_path: {exc}"]}
        if isinstance(loaded, dict):
            steps = loaded.get("steps") or loaded.get("mechanism") or []
            reactants = reactants if reactants is not None else (loaded.get("reactants") or loaded.get("reactants_smiles"))
            products = products if products is not None else (loaded.get("products") or loaded.get("products_smiles"))
            conditions = conditions if conditions is not None else loaded.get("conditions")
            reaction = reaction if reaction is not None else loaded.get("reaction")
        elif isinstance(loaded, list):
            steps = loaded
    report = _mechanism_steps.mechanism_steps_check(
        steps=steps or [],
        reactants=reactants,
        products=products,
        conditions=conditions,
        reaction=reaction,
    )
    if issues_only and isinstance(report, dict):
        report = {k: v for k, v in report.items() if k != "steps"}
    return report


def chem_route_check(
    route: dict[str, Any] | None = None,
    figure_output_dir: str | None = None,
    mcs_timeout_s: int = 6,
    route_path: str | None = None,
    figure_dir: str | None = None,
    issues_only: bool = False,
) -> dict[str, Any]:
    """Audit a route step by step; optionally draw every numbered compound as a scheme."""
    try:
        route = _load_route(route, route_path)
    except (ValueError, OSError) as exc:
        return {"status": "error", "warnings": [str(exc)]}
    report = _guard(_route_check.route_check, route, mcs_timeout_s=mcs_timeout_s)
    out_dir = figure_output_dir or figure_dir
    if out_dir and report.get("status") not in {"error", "unavailable"}:
        try:
            spec = _route_check.route_to_figure_spec(route)
            figure = render_route_figure(spec=spec, output_dir=out_dir, basename="route", formats=["svg", "png"])
            report["figure"] = {
                "status": figure.get("status"),
                "files": figure.get("files") or figure.get("output_files"),
                "warnings": figure.get("warnings", []),
                "note": "Open the image and compare every drawn structure with what you intended.",
            }
        except Exception as exc:
            report["figure"] = {"status": "error", "warnings": [f"Route figure failed: {exc}"]}
    if issues_only and isinstance(report, dict):
        report = {k: v for k, v in report.items() if k != "steps"}
    return report


def chem_route_check_example() -> dict[str, Any]:
    return _route_check.route_check_example()


PYTHON_ENV_FLAG = "CODEX_CHEM_ENABLE_PYTHON"


def python_tool_enabled() -> bool:
    import os

    return os.environ.get(PYTHON_ENV_FLAG, "").strip().lower() in {"1", "true", "yes"}


def _scrubbed_python_env() -> dict[str, str]:
    import os
    import re

    secret_pattern = re.compile(r"(API_KEY|SECRET|TOKEN|PASSWORD|CREDENTIAL|AUTH|COOKIE|SESSION|PRIVATE_KEY)", re.IGNORECASE)
    blocked_exact = {"SSH_AUTH_SOCK", "SSH_AGENT_PID", "AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"}
    clean: dict[str, str] = {}
    for key, value in os.environ.items():
        if key in blocked_exact or secret_pattern.search(key):
            continue
        clean[key] = value
    return clean


def _clip_output(text: str, limit: int) -> tuple[str, bool]:
    if len(text) <= limit:
        return text, False
    half = max(1, (limit - 40) // 2)
    omitted = len(text) - 2 * half
    return f"{text[:half]}\n...[truncated {omitted} chars]...\n{text[-half:]}", True


def chem_python(code: str, timeout_s: int = 120) -> dict[str, Any]:
    """Run a Python snippet with RDKit and this toolkit importable (opt-in).

    Meant for hosts without a shell. It executes arbitrary code as the current
    user, so it stays disabled unless CODEX_CHEM_ENABLE_PYTHON=1 is set for the
    MCP server process.
    """
    import os
    import signal
    import subprocess
    import sys
    import tempfile

    if not python_tool_enabled():
        return {
            "status": "disabled",
            "warnings": [
                f"Code execution is off. Set {PYTHON_ENV_FLAG}=1 for the MCP server to enable it, or run "
                "`uv run python` from a shell: `from codex_organic_chem import kit`."
            ],
        }
    limit = 20000
    bounded_timeout = max(1, min(int(timeout_s), 900))
    with tempfile.TemporaryDirectory(prefix="codex-chem-python-") as work:
        proc = subprocess.Popen(
            [sys.executable, "-I", "-c", code],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=work,
            env=_scrubbed_python_env(),
            start_new_session=True,
        )
        try:
            out_raw, err_raw = proc.communicate(timeout=bounded_timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except OSError:
                proc.kill()
            proc.communicate()
            return {"status": "timeout", "warnings": [f"Snippet exceeded {bounded_timeout} s."]}
    stdout, out_trunc = _clip_output(out_raw or "", limit)
    stderr, err_trunc = _clip_output(err_raw or "", limit)
    return {
        "status": "ok" if proc.returncode == 0 else "error",
        "exit_code": proc.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "truncated": out_trunc or err_trunc,
    }
