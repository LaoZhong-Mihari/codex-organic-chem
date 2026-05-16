from __future__ import annotations

import json
import os
import shlex
from pathlib import Path
from typing import Any

from .image_preprocessing import OcsrImageVariant, build_ocsr_image_variants
from .ocsr_adapters import OcsrAdapterSpec, adapter_specs_for_kind, default_adapter_command
from .rdkit_tools import normalize_structure, reaction_to_svg
from .scheme_ocsr import rank_candidates
from .utils import executable_status, run_command


def _candidate_from_payload(payload: dict[str, Any], default_tool: str) -> tuple[list[dict[str, Any]], list[str]]:
    warnings: list[str] = []
    candidates: list[dict[str, Any]] = []
    nested = payload.get("candidates")
    if isinstance(nested, list):
        for item in nested:
            if isinstance(item, dict):
                child_candidates, child_warnings = _candidate_from_payload(item, default_tool)
                candidates.extend(child_candidates)
                warnings.extend(child_warnings)
            elif item:
                candidates.append({"tool": default_tool, "smiles": str(item), "warnings": []})
    for key in ("warnings", "adapter_warnings"):
        payload_warnings = payload.get(key)
        if isinstance(payload_warnings, list):
            warnings.extend(str(item) for item in payload_warnings)
        elif payload_warnings:
            warnings.append(str(payload_warnings))
    value = payload.get("smiles") or payload.get("canonical_smiles") or payload.get("reaction_smiles")
    molfile = payload.get("molfile") or payload.get("molblock")
    if value or molfile:
        candidate = {
            "tool": str(payload.get("tool") or payload.get("source") or default_tool),
            "smiles": str(value) if value and not payload.get("reaction_smiles") else None,
            "reaction_smiles": str(payload.get("reaction_smiles")) if payload.get("reaction_smiles") else None,
            "molfile": str(molfile) if molfile else None,
            "confidence": payload.get("confidence"),
            "warnings": warnings.copy(),
        }
        for key in ("atoms", "bonds", "bboxes", "boxes", "labels"):
            if key in payload:
                candidate[key] = payload[key]
        candidates.append(candidate)
    return candidates, warnings


def _parse_candidate_lines(text: str, default_tool: str) -> tuple[list[dict[str, Any]], list[str]]:
    candidates: list[dict[str, Any]] = []
    warnings: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("{"):
            try:
                payload = json.loads(stripped)
                parsed_candidates, parse_warnings = _candidate_from_payload(payload, default_tool)
                warnings.extend(parse_warnings)
                if parsed_candidates:
                    candidates.extend(parsed_candidates)
                continue
            except json.JSONDecodeError:
                pass
        first = stripped.split()[0]
        if any(ch.isalpha() for ch in first) and not first.lower().startswith("error"):
            candidates.append({"tool": default_tool, "smiles": first, "warnings": []})
    return candidates, warnings


def _run_custom_adapter(
    spec: OcsrAdapterSpec,
    image_variant: OcsrImageVariant,
) -> tuple[list[dict[str, Any]], list[str]]:
    env_command = os.environ.get(spec.env_var)
    command = env_command or default_adapter_command(spec)
    if not command:
        return [], []
    input_path = str(image_variant.path)
    quoted_input = shlex.quote(input_path)
    rendered = command.format(input=quoted_input, input_raw=input_path, input_quoted=quoted_input)
    timeout_s = int(os.environ.get("CODEX_CHEM_OCR_TIMEOUT_S", "300"))
    try:
        code, stdout, stderr = run_command(shlex.split(rendered), timeout_s=timeout_s)
    except Exception as exc:
        return [], [f"{spec.env_var} adapter failed: {exc}"]
    warnings = []
    if code != 0:
        warnings.append(f"{spec.env_var} exited with code {code}: {stderr.strip()}")
    values, parse_warnings = _parse_candidate_lines(stdout, spec.name)
    warnings.extend(parse_warnings)
    for value in values:
        value["adapter_env"] = spec.env_var if env_command else f"default:{spec.name}"
        value["image_variant"] = image_variant.name
    return values, warnings


def _run_osra(image_variant: OcsrImageVariant) -> tuple[list[dict[str, Any]], list[str]]:
    status = executable_status("osra", ("--version",))
    if status.status != "available":
        return [], ["OSRA is unavailable; install osra or configure an OCSR adapter command."]
    try:
        code, stdout, stderr = run_command(["osra", "-f", "smi", str(image_variant.path)], timeout_s=120)
    except Exception as exc:
        return [], [f"OSRA execution failed: {exc}"]
    warnings = []
    if code != 0:
        warnings.append(f"OSRA exited with code {code}: {stderr.strip()}")
    values, parse_warnings = _parse_candidate_lines(stdout, "osra")
    warnings.extend(parse_warnings)
    for value in values:
        value["image_variant"] = image_variant.name
    return values, warnings


def _tools_payload(kind: str) -> dict[str, Any]:
    tools: dict[str, Any] = {}
    for spec in adapter_specs_for_kind(kind):
        command = os.environ.get(spec.env_var) or default_adapter_command(spec)
        tools[f"{spec.name}_command"] = {
            "name": f"{spec.name}_command",
            "status": "available" if command else "unavailable",
            "env_var": spec.env_var,
            "command": command,
            "purpose": spec.purpose,
        }
    tools["osra"] = executable_status("osra", ("--version",)).to_dict()
    return tools


def _dedupe_candidates(raw_candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[str, str, str]] = set()
    deduped: list[dict[str, Any]] = []
    for candidate in raw_candidates:
        tool = str(candidate.get("tool") or "unknown")
        value = str(candidate.get("reaction_smiles") or candidate.get("smiles") or candidate.get("molfile") or "")
        key = (tool, value, str(candidate.get("adapter_env") or ""))
        if value and key not in seen:
            deduped.append(candidate)
            seen.add(key)
    return deduped


def _coerce_confidence(value: Any) -> float | None:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_image(path: str, kind: str = "auto") -> dict:
    image_path = Path(path).expanduser().resolve()
    warnings: list[str] = []
    if not image_path.exists():
        return {
            "kind": kind,
            "path": str(image_path),
            "candidates": [],
            "tools": {},
            "warnings": [f"Image path does not exist: {image_path}"],
        }
    raw_candidates: list[dict[str, Any]] = []
    variant_set = build_ocsr_image_variants(image_path)
    warnings.extend(variant_set.warnings)
    try:
        for spec in adapter_specs_for_kind(kind):
            for image_variant in variant_set.variants:
                values, adapter_warnings = _run_custom_adapter(spec, image_variant)
                raw_candidates.extend(values)
                warnings.extend(adapter_warnings)
        if not raw_candidates:
            for image_variant in variant_set.variants:
                values, adapter_warnings = _run_osra(image_variant)
                raw_candidates.extend(values)
                warnings.extend(adapter_warnings)
    finally:
        variant_set.cleanup()
    candidates = []
    for candidate in _dedupe_candidates(raw_candidates):
        value = candidate.get("reaction_smiles") or candidate.get("smiles")
        molfile = candidate.get("molfile")
        adapter_confidence = _coerce_confidence(candidate.get("confidence"))
        candidate_warnings = [str(item) for item in candidate.get("warnings", [])]
        value_text = str(value or "")
        if value_text and (">>" in value_text or value_text.count(">") == 2):
            svg, draw_warnings = reaction_to_svg(value_text)
            candidates.append(
                {
                    "source": "ocsr",
                    "reaction_smiles": value_text,
                    "svg": svg,
                    "confidence": adapter_confidence if adapter_confidence is not None else 0.55,
                    "metadata": {
                        "ocsr_tool": candidate.get("tool"),
                        "adapter_env": candidate.get("adapter_env"),
                        "adapter_confidence": adapter_confidence,
                        "image_variant": candidate.get("image_variant"),
                    },
                    "warnings": [
                        "Reaction OCSR output is not fully validated; verify atom mapping, reagents, and stoichiometry manually.",
                        *candidate_warnings,
                        *draw_warnings,
                    ],
                }
            )
            continue
        if molfile:
            record = normalize_structure(molfile=molfile, source="ocsr").to_dict()
        else:
            record = normalize_structure(smiles=value_text, source="ocsr").to_dict()
        metadata = record.setdefault("metadata", {})
        metadata["ocsr_tool"] = candidate.get("tool")
        metadata["adapter_env"] = candidate.get("adapter_env")
        metadata["adapter_confidence"] = adapter_confidence
        metadata["image_variant"] = candidate.get("image_variant")
        metadata["raw_smiles"] = value_text
        for key in ("atoms", "bonds", "bboxes", "boxes", "labels"):
            if key in candidate:
                metadata[key] = candidate[key]
        if adapter_confidence is not None:
            record["confidence"] = min(float(record.get("confidence", adapter_confidence)), adapter_confidence)
        if candidate_warnings:
            record.setdefault("warnings", [])
            record["warnings"].extend(candidate_warnings)
        candidates.append(record)
    candidates = rank_candidates(candidates)
    if not candidates:
        warnings.append(
            "No structure candidates were produced. For complex scans or hand drawings, use Ketcher/ChemDraw-style manual correction."
        )
    return {
        "status": "awaiting_user_confirmation" if candidates else "no_candidates",
        "kind": kind,
        "path": str(image_path),
        "candidates": candidates,
        "confirmation_required": bool(candidates),
        "next_action": (
            "Review ranked OCSR candidates, resolve any dummy atoms or route-consistency warnings, then render the selected "
            "structure and wait for explicit confirmation before mechanism, synthesis, calculation, or literature workflows."
            if candidates
            else "Provide a clearer image or corrected SMILES/Molfile."
        ),
        "ranked_candidates": candidates,
        "tools": _tools_payload(kind),
        "warnings": warnings,
    }
