from __future__ import annotations

import json
import os
import shlex
from pathlib import Path

from .rdkit_tools import normalize_structure, reaction_to_svg
from .utils import executable_status, run_command


def _parse_smiles_lines(text: str) -> list[str]:
    smiles: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("{"):
            try:
                payload = json.loads(stripped)
                value = payload.get("smiles") or payload.get("canonical_smiles") or payload.get("reaction_smiles")
                if value:
                    smiles.append(str(value))
                    continue
            except json.JSONDecodeError:
                pass
        first = stripped.split()[0]
        if any(ch.isalpha() for ch in first) and not first.lower().startswith("error"):
            smiles.append(first)
    return smiles


def _run_custom_adapter(env_name: str, image_path: Path) -> tuple[list[str], list[str]]:
    command = os.environ.get(env_name)
    if not command:
        return [], [f"{env_name} is not set."]
    input_path = str(image_path)
    quoted_input = shlex.quote(input_path)
    rendered = command.format(input=quoted_input, input_raw=input_path, input_quoted=quoted_input)
    timeout_s = int(os.environ.get("CODEX_CHEM_OCR_TIMEOUT_S", "300"))
    try:
        code, stdout, stderr = run_command(shlex.split(rendered), timeout_s=timeout_s)
    except Exception as exc:
        return [], [f"{env_name} adapter failed: {exc}"]
    warnings = []
    if code != 0:
        warnings.append(f"{env_name} exited with code {code}: {stderr.strip()}")
    return _parse_smiles_lines(stdout), warnings


def _run_osra(image_path: Path) -> tuple[list[str], list[str]]:
    status = executable_status("osra", ("--version",))
    if status.status != "available":
        return [], ["OSRA is unavailable; install osra or configure MolScribe/RxnScribe adapter."]
    try:
        code, stdout, stderr = run_command(["osra", "-f", "smi", str(image_path)], timeout_s=120)
    except Exception as exc:
        return [], [f"OSRA execution failed: {exc}"]
    warnings = []
    if code != 0:
        warnings.append(f"OSRA exited with code {code}: {stderr.strip()}")
    return _parse_smiles_lines(stdout), warnings


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
    raw_candidates: list[str] = []
    if kind in {"auto", "molecule"}:
        values, adapter_warnings = _run_custom_adapter("CODEX_CHEM_MOLSCRIBE_CMD", image_path)
        raw_candidates.extend(values)
        warnings.extend(adapter_warnings)
    if kind in {"auto", "reaction"}:
        values, adapter_warnings = _run_custom_adapter("CODEX_CHEM_RXNSCRIBE_CMD", image_path)
        raw_candidates.extend(values)
        warnings.extend(adapter_warnings)
    if not raw_candidates:
        values, adapter_warnings = _run_osra(image_path)
        raw_candidates.extend(values)
        warnings.extend(adapter_warnings)
    deduped: list[str] = []
    for value in raw_candidates:
        if value not in deduped:
            deduped.append(value)
    candidates = []
    for value in deduped:
        if ">>" in value or value.count(">") == 2:
            svg, draw_warnings = reaction_to_svg(value)
            candidates.append(
                {
                    "source": "ocsr",
                    "reaction_smiles": value,
                    "svg": svg,
                    "confidence": 0.55,
                    "warnings": [
                        "Reaction OCSR output is not fully validated; verify atom mapping, reagents, and stoichiometry manually.",
                        *draw_warnings,
                    ],
                }
            )
            continue
        candidates.append(normalize_structure(smiles=value, source="ocsr").to_dict())
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
            "Render the candidate SVG/MolBlock to the user and wait for explicit confirmation or corrected SMILES/Molfile before "
            "running mechanism, synthesis, calculation, or literature workflows."
            if candidates
            else "Provide a clearer image or corrected SMILES/Molfile."
        ),
        "tools": {
            "molscribe_command": "configured" if os.environ.get("CODEX_CHEM_MOLSCRIBE_CMD") else "unavailable",
            "rxnscribe_command": "configured" if os.environ.get("CODEX_CHEM_RXNSCRIBE_CMD") else "unavailable",
            "osra": executable_status("osra", ("--version",)).to_dict(),
        },
        "warnings": warnings,
    }
