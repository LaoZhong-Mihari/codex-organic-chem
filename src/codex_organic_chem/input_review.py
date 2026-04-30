from __future__ import annotations

from typing import Any

from .ocsr import parse_image
from .rdkit_tools import normalize_structure, reaction_to_svg


def _confirmation_payload(kind: str, preview: dict[str, Any], warnings: list[str] | None = None) -> dict[str, Any]:
    return {
        "status": "awaiting_user_confirmation",
        "kind": kind,
        "preview": preview,
        "confirmation_required": True,
        "user_prompt": (
            "Please verify that the rendered structure/reaction matches the intended input. "
            "Reply with explicit confirmation, or provide corrected SMILES/Molfile/Rxnfile before continuing."
        ),
        "blocked_until_confirmed": [
            "mechanism_draft",
            "publication_figure_package",
            "synthesis_suggest",
            "literature-backed route expansion",
            "expensive quantum or conformer calculations",
        ],
        "warnings": warnings or [],
    }


def prepare_input_review(
    smiles: str | None = None,
    reaction_smiles: str | None = None,
    molfile: str | None = None,
    image_path: str | None = None,
    kind: str = "auto",
) -> dict[str, Any]:
    if image_path:
        parsed = parse_image(image_path, kind=kind)
        return {
            **parsed,
            "workflow_gate": {
                "must_render_before_continuing": True,
                "requires_explicit_user_confirmation": bool(parsed.get("candidates")),
                "acceptable_confirmation": "The user confirms candidate id/SMILES or supplies corrected machine-readable structure.",
            },
        }
    if reaction_smiles:
        svg, warnings = reaction_to_svg(reaction_smiles)
        return _confirmation_payload(
            "reaction",
            {"reaction_smiles": reaction_smiles, "svg": svg},
            warnings,
        )
    if smiles or molfile:
        record = normalize_structure(smiles=smiles, molfile=molfile, source="input_review")
        return _confirmation_payload(
            "molecule",
            record.to_dict(),
            record.warnings,
        )
    return {
        "status": "error",
        "warnings": ["Provide image_path, smiles, reaction_smiles, or molfile for review."],
    }

