"""Human confirmation gate for structures, reactions, and OCSR results.

Every review here is interactive by default: it opens the local editor UI
(Ketcher when built, the built-in fallback otherwise) so the user can *fix* a
structure, not just look at a picture of it. The static SVG preview is kept in
the payload for transcript context, but the confirmation itself flows through
the editor session — callers should wait on ``chem_structure_review_result``.
"""

from __future__ import annotations

from typing import Any

from .ocsr import parse_image
from .rdkit_tools import normalize_structure, reaction_to_svg
from .structure_review import start_structure_review_batch


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


def _attach_review_session(
    payload: dict[str, Any],
    items: list[dict[str, Any]],
    timeout_s: int = 1800,
) -> dict[str, Any]:
    """Start an interactive editor session and merge its handles into ``payload``.

    On any failure the static payload still stands on its own, with a warning
    explaining that confirmation has to happen by reply instead of in-editor.
    """
    try:
        session = start_structure_review_batch(items=items, wait=False, timeout_s=timeout_s)
    except Exception as exc:
        payload["warnings"] = [*payload.get("warnings", []), f"Interactive review session unavailable: {exc}"]
        return payload
    payload.update(
        {
            "session_id": session["session_id"],
            "review_token": session["review_token"],
            "review_url": session["review_url"],
            "editor": session.get("editor"),
            "browser_opened": session.get("browser_opened", False),
            "items": session["items"],
            "next_action": (
                "The structure editor has been opened for the user (see review_url). Wait for their "
                "confirmation or corrections with chem_structure_review_result(session_id=..., wait=true), "
                "then continue only with the returned canonical structures."
            ),
        }
    )
    payload["warnings"] = [*payload.get("warnings", []), *session.get("warnings", [])]
    return payload


def _reaction_component_items(reaction_smiles: str) -> list[dict[str, Any]] | None:
    """Split a reaction into per-component editor items.

    The review queue edits one molecule per card, so each reactant, agent, and
    product becomes its own item; the caller reassembles the corrected reaction
    from the confirmed parts. Returns None when the string is not a 3-part
    reaction SMILES.
    """
    parts = reaction_smiles.split(">")
    if len(parts) != 3:
        return None
    roles = ("reactant", "agent", "product")
    items: list[dict[str, Any]] = []
    for role, part in zip(roles, parts):
        components = [component for component in part.split(".") if component.strip()]
        for index, component in enumerate(components, start=1):
            suffix = f" {index}" if len(components) > 1 else ""
            items.append(
                {
                    "id": f"{role}{'-' + str(index) if len(components) > 1 else ''}",
                    "label": f"{role.capitalize()}{suffix}: {component}",
                    "smiles": component,
                }
            )
    return items or None


def prepare_input_review(
    smiles: str | None = None,
    reaction_smiles: str | None = None,
    molfile: str | None = None,
    image_path: str | None = None,
    kind: str = "auto",
    interactive: bool = True,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    if image_path:
        parsed = parse_image(image_path, kind=kind)
        result = {
            **parsed,
            "workflow_gate": {
                "must_render_before_continuing": True,
                "requires_explicit_user_confirmation": bool(parsed.get("candidates")),
                "acceptable_confirmation": "The user confirms candidate id/SMILES or supplies corrected machine-readable structure.",
            },
        }
        top = (parsed.get("candidates") or [None])[0]
        if interactive and isinstance(top, dict):
            top_molecule = top.get("isomeric_smiles") or top.get("canonical_smiles")
            top_reaction = top.get("reaction_smiles")
            if top_molecule:
                result = _attach_review_session(
                    result,
                    [{"id": "ocsr-top", "label": "Top OCSR candidate (edit to correct)", "smiles": top_molecule}],
                    timeout_s=timeout_s,
                )
            elif top_reaction:
                items = _reaction_component_items(str(top_reaction))
                if items:
                    result = _attach_review_session(result, items, timeout_s=timeout_s)
        return result
    if reaction_smiles:
        svg, warnings = reaction_to_svg(reaction_smiles)
        payload = _confirmation_payload(
            "reaction",
            {"reaction_smiles": reaction_smiles, "svg": svg},
            warnings,
        )
        if interactive:
            items = _reaction_component_items(reaction_smiles)
            if items:
                payload = _attach_review_session(payload, items, timeout_s=timeout_s)
        return payload
    if smiles or molfile:
        record = normalize_structure(smiles=smiles, molfile=molfile, source="input_review")
        payload = _confirmation_payload(
            "molecule",
            record.to_dict(),
            record.warnings,
        )
        review_smiles = record.isomeric_smiles or record.canonical_smiles or smiles
        if interactive and review_smiles:
            payload = _attach_review_session(
                payload,
                [{"id": "input", "label": "Input structure (edit to correct)", "smiles": review_smiles}],
                timeout_s=timeout_s,
            )
        return payload
    return {
        "status": "error",
        "warnings": ["Provide image_path, smiles, reaction_smiles, or molfile for review."],
    }
