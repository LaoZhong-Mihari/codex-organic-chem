"""Human confirmation gate for structures, reactions, and OCSR results.

A review here is interactive: it opens the local review page so the user can
*fix* a structure (flip a wedge, pick another OCSR candidate, redraw it), not
just look at a picture of it. The static SVG preview stays in the payload for
transcript context, but the confirmation itself flows through the review
session; callers collect it with ``chem_structure_review_result``.

Whether the page opens is the user's setting (``review_mode``): ``auto`` opens
it, ``ask`` returns a one-line question to put to the user first, ``off``
returns the static payload only.
"""

from __future__ import annotations

from typing import Any

from .ocsr import parse_image
from .rdkit_tools import normalize_structure, reaction_to_svg
from .structure_review import open_review, review_gate


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
    *,
    confirmed: bool = False,
    detached: bool = False,
    **options: Any,
) -> dict[str, Any]:
    """Start an interactive review session and merge its handles into ``payload``.

    On any failure the static payload still stands on its own, with a warning
    explaining that confirmation has to happen by reply instead of in the page.
    """
    gate = review_gate(confirmed, what="this structure")
    if gate is not None:
        payload["review_mode"] = gate["review_mode"]
        payload["next_action"] = gate["next_action"]
        if gate["status"] == "ask_user_first":
            payload["user_prompt"] = gate["user_prompt"]
        else:
            payload["confirmation_required"] = False
        return payload
    try:
        session = open_review(items=items, wait=False, timeout_s=timeout_s, detached=detached, **options)
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
                "The review page has been opened for the user (see review_url). Wait for their "
                "confirmation or corrections with chem_structure_review_result(session_id=..., wait=true), "
                "then continue only with the returned canonical structures."
            ),
        }
    )
    payload["warnings"] = [*payload.get("warnings", []), *session.get("warnings", [])]
    return payload


def _reaction_component_items(reaction_smiles: str) -> list[dict[str, Any]] | None:
    """Split a reaction into per-component review items.

    The review queue shows one molecule per card, so each reactant, agent, and
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


def _other_candidates(candidates: list[Any], chosen: str) -> list[dict[str, str]]:
    """Distinct OCSR readings other than the top one, offered as alternatives."""
    seen = {chosen}
    alternatives: list[dict[str, str]] = []
    for candidate in candidates[1:]:
        if not isinstance(candidate, dict):
            continue
        smiles = candidate.get("isomeric_smiles") or candidate.get("canonical_smiles")
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)
        engine = candidate.get("engine") or candidate.get("source") or "OCSR"
        alternatives.append({"label": f"{engine} reading {len(alternatives) + 2}", "smiles": smiles})
        if len(alternatives) == 4:
            break
    return alternatives


def prepare_input_review(
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
    """``image_path`` runs OCSR on an image and reviews its top reading.
    ``smiles`` with ``source_image`` reviews *your* reading beside the image it
    came from, which is the usual case after transcribing a drawing."""
    session_options = {"confirmed": confirmed, "detached": detached, "language": language}
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
        candidates = parsed.get("candidates") or []
        top = candidates[0] if candidates else None
        if interactive and isinstance(top, dict):
            top_molecule = top.get("isomeric_smiles") or top.get("canonical_smiles")
            top_reaction = top.get("reaction_smiles")
            if top_molecule:
                item = {
                    "id": "ocsr-top",
                    "label": "Structure read from the image",
                    "smiles": top_molecule,
                    "image_path": image_path,
                    "candidates": _other_candidates(candidates, top_molecule),
                    "note": note,
                }
                result = _attach_review_session(result, [item], timeout_s=timeout_s, **session_options)
            elif top_reaction:
                items = _reaction_component_items(str(top_reaction))
                if items:
                    result = _attach_review_session(
                        result, items, timeout_s=timeout_s, images=[image_path], **session_options
                    )
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
                payload = _attach_review_session(payload, items, timeout_s=timeout_s, **session_options)
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
            item = {"id": "input", "label": "Input structure", "smiles": review_smiles, "note": note}
            if source_image:
                item["image_path"] = source_image
                item["label"] = "Structure read from the image"
            payload = _attach_review_session(payload, [item], timeout_s=timeout_s, **session_options)
        return payload
    return {
        "status": "error",
        "warnings": ["Provide image_path, smiles, reaction_smiles, or molfile for review."],
    }
