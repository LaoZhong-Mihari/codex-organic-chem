from __future__ import annotations

from typing import Any

from .input_review import prepare_input_review
from .literature import literature_search
from .rdkit_tools import descriptors_record, normalize_structure
from .reaction import retrosynthesis_hints


def _tool_summary(smiles: str, record: dict[str, Any]) -> dict[str, Any]:
    descriptor = descriptors_record(smiles).to_dict()
    return {
        "normalized_target": record,
        "descriptor_record": descriptor,
        "tool_facts": [
            "RDKit sanitized and canonicalized the target structure.",
            "Functional-group labels and descriptors are local RDKit/tool facts.",
            "No commercial retrosynthesis database was queried.",
        ],
    }


def _query_for_hint(hint: dict[str, str], target_record: dict[str, Any]) -> str:
    groups = " ".join(target_record.get("metadata", {}).get("functional_groups", []))
    strategy = hint.get("strategy", "")
    disconnection = hint.get("disconnection", "")
    formula = target_record.get("metadata", {}).get("formula", "")
    return f"organic synthesis {groups} {disconnection} {strategy} {formula}".strip()


def suggest_synthesis_path(
    target_smiles: str,
    stage: str = "first_disconnection",
    selected_option: str | None = None,
    confirmed: bool = False,
    literature: bool = True,
    literature_rows: int = 4,
) -> dict[str, Any]:
    target = normalize_structure(smiles=target_smiles, source="synthesis_target")
    target_dict = target.to_dict()
    if not confirmed:
        review = prepare_input_review(smiles=target_smiles)
        return {
            "status": "awaiting_user_confirmation",
            "workflow": "stepwise_synthesis_planning",
            "message": "Target structure must be rendered and confirmed before suggesting a route.",
            "review": review,
            "next_action": "After confirmation, rerun with confirmed=true/--confirmed. Do not expand a route yet.",
        }
    hints = retrosynthesis_hints([target.canonical_smiles or target_smiles])
    tool_data = _tool_summary(target.canonical_smiles or target_smiles, target_dict)
    options = []
    for idx, hint in enumerate(hints, start=1):
        option_id = f"disconnection-{idx}"
        lit_result = None
        if literature:
            lit_result = literature_search(_query_for_hint(hint, target_dict), rows=literature_rows)
        options.append(
            {
                "option_id": option_id,
                "stage": stage,
                "single_step_only": True,
                "proposal": hint,
                "tool_data_used": {
                    "functional_groups": target_dict.get("metadata", {}).get("functional_groups", []),
                    "formula": target_dict.get("metadata", {}).get("formula"),
                    "heavy_atom_count": target_dict.get("metadata", {}).get("heavy_atom_count"),
                    "descriptor_method": "rdkit_descriptors",
                },
                "literature_evidence": lit_result,
                "decision_needed_before_next_step": (
                    "Select this option, reject it, or provide constraints such as available starting materials, protected groups, "
                    "scale, forbidden reagents, or stereochemical requirements."
                ),
            }
        )
    return {
        "status": "awaiting_route_decision",
        "workflow": "stepwise_synthesis_planning",
        "stage": stage,
        "selected_option": selected_option,
        "target": target_dict,
        "tool_data": tool_data,
        "options": options,
        "next_action": (
            "Choose exactly one disconnection option or ask for a refined first-step comparison. "
            "The assistant must not output a full multistep route until the previous decision is confirmed."
        ),
        "evidence_policy": {
            "tool_facts": "RDKit descriptors, functional groups, structure warnings, and optional calculation outputs.",
            "literature_facts": "Crossref metadata returned with DOI/URL when available; user must inspect full papers for exact procedures/yields.",
            "llm_assumptions": "Any reagent/condition idea not directly present in tool or literature metadata is hypothesis only.",
        },
        "warnings": [
            "This is not a complete route. It is the first decision layer in a stepwise retrosynthesis workflow.",
            "Literature metadata is not equivalent to verified experimental scope, yield, or safety.",
            *target.warnings,
        ],
    }

