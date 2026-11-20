"""Route review: the route audit, laid out for a person to check and correct.

The state keeps the route as the assistant proposed it next to the route as
the reviewer left it. Every edit (a structure, the conditions of a step) re-runs
``route_check`` so the reviewer sees at once what their correction does to atom
balance, configuration tracking and target identity. What goes back to the
assistant is the corrected route plus a decision and comment per step.
"""

from __future__ import annotations

import copy
from typing import Any

from . import review_chem as rc
from . import structure_compare as sc
from .route_check import route_check

TARGET_ID = "__target__"
STEP_DECISIONS = ("approved", "flagged")
ISSUE_ANSWERS = ("fine", "problem")
VERDICTS = ("approved", "needs_changes")

_INVERTED = (1.0, 0.72, 0.30)
_CREATED = (0.62, 0.88, 0.62)
_PROBLEM = (1.0, 0.62, 0.62)

# stereo event lists of a step audit -> (kind shown on the page, highlight)
_STEREO_EVENTS = (
    ("inverted_with_substituent_replaced", "inverted", _INVERTED),
    ("inverted_where_only_heteroatom_derivatised", "inverted_unexpected", _PROBLEM),
    ("epimerized_without_bond_change", "epimerized", _PROBLEM),
    ("created_defined", "created", _CREATED),
    ("created_undefined", "created_undefined", _PROBLEM),
    ("lost_definition", "lost_definition", _PROBLEM),
    ("destroyed", "destroyed", None),
    ("ambiguous", "ambiguous", None),
)


def _smiles_of(entry: Any) -> str:
    if isinstance(entry, dict):
        return str(entry.get("smiles") or "")
    return str(entry or "")


def _name_of(entry: Any) -> str:
    return str(entry.get("name") or "") if isinstance(entry, dict) else ""


def _step_key(step: dict[str, Any], position: int) -> str:
    return str(step.get("step", position + 1))


class RouteReview:
    def __init__(self, route: dict[str, Any], mcs_timeout_s: int = 4) -> None:
        if not isinstance(route, dict) or not isinstance(route.get("compounds"), dict) or not route.get("steps"):
            raise ValueError("A route needs 'compounds' (id -> SMILES) and a non-empty 'steps' list.")
        self.original = copy.deepcopy(route)
        self.current = copy.deepcopy(route)
        self.mcs_timeout_s = mcs_timeout_s
        self.step_decisions: dict[str, dict[str, str]] = {}
        self.issue_answers: dict[str, dict[str, str]] = {}
        self.verdict: str | None = None
        self.comment = ""
        self._audit: dict[str, Any] | None = None

    # ------------------------------------------------------------- edits

    def _invalidate(self) -> None:
        self._audit = None

    def set_compound(self, compound_id: str, smiles: str | None = None, molfile: str | None = None) -> None:
        described = rc.describe(smiles, molfile)
        if described["status"] != "ok":
            raise ValueError("; ".join(described["warnings"]) or "The structure could not be read.")
        new_smiles = described["smiles"]
        if compound_id == TARGET_ID:
            self.current["target"] = new_smiles
        else:
            compounds = self.current["compounds"]
            if compound_id not in compounds:
                raise KeyError(f"Unknown compound: {compound_id}")
            if isinstance(compounds[compound_id], dict):
                compounds[compound_id]["smiles"] = new_smiles
            else:
                compounds[compound_id] = new_smiles
        self._invalidate()

    def reset_compound(self, compound_id: str | None = None) -> None:
        if compound_id is None:
            self.current = copy.deepcopy(self.original)
        elif compound_id == TARGET_ID:
            self.current["target"] = self.original.get("target")
        elif compound_id in self.original["compounds"]:
            self.current["compounds"][compound_id] = copy.deepcopy(self.original["compounds"][compound_id])
        else:
            raise KeyError(f"Unknown compound: {compound_id}")
        self._invalidate()

    def _step(self, key: str) -> dict[str, Any]:
        for position, step in enumerate(self.current["steps"]):
            if _step_key(step, position) == str(key):
                return step
        raise KeyError(f"Unknown step: {key}")

    def update_step(
        self,
        key: str,
        conditions: str | None = None,
        decision: str | None = None,
        comment: str | None = None,
    ) -> None:
        step = self._step(key)
        if conditions is not None and conditions != step.get("conditions"):
            step["conditions"] = str(conditions)
            self._invalidate()
        entry = self.step_decisions.setdefault(str(key), {"decision": "", "comment": ""})
        if decision is not None:
            if decision and decision not in STEP_DECISIONS:
                raise ValueError(f"Step decision must be one of {', '.join(STEP_DECISIONS)} or empty.")
            entry["decision"] = decision
        if comment is not None:
            entry["comment"] = str(comment)

    def answer_issue(self, key: str, answer: str | None, comment: str | None = None) -> None:
        if answer and answer not in ISSUE_ANSWERS:
            raise ValueError(f"Issue answer must be one of {', '.join(ISSUE_ANSWERS)} or empty.")
        entry = self.issue_answers.setdefault(str(key), {"answer": "", "comment": ""})
        entry["answer"] = answer or ""
        if comment is not None:
            entry["comment"] = str(comment)

    def finish(self, verdict: str | None, comment: str | None = None) -> None:
        if verdict and verdict not in VERDICTS:
            raise ValueError(f"Verdict must be one of {', '.join(VERDICTS)}.")
        self.verdict = verdict or self.verdict
        if comment is not None:
            self.comment = str(comment)

    # ------------------------------------------------------------- audit

    def audit(self) -> dict[str, Any]:
        if self._audit is None:
            try:
                self._audit = route_check(self.current, mcs_timeout_s=self.mcs_timeout_s)
            except Exception as exc:  # a reviewer's edit must never take the page down
                self._audit = {
                    "status": "error",
                    "issues": [{"severity": "blocking", "code": "audit_failed", "message": str(exc)}],
                    "steps": [],
                }
        return self._audit

    def _changes(self) -> tuple[dict[str, Any], dict[str, Any]]:
        compounds: dict[str, Any] = {}
        pairs = [(cid, _smiles_of(entry)) for cid, entry in self.original["compounds"].items()]
        pairs.append((TARGET_ID, str(self.original.get("target") or "")))
        for cid, before in pairs:
            after = (
                str(self.current.get("target") or "")
                if cid == TARGET_ID
                else _smiles_of(self.current["compounds"].get(cid))
            )
            if not before or not after or before == after:
                continue
            try:
                relation = sc.compare_structures(before, after)
            except sc.StructureError:
                relation = {"relation": "unparseable"}
            if relation["relation"] == "identical":
                continue
            entry = {"from": before, "to": after, "relation": relation["relation"]}
            for field in ("formula_delta_b_minus_a", "n_differing", "differing_elements", "note"):
                if field in relation:
                    entry[field] = relation[field]
            compounds[cid] = entry
        conditions: dict[str, Any] = {}
        for position, step in enumerate(self.original["steps"]):
            key = _step_key(step, position)
            try:
                now = self._step(key).get("conditions")
            except KeyError:
                continue
            if (step.get("conditions") or "") != (now or ""):
                conditions[key] = {"from": step.get("conditions") or "", "to": now or ""}
        return compounds, conditions

    # -------------------------------------------------------------- view

    def _compound_view(self, cid: str, entry: Any, role: str, changed: dict[str, Any]) -> dict[str, Any]:
        smiles = _smiles_of(entry)
        described = rc.describe(smiles, width=250, height=180)
        view = {
            "id": cid,
            "name": _name_of(entry),
            "role": role,
            "smiles": smiles,
            "original_smiles": changed.get(cid, {}).get("from", smiles),
            "changed": cid in changed,
            "relation_to_original": changed.get(cid, {}).get("relation"),
            "status": described["status"],
            "warnings": described.get("warnings", []),
        }
        for field in ("svg", "formula", "mol_weight", "n_stereo_elements", "n_undefined"):
            view[field] = described.get(field)
        return view

    def _product_svg(self, smiles: str, stereo: dict[str, Any]) -> str | None:
        highlight: dict[int, tuple[float, float, float]] = {}
        for field, _kind, colour in _STEREO_EVENTS:
            if colour is None:
                continue
            for event in stereo.get(field) or []:
                idx = event.get("product_atom") if isinstance(event, dict) else None
                if isinstance(idx, int):
                    highlight[idx] = colour
        if not highlight:
            return None
        try:
            mol = sc.parse(smiles)
            highlight = {idx: colour for idx, colour in highlight.items() if idx < mol.GetNumAtoms()}
            return rc.depict(mol, highlight=highlight, width=250, height=180)
        except Exception:
            return None

    @staticmethod
    def _stereo_chips(stereo: dict[str, Any]) -> list[dict[str, Any]]:
        chips: list[dict[str, Any]] = []
        if stereo.get("n_retained"):
            chips.append({"kind": "retained", "n": stereo["n_retained"]})
        for field, kind, _colour in _STEREO_EVENTS:
            events = stereo.get(field) or []
            if not events:
                continue
            details = []
            for event in events:
                if not isinstance(event, dict):
                    continue
                atom = event.get("product_atom", event.get("reactant_atom"))
                name = f"{event.get('element', 'C')}{atom + 1}" if isinstance(atom, int) else event.get("element", "")
                before, after = event.get("cip_before"), event.get("cip_after") or event.get("cip")
                if before and after:
                    details.append(f"{name} {before}→{after}")
                elif after or before:
                    details.append(f"{name} ({after or before})")
                else:
                    details.append(str(name))
            chips.append({"kind": kind, "n": len(events), "details": details})
        return chips

    def view(self) -> dict[str, Any]:
        audit = self.audit()
        changed_compounds, changed_conditions = self._changes()
        compounds = self.current["compounds"]
        produced = {str(step.get("product")) for step in self.current["steps"]}
        final = audit.get("final_product")
        compound_views: dict[str, Any] = {}
        for cid, entry in compounds.items():
            role = "final" if cid == final else ("intermediate" if cid in produced else "starting_material")
            compound_views[cid] = self._compound_view(cid, entry, role, changed_compounds)
        audit_steps = {str(s.get("step")): s for s in audit.get("steps", [])}
        steps = []
        for position, step in enumerate(self.current["steps"]):
            key = _step_key(step, position)
            audited = audit_steps.get(key, {})
            issues = []
            seen: dict[str, int] = {}
            for issue in audited.get("issues", []):
                seen[issue["code"]] = seen.get(issue["code"], 0) + 1
                issue_key = f"{key}:{issue['code']}" + (f":{seen[issue['code']]}" if seen[issue["code"]] > 1 else "")
                issues.append(
                    {
                        "key": issue_key,
                        "severity": issue["severity"],
                        "code": issue["code"],
                        "message": issue["message"],
                        **self.issue_answers.get(issue_key, {"answer": "", "comment": ""}),
                    }
                )
            product_id = str(step.get("product"))
            product_smiles = _smiles_of(compounds.get(product_id)) if product_id in compounds else product_id
            balance = audited.get("atom_balance") or {}
            steps.append(
                {
                    "key": key,
                    "reactants": [str(r) for r in step.get("reactants", [])],
                    "product": product_id,
                    "conditions": step.get("conditions") or "",
                    "conditions_changed": key in changed_conditions,
                    "original_conditions": changed_conditions.get(key, {}).get("from"),
                    "why": step.get("why") or "",
                    "reagents": step.get("reagents") or [],
                    "byproducts": step.get("byproducts") or [],
                    "balance": {
                        "delta": balance.get("delta_product_minus_reactants"),
                        "status": balance.get("status"),
                        "matches": balance.get("matches") or [],
                    },
                    "stereo": self._stereo_chips(audited.get("stereo") or {}),
                    "product_svg": self._product_svg(product_smiles, audited.get("stereo") or {}),
                    "bystanders": (audited.get("functional_groups") or {}).get("bystanders_that_must_survive") or {},
                    "issues": issues,
                    **self.step_decisions.get(key, {"decision": "", "comment": ""}),
                }
            )
        target_smiles = str(self.current.get("target") or "")
        target = None
        if target_smiles:
            target = self._compound_view(TARGET_ID, {"smiles": target_smiles, "name": "target"}, "target", changed_compounds)
        route_issues = [
            {"severity": i["severity"], "code": i["code"], "message": i["message"]}
            for i in audit.get("issues", [])
            if "step" not in i
        ]
        return {
            "title": self.current.get("title") or "",
            "status": audit.get("status"),
            "n_blocking": audit.get("n_blocking", 0),
            "n_review": audit.get("n_review", 0),
            "n_steps": len(steps),
            "longest_linear_sequence": audit.get("longest_linear_sequence"),
            "final_product": final,
            "final_vs_target": (audit.get("final_vs_target") or {}).get("relation"),
            "target": target,
            "compounds": compound_views,
            "steps": steps,
            "route_issues": route_issues,
            "scope": audit.get("scope"),
            "n_edited_compounds": len(changed_compounds),
            "n_edited_conditions": len(changed_conditions),
            "verdict": self.verdict,
            "comment": self.comment,
        }

    # ------------------------------------------------------------ result

    def result(self) -> dict[str, Any]:
        """What the assistant gets back: the corrected route and every decision."""
        audit = self.audit()
        changed_compounds, changed_conditions = self._changes()
        decisions = [
            {"step": key, **value}
            for key, value in self.step_decisions.items()
            if value.get("decision") or value.get("comment")
        ]
        answers = [
            {"issue": key, **value}
            for key, value in self.issue_answers.items()
            if value.get("answer") or value.get("comment")
        ]
        flagged = [d["step"] for d in decisions if d.get("decision") == "flagged"]
        problems = [a["issue"] for a in answers if a.get("answer") == "problem"]
        return {
            "verdict": self.verdict,
            "comment": self.comment,
            "route": copy.deepcopy(self.current),
            "changed_compounds": changed_compounds,
            "changed_conditions": changed_conditions,
            "step_decisions": decisions,
            "issue_answers": answers,
            "needs_revision": bool(
                self.verdict == "needs_changes" or flagged or problems or audit.get("n_blocking")
            ),
            "audit_after_review": {
                "status": audit.get("status"),
                "n_blocking": audit.get("n_blocking", 0),
                "n_review": audit.get("n_review", 0),
                "final_vs_target": (audit.get("final_vs_target") or {}).get("relation"),
                "issues": audit.get("issues", []),
            },
        }


def strip_drawings(payload: Any) -> Any:
    """Copy of a payload without SVG strings (for transcripts and disk)."""
    if isinstance(payload, dict):
        return {k: strip_drawings(v) for k, v in payload.items() if k not in {"svg", "product_svg"}}
    if isinstance(payload, list):
        return [strip_drawings(v) for v in payload]
    return payload
