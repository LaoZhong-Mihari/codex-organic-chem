"""Check electron bookkeeping independently of drawing coordinates.

This is a bounded main-group graph checker, not a mechanism predictor or an
energetic/kinetic calculation. Arrows in one step act on the SAME starting
state; their edits are committed together, never as fictitious intermediates.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from .rdkit_tools import RDKIT_AVAILABLE

if RDKIT_AVAILABLE:
    from rdkit import Chem


def parse_molecule(spec: dict[str, Any]) -> Any:
    """Preserve explicit mapped H: it may be an electron-arrow endpoint."""
    if not RDKIT_AVAILABLE:
        return None
    if spec.get("smiles"):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        return Chem.MolFromSmiles(spec["smiles"], params)
    if spec.get("molfile"):
        return Chem.MolFromMolBlock(spec["molfile"], sanitize=True, removeHs=False)
    return None


def state_selection(panel: dict, key: str) -> list[int] | None:
    """Return original ONE-based molecule indices, scoped to a state."""
    declared = panel.get(key, {})
    if "molecule_indices" in declared:
        return declared["molecule_indices"]
    molecules = panel.get("molecules", [])
    split = panel.get("reaction_arrow_after")
    if split is not None:
        return (
            list(range(1, int(split) + 1))
            if key == "starting_state"
            else list(range(int(split) + 1, len(molecules) + 1))
        )
    side = "reactant" if key == "starting_state" else "product"
    indices = [i for i, mol in enumerate(molecules, 1) if mol.get("side") == side]
    return indices or None


def resolve_states(panel: dict, next_panel: dict | None) -> tuple[list, list | None]:
    def select(owner: dict, indices: list[int] | None) -> list:
        molecules = owner.get("molecules", [])
        if indices is None:
            indices = list(range(1, len(molecules) + 1))
        if not indices or len(indices) != len(set(indices)):
            raise ValueError("state molecule_indices must be nonempty and unique")
        if any(type(i) is not int or i < 1 or i > len(molecules) for i in indices):
            raise ValueError("state molecule_indices are outside the panel")
        return [(i, molecules[i - 1]) for i in indices]

    start = select(panel, state_selection(panel, "starting_state"))
    end_indices = state_selection(panel, "expected_state")
    if end_indices is not None:
        end = select(panel, end_indices)
        if {i for i, _ in start} & set(end_indices):
            raise ValueError(
                "starting_state and expected_state must select disjoint molecule entries"
            )
        split = panel.get("reaction_arrow_after")
        if split is not None and (
            [i for i, _ in start] != list(range(1, int(split) + 1))
            or end_indices != list(range(int(split) + 1, len(panel["molecules"]) + 1))
        ):
            raise ValueError(
                "state selections disagree with the visible reaction_arrow_after split"
            )
    elif next_panel is not None:
        end = select(next_panel, state_selection(next_panel, "starting_state"))
    else:
        end = None
    return start, end


def _state(entries: list, errors: list[str]) -> dict:
    atoms, bonds, inventory = {}, {}, Counter()
    unsupported = []
    for mol_index, mol_spec in entries:
        mol = parse_molecule(mol_spec)
        if mol is None:
            errors.append(f"molecule {mol_index}: structure cannot be parsed/sanitized")
            continue
        for atom in mol.GetAtoms():
            n = atom.GetAtomMapNum()
            if not n or n in atoms:
                errors.append(
                    f"molecule {mol_index}: missing or duplicate atom map {n} within state"
                )
            h = (
                atom.GetTotalNumHs()
            )  # Excludes separate H vertices; do not count them twice.
            inventory[(atom.GetSymbol(), atom.GetIsotope())] += 1
            inventory[("H", 0)] += h
            atoms[n] = {
                "symbol": atom.GetSymbol(),
                "isotope": atom.GetIsotope(),
                "charge": atom.GetFormalCharge(),
                "hydrogens": h,
                "radicals": atom.GetNumRadicalElectrons(),
                "molecule_index": mol_index,
                "lone_pairs": (
                    Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())
                    - atom.GetFormalCharge()
                    - atom.GetTotalValence()
                    - atom.GetNumRadicalElectrons()
                )
                / 2,
            }
            if atom.GetSymbol() not in {
                "H",
                "B",
                "C",
                "N",
                "O",
                "F",
                "Si",
                "P",
                "S",
                "Cl",
                "Br",
                "I",
            }:
                unsupported.append("metal/unsupported element electron bookkeeping")
            if atom.GetIsAromatic():
                unsupported.append(
                    "aromatic electron replay (supply an explicit localized model)"
                )
        for bond in mol.GetBonds():
            key = tuple(
                sorted(
                    (
                        bond.GetBeginAtom().GetAtomMapNum(),
                        bond.GetEndAtom().GetAtomMapNum(),
                    )
                )
            )
            bonds[key] = bond.GetBondTypeAsDouble()
    return {
        "atoms": atoms,
        "bonds": bonds,
        "inventory": +inventory,
        "unsupported": sorted(set(unsupported)),
    }


def _pair(value: Any) -> tuple[int, int]:
    if (
        not isinstance(value, (list, tuple))
        or len(value) != 2
        or any(type(i) is not int or i <= 0 for i in value)
        or value[0] == value[1]
    ):
        raise ValueError("bond anchors require two distinct positive atom-map integers")
    return tuple(sorted(value))


def _signature(state: dict) -> tuple:
    fields = ("symbol", "isotope", "charge", "hydrogens", "radicals")
    return (
        {k: {f: a[f] for f in fields} for k, a in state["atoms"].items()},
        state["bonds"],
    )


def _changes(before: dict, after: dict) -> dict:
    bonds = [
        {
            "atom_maps": list(k),
            "from": before["bonds"].get(k, 0),
            "to": after["bonds"].get(k, 0),
        }
        for k in sorted(before["bonds"].keys() | after["bonds"].keys())
        if before["bonds"].get(k, 0) != after["bonds"].get(k, 0)
    ]
    changes = {"bonds": bonds}
    for field in ("charge", "radicals"):
        changes[field] = [
            {
                "atom_map": k,
                "from": before["atoms"][k][field],
                "to": after["atoms"][k][field],
            }
            for k in sorted(before["atoms"].keys() & after["atoms"].keys())
            if before["atoms"][k][field] != after["atoms"][k][field]
        ]
    return changes


def _replay(
    panel: dict, before: dict, errors: list[str], unsupported: list[str]
) -> dict:
    atoms, bonds = before["atoms"], before["bonds"]
    bond_delta, charge_delta, radical_delta = Counter(), Counter(), Counter()
    bond_spent, lp_spent = Counter(), Counter()
    for number, arrow in enumerate(panel.get("arrows", []), 1):
        label = f"arrow {number}"
        try:
            sources = [
                k
                for k in (
                    "from_lone_pair_atom_map",
                    "from_bond",
                    "from_atom_map",
                    "from_xy",
                )
                if k in arrow
            ]
            sinks = [k for k in ("to_atom_map", "to_bond", "to_xy") if k in arrow]
            if len(sources) != 1 or len(sinks) != 1:
                raise ValueError("requires exactly one source and one sink")
            kind = arrow.get("kind")
            electrons = arrow.get(
                "electron_count", 1 if kind == "single_electron" else 2
            )
            if (
                kind not in {"electron_pair", "single_electron"}
                or type(electrons) is not int
                or electrons != (1 if kind == "single_electron" else 2)
            ):
                raise ValueError("arrow kind and electron_count disagree")
            if sources[0] == "from_xy" or sinks[0] == "to_xy":
                unsupported.append("free-coordinate electron anchors")
                continue
            src = (
                _pair(arrow["from_bond"])
                if sources[0] == "from_bond"
                else (arrow[sources[0]],)
            )
            dst = (
                _pair(arrow["to_bond"])
                if sinks[0] == "to_bond"
                else (arrow["to_atom_map"],)
            )
            if any(type(n) is not int or n not in atoms for n in (*src, *dst)):
                raise ValueError("anchor atom map absent from starting state")
            if len(src) == 2 and src not in bonds:
                raise ValueError(f"source bond {src} does not exist in starting state")
            for ns, hint in ((src, "from_molecule_index"), (dst, "to_molecule_index")):
                if hint in arrow and any(
                    atoms[n]["molecule_index"] != arrow[hint] for n in ns
                ):
                    raise ValueError(
                        f"{hint} does not refer to the starting-state anchor"
                    )
            if len(src) == 2:
                bond_spent[src] += electrons
            if electrons == 1:
                if len(src) == 2 and len(dst) == 1 and dst[0] in src:
                    bond_delta[src] -= 0.5
                    radical_delta[dst[0]] += 1
                else:
                    unsupported.append(
                        "single-electron transfer other than paired bond homolysis"
                    )
                continue
            if sources[0] == "from_atom_map":
                unsupported.append(
                    "electron-pair source at atom center instead of lone pair"
                )
                continue
            if len(src) == 1:
                donor = src[0]
                lp_spent[donor] += 1
                if len(dst) == 2 and donor not in dst:
                    raise ValueError("lone-pair to_bond must contain the donating atom")
                acceptor = (
                    dst[0] if len(dst) == 1 else next(n for n in dst if n != donor)
                )
                if acceptor == donor:
                    raise ValueError("lone pair cannot attack its own atom")
                bond_delta[_pair([donor, acceptor])] += 1
                charge_delta[donor] += 1
                charge_delta[acceptor] -= 1
            elif len(dst) == 1:
                if dst[0] not in src:
                    raise ValueError(
                        "bond-to-atom transfer must end on a bond endpoint; use to_bond for migration"
                    )
                acceptor = dst[0]
                donor = next(n for n in src if n != acceptor)
                bond_delta[src] -= 1
                charge_delta[donor] += 1
                charge_delta[acceptor] -= 1
            else:
                shared = set(src) & set(dst)
                if len(shared) != 1:
                    raise ValueError(
                        "bond-to-bond pair transfer must share exactly one atom"
                    )
                donor = next(n for n in src if n not in shared)
                acceptor = next(n for n in dst if n not in shared)
                bond_delta[src] -= 1
                bond_delta[dst] += 1
                charge_delta[donor] += 1
                charge_delta[acceptor] -= 1
        except (ValueError, TypeError, KeyError) as exc:
            errors.append(f"{label}: {exc}")
    for k, count in bond_spent.items():
        if count > 2 * bonds[k]:
            errors.append(
                f"source bond {k}: arrows spend more electrons than are present"
            )
    for k, count in lp_spent.items():
        if count > atoms[k]["lone_pairs"]:
            errors.append(f"atom {k}: arrows spend more lone pairs than are present")
    predicted_bonds = {
        k: bonds.get(k, 0) + bond_delta[k] for k in bonds.keys() | bond_delta.keys()
    }
    if any(order < 0 or order != int(order) for order in predicted_bonds.values()):
        errors.append(
            "electron moves leave a negative/fractional bond order; homolysis requires paired fishhooks"
        )
    predicted_atoms = {
        k: {
            **a,
            "charge": a["charge"] + charge_delta[k],
            "radicals": a["radicals"] + radical_delta[k],
        }
        for k, a in atoms.items()
    }
    return {
        "atoms": predicted_atoms,
        "bonds": {k: v for k, v in predicted_bonds.items() if v},
    }


def _check_edits(panel: dict, before: dict, after: dict, errors: list[str]) -> None:
    """Compare declarations to independently computed state differences."""
    if not panel.get("graph_edits"):
        return
    declared_bonds, declared_atoms = {}, {"charge": {}, "radicals": {}}
    for number, edit in enumerate(panel["graph_edits"], 1):
        try:
            kind = edit.get("type", edit.get("operation"))
            if kind in {"form_bond", "break_bond", "change_bond"}:
                key = _pair(edit["atom_maps"])
                old = before["bonds"].get(key, 0)
                new = 0 if kind == "break_bond" else edit.get("to", edit.get("order"))
                if new is None or type(new) not in (int, float):
                    raise ValueError("bond edit requires numeric order/to")
                if (
                    key in declared_bonds
                    or (kind == "form_bond" and old)
                    or (kind != "form_bond" and not old)
                ):
                    raise ValueError(
                        "bond edit duplicates or contradicts starting bond"
                    )
                if "from" in edit and edit["from"] != old:
                    raise ValueError(
                        "bond edit from value disagrees with starting state"
                    )
                declared_bonds[key] = (old, new)
            elif kind in {"formal_charge", "radical"}:
                field = "charge" if kind == "formal_charge" else "radicals"
                key = edit["atom_map"]
                if (
                    key in declared_atoms[field]
                    or before["atoms"][key][field] != edit["from"]
                ):
                    raise ValueError("atom edit duplicates or has wrong from value")
                declared_atoms[field][key] = (edit["from"], edit["to"])
            elif kind not in {
                "partial_charge",
                "proton_transfer",
                "counterion",
                "resonance",
            }:
                raise ValueError(f"unknown graph edit {kind!r}")
            if "source_arrow" in edit and (
                type(edit["source_arrow"]) is not int
                or not 1 <= edit["source_arrow"] <= len(panel.get("arrows", []))
            ):
                raise ValueError("source_arrow is outside the step's arrows")
            if "source_arrow" in edit and kind in {
                "form_bond",
                "break_bond",
                "change_bond",
            }:
                arrow_errors, arrow_unsupported = [], []
                contribution = _replay(
                    {"arrows": [panel["arrows"][edit["source_arrow"] - 1]]},
                    before,
                    arrow_errors,
                    arrow_unsupported,
                )
                delta = contribution["bonds"].get(key, 0) - old
                if not arrow_unsupported and (not delta or delta * (new - old) <= 0):
                    raise ValueError(
                        "source_arrow does not contribute the declared bond change"
                    )
        except (ValueError, TypeError, KeyError) as exc:
            errors.append(f"graph_edit {number}: {exc}")
    changes = _changes(before, after)
    expected_bonds = {
        tuple(c["atom_maps"]): (c["from"], c["to"]) for c in changes["bonds"]
    }
    if declared_bonds != expected_bonds:
        errors.append("graph_edits do not match the actual bond changes between states")
    for field in ("charge", "radicals"):
        expected = {c["atom_map"]: (c["from"], c["to"]) for c in changes[field]}
        if declared_atoms[field] != expected:
            errors.append(
                f"graph_edits do not match the actual {field} changes between states"
            )


def validate_semantics(spec: dict[str, Any]) -> dict[str, Any]:
    checks, all_errors, all_unsupported = [], [], []
    panels = spec.get("panels", [])
    previous_after = None
    for index, panel in enumerate(panels, 1):
        if not panel.get("arrows") and not panel.get("graph_edits"):
            if str(spec.get("spec_version")) == "2.1":
                error = "MechanismSpec 2.1 requires electron moves and graph edits in every step panel"
                checks.append(
                    {
                        "panel": index,
                        "status": "invalid",
                        "errors": [error],
                        "unsupported": [],
                    }
                )
                all_errors.append(f"Panel {index}: {error}")
            continue
        errors, unsupported = [], []
        result = {
            "panel": index,
            "status": "incomplete",
            "errors": errors,
            "unsupported": unsupported,
        }
        try:
            if str(spec.get("spec_version")) == "2.1" and any(
                "molecule_indices" not in panel.get(k, {})
                for k in ("starting_state", "expected_state")
            ):
                errors.append(
                    "MechanismSpec 2.1 requires explicit starting_state and expected_state molecule_indices"
                )
            start, end = resolve_states(
                panel, panels[index] if index < len(panels) else None
            )
            before = _state(start, errors)
            if previous_after is not None and _signature(previous_after) != _signature(
                before
            ):
                errors.append(
                    "previous expected state does not match this starting state; include all participants consistently"
                )
            if end is None:
                unsupported.append(
                    "explicit expected state is missing; electron moves cannot be checked against a product"
                )
                previous_after = None
            else:
                after = _state(end, errors)
                for key, actual in (
                    ("starting_state", before),
                    ("expected_state", after),
                ):
                    declared_smiles = panel.get(key, {}).get("reaction_smiles")
                    if declared_smiles and _signature(
                        _state([(1, {"smiles": declared_smiles})], errors)
                    ) != _signature(actual):
                        errors.append(
                            f"{key}: declared reaction_smiles disagrees with selected molecules"
                        )
                if str(spec.get("spec_version")) == "2.1" and not panel.get(
                    "graph_edits"
                ):
                    errors.append(
                        "MechanismSpec 2.1 requires explicit graph_edits for each elementary step"
                    )
                previous_after = after
                unsupported.extend(before["unsupported"] + after["unsupported"])
                left, right = before["atoms"], after["atoms"]
                if left.keys() != right.keys():
                    errors.append("atom maps are not conserved between states")
                for n in left.keys() & right.keys():
                    if (left[n]["symbol"], left[n]["isotope"]) != (
                        right[n]["symbol"],
                        right[n]["isotope"],
                    ):
                        errors.append(f"atom {n}: element/isotope identity changed")
                    if left[n]["hydrogens"] != right[n]["hydrogens"]:
                        errors.append(
                            f"atom {n}: implicit/bracket H count changed; transfer an explicit mapped H vertex"
                        )
                if before["inventory"] != after["inventory"]:
                    errors.append(
                        "element/isotope inventory (including hydrogens) is not conserved"
                    )
                if sum(a["charge"] for a in left.values()) != sum(
                    a["charge"] for a in right.values()
                ):
                    errors.append(
                        "total formal charge is not conserved; prose cannot replace missing species"
                    )
                result["observed_changes"] = _changes(before, after)
                if not errors and not unsupported:
                    predicted = _replay(panel, before, errors, unsupported)
                    result["arrow_predicted_changes"] = _changes(before, predicted)
                    if not unsupported and _changes(before, predicted) != _changes(
                        before, after
                    ):
                        errors.append(
                            "electron moves do not generate the expected state's bonds/charges/radicals"
                        )
                _check_edits(panel, before, after, errors)
                result["status"] = (
                    "invalid" if errors else "unsupported" if unsupported else "valid"
                )
        except (ValueError, TypeError, KeyError, IndexError, AttributeError) as exc:
            errors.append(f"malformed state/arrow data: {exc}")
        if errors:
            result["status"] = "invalid"
        result["unsupported"] = sorted(set(unsupported))
        checks.append(result)
        all_errors.extend(f"Panel {index}: {e}" for e in errors)
        all_unsupported.extend(f"Panel {index}: {e}" for e in result["unsupported"])
    return {
        "status": "invalid"
        if all_errors
        else "incomplete"
        if not checks or any(c["status"] == "incomplete" for c in checks)
        else "unsupported"
        if all_unsupported
        else "valid",
        "scope": "main-group electron-pair graph bookkeeping and paired homolysis; not kinetic, energetic, stereo or experimental proof",
        "steps": checks,
        "errors": all_errors,
        "unsupported": all_unsupported,
    }
