from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any

from .models import ReactionRecord
from .rdkit_tools import RDKIT_AVAILABLE, detect_functional_groups, normalize_structure, reaction_to_svg
from .utils import split_nonempty

if RDKIT_AVAILABLE:  # pragma: no cover - imported through service tests
    from rdkit import Chem
else:  # pragma: no cover
    Chem = None


@dataclass
class ReactionParts:
    reactants: list[str]
    reagents: list[str]
    products: list[str]


def split_reaction_smiles(reaction_smiles: str) -> ReactionParts:
    pieces = reaction_smiles.split(">")
    if len(pieces) == 2:
        left, right = pieces
        middle = ""
    elif len(pieces) == 3:
        left, middle, right = pieces
    else:
        return ReactionParts([], [], [])
    return ReactionParts(
        reactants=split_nonempty(left),
        reagents=split_nonempty(middle),
        products=split_nonempty(right),
    )


def _records_for(smiles_values: list[str], role: str) -> list[dict]:
    return [normalize_structure(smiles=value, source=role).to_dict() for value in smiles_values]


def _mol_from_smiles(smiles: str) -> Any:
    if not RDKIT_AVAILABLE:
        return None
    return Chem.MolFromSmiles(smiles)


def _atom_counter(smiles_values: list[str]) -> Counter:
    """Count every atom, hydrogens included.

    GetAtoms() on a SMILES-parsed mol does not yield implicit hydrogens, so
    counting symbols alone calls CC(=O)C >> CC(O)C balanced even though the
    product gained two H. GetTotalNumHs() restores them.
    """
    counter: Counter = Counter()
    if not RDKIT_AVAILABLE:
        return counter
    for smi in smiles_values:
        mol = _mol_from_smiles(smi)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            counter[atom.GetSymbol()] += 1
            counter["H"] += atom.GetTotalNumHs()
    return +counter


def _charge_total(smiles_values: list[str]) -> int:
    total = 0
    if not RDKIT_AVAILABLE:
        return total
    for smi in smiles_values:
        mol = _mol_from_smiles(smi)
        if mol is None:
            continue
        total += sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    return total


def _counter_delta(left: Counter, right: Counter) -> dict[str, int]:
    return {
        key: right.get(key, 0) - left.get(key, 0)
        for key in sorted(set(left) | set(right))
        if right.get(key, 0) != left.get(key, 0)
    }


def _formula(counter: Counter) -> str:
    if not counter:
        return ""
    pieces = []
    for symbol in ("C", "H"):
        count = counter.get(symbol, 0)
        if count:
            pieces.append(f"{symbol}{count if count > 1 else ''}")
    for symbol in sorted(key for key in counter if key not in {"C", "H"} and counter[key]):
        count = counter[symbol]
        pieces.append(f"{symbol}{count if count > 1 else ''}")
    return "".join(pieces)


def _delta_text(delta: dict[str, int]) -> str:
    return ", ".join(f"{symbol}{value:+d}" for symbol, value in delta.items())


def mass_charge_balance(parts: ReactionParts) -> dict:
    """Full mass and charge accounting for a written reaction.

    Reported deltas are product minus reactant. Two variants are computed
    because reagents written between the '>' separators may be catalysts
    (excluding them balances) or stoichiometric partners (including them
    balances); the reaction is only flagged when neither variant closes.
    """
    if not RDKIT_AVAILABLE:
        return {
            "status": "unavailable",
            "messages": ["RDKit is unavailable; mass and charge balance were not checked."],
        }
    reactants = _atom_counter(parts.reactants)
    reagents = _atom_counter(parts.reagents)
    products = _atom_counter(parts.products)
    without_reagents = _counter_delta(reactants, products)
    with_reagents = _counter_delta(reactants + reagents, products)

    def magnitude(delta: dict[str, int]) -> int:
        return sum(abs(value) for value in delta.values())

    # Charge must use the same reagent convention as the atom count, or a
    # catalytic proton written in the reagent slot reads as a lost charge.
    reactant_charge = _charge_total(parts.reactants)
    reagent_charge = _charge_total(parts.reagents)
    product_charge = _charge_total(parts.products)
    if not without_reagents and with_reagents:
        include_reagents = False
    elif not with_reagents and without_reagents:
        include_reagents = True
    else:
        include_reagents = magnitude(with_reagents) < magnitude(without_reagents)
    charge_left = reactant_charge + (reagent_charge if include_reagents else 0)
    charge_delta = product_charge - charge_left
    if charge_delta and not include_reagents and product_charge == reactant_charge + reagent_charge:
        # Reagents are stoichiometric for charge even though atoms balance
        # without them (e.g. a counterion carried through unchanged).
        charge_left = reactant_charge + reagent_charge
        charge_delta = 0

    residual = with_reagents if include_reagents else without_reagents
    atoms_balanced = not without_reagents or not with_reagents
    messages: list[str] = []
    if not atoms_balanced:
        if set(residual) == {"H"}:
            gained = residual["H"]
            if gained > 0:
                messages.append(
                    f"Hydrogen count changes by {gained:+d} with no heavy-atom change: the written reaction implies a "
                    "net reduction or protonation, so a hydride/proton source is missing."
                )
            else:
                messages.append(
                    f"Hydrogen count changes by {gained:+d} with no heavy-atom change: the written reaction implies a "
                    "net oxidation or deprotonation, so an oxidant/base is missing."
                )
        else:
            messages.append(
                f"Atom balance does not close ({_delta_text(residual)}); reagents, leaving groups, solvent-derived "
                "atoms, or stoichiometry are missing from the written reaction."
            )
    if charge_delta:
        messages.append(
            f"Formal charge is not conserved ({charge_left:+d} -> {product_charge:+d}); a counterion, proton, or "
            "electron transfer is missing from the written reaction."
        )
    if atoms_balanced and not charge_delta:
        status = "balanced"
    else:
        status = "unbalanced"
    return {
        "status": status,
        "reactant_formula": _formula(reactants),
        "reagent_formula": _formula(reagents),
        "product_formula": _formula(products),
        "element_delta_without_reagents": dict(without_reagents),
        "element_delta_with_reagents": dict(with_reagents),
        "hydrogen_delta": residual.get("H", 0) if not atoms_balanced else 0,
        "reagents_counted_in_balance": include_reagents,
        "reactant_charge": charge_left,
        "product_charge": product_charge,
        "charge_delta": charge_delta,
        "messages": messages,
    }


def _functional_groups(smiles_values: list[str]) -> set[str]:
    groups: set[str] = set()
    if not RDKIT_AVAILABLE:
        return groups
    for smi in smiles_values:
        mol = _mol_from_smiles(smi)
        if mol is not None:
            groups.update(detect_functional_groups(mol))
    return groups


def _has_substructure(smiles_values: list[str], smarts: str) -> bool:
    if not RDKIT_AVAILABLE:
        return False
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return False
    return any((mol := _mol_from_smiles(smi)) is not None and mol.HasSubstructMatch(patt) for smi in smiles_values)


def _mapped_bonds(smiles_values: list[str]) -> dict[tuple[int, int], str]:
    bonds: dict[tuple[int, int], str] = {}
    if not RDKIT_AVAILABLE:
        return bonds
    for smi in smiles_values:
        mol = _mol_from_smiles(smi)
        if mol is None:
            continue
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom().GetAtomMapNum()
            b = bond.GetEndAtom().GetAtomMapNum()
            if a and b:
                bonds[tuple(sorted((a, b)))] = str(bond.GetBondType())
    return bonds


def atom_mapping_delta(parts: ReactionParts) -> dict:
    reactant_bonds = _mapped_bonds(parts.reactants)
    product_bonds = _mapped_bonds(parts.products)
    if not reactant_bonds and not product_bonds:
        return {
            "status": "unmapped",
            "formed_bonds": [],
            "broken_bonds": [],
            "changed_bonds": [],
            "message": "No atom-map numbers detected; bond-change analysis is heuristic only.",
        }
    formed = [f"{a}-{b}:{product_bonds[(a, b)]}" for a, b in sorted(product_bonds.keys() - reactant_bonds.keys())]
    broken = [f"{a}-{b}:{reactant_bonds[(a, b)]}" for a, b in sorted(reactant_bonds.keys() - product_bonds.keys())]
    changed = []
    for key in sorted(reactant_bonds.keys() & product_bonds.keys()):
        if reactant_bonds[key] != product_bonds[key]:
            changed.append(f"{key[0]}-{key[1]}:{reactant_bonds[key]}->{product_bonds[key]}")
    return {
        "status": "mapped",
        "formed_bonds": formed,
        "broken_bonds": broken,
        "changed_bonds": changed,
    }


def infer_reaction_classes(parts: ReactionParts) -> list[str]:
    classes: list[str] = []
    reactant_groups = _functional_groups(parts.reactants)
    product_groups = _functional_groups(parts.products)
    all_reactants = parts.reactants + parts.reagents
    if _has_substructure(parts.reactants, "[CX4][Br,Cl,I]") and (
        _has_substructure(all_reactants, "[O-]") or _has_substructure(all_reactants, "[N-]") or _has_substructure(all_reactants, "[S-]")
    ):
        classes.append("possible_SN2_substitution")
    if "carboxylic_acid" in reactant_groups and "alcohol" in reactant_groups and "ester" in product_groups:
        classes.append("possible_fischer_esterification")
    if "aryl_halide" in reactant_groups and "boronic_acid" in reactant_groups:
        classes.append("possible_suzuki_coupling")
    if {"aldehyde", "ketone"} & reactant_groups and ("alcohol" in product_groups or "alkene" in product_groups):
        classes.append("possible_aldol_family")
    if _has_substructure(parts.reactants, "c1ccccc1") and (
        _has_substructure(all_reactants, "[N+](=O)[O-]") or _has_substructure(all_reactants, "BrBr") or _has_substructure(all_reactants, "ClCl")
    ):
        classes.append("possible_electrophilic_aromatic_substitution")
    if not classes:
        classes.append("unclassified_rule_based")
    return classes


CONDITION_FAMILIES = {
    "possible_SN2_substitution": [
        "Polar aprotic solvent family; choose nucleophile/base strength against substrate sensitivity.",
        "Check leaving group, steric hindrance, beta-elimination risk, and competing solvolysis.",
    ],
    "possible_fischer_esterification": [
        "Acid catalysis with removal/excess of alcohol or water management.",
        "For sensitive substrates, compare coupling reagents or acid chloride alternatives in literature.",
    ],
    "possible_suzuki_coupling": [
        "Pd catalyst, base, and mixed aqueous/organic solvent family.",
        "Check oxidative-addition difficulty, boronic acid stability, protodeboronation, and ligand precedents.",
    ],
    "possible_aldol_family": [
        "Base or Lewis-acid/enamine conditions depending on desired regio- and stereocontrol.",
        "Check self-condensation, dehydration, and enolate geometry assumptions.",
    ],
    "possible_electrophilic_aromatic_substitution": [
        "Electrophile generation plus acid/Lewis-acid activation family.",
        "Check directing effects, regioisomer mixtures, overreaction, and functional-group tolerance.",
    ],
}


def retrosynthesis_hints(product_smiles: list[str]) -> list[dict]:
    hints: list[dict] = []
    groups = _functional_groups(product_smiles)
    if "ester" in groups:
        hints.append(
            {
                "disconnection": "acyl-O bond of ester",
                "strategy": "carboxylic acid derivative plus alcohol; compare Fischer, acid chloride, anhydride, or coupling reagent options.",
                "risk": "acid/base sensitivity, transesterification, and chemoselectivity need literature support.",
            }
        )
    if "amide" in groups:
        hints.append(
            {
                "disconnection": "acyl-N bond of amide",
                "strategy": "carboxylic acid derivative plus amine or peptide-coupling style activation.",
                "risk": "amine nucleophilicity, racemization, and overacylation need review.",
            }
        )
    if _has_substructure(product_smiles, "c-c"):
        hints.append(
            {
                "disconnection": "aryl-aryl or aryl-alkyl C-C bond",
                "strategy": "cross-coupling or directed functionalization candidate.",
                "risk": "regioselectivity and catalyst/ligand precedents are decisive.",
            }
        )
    if "alcohol" in groups and ({"aldehyde", "ketone"} & groups):
        hints.append(
            {
                "disconnection": "C-C bond adjacent to hydroxy/carbonyl motif",
                "strategy": "aldol, Reformatsky, Grignard, or related carbonyl addition family.",
                "risk": "diastereoselectivity and dehydration may dominate.",
            }
        )
    if not hints:
        hints.append(
            {
                "disconnection": "no high-confidence rule-based disconnection",
                "strategy": "run literature search/retrosynthesis engine and inspect functional groups manually.",
                "risk": "MVP rule library is intentionally conservative.",
            }
        )
    return hints


def analyze_reaction(reaction_smiles: str, mode: str = "sanity_check") -> ReactionRecord:
    parts = split_reaction_smiles(reaction_smiles)
    warnings: list[str] = []
    if not parts.reactants or not parts.products:
        warnings.append("Reaction SMILES should contain reactants and products separated by '>' or '>>'.")
    reactants = _records_for(parts.reactants, "reactant")
    reagents = _records_for(parts.reagents, "reagent")
    products = _records_for(parts.products, "product")
    atom_delta = atom_mapping_delta(parts)
    if atom_delta.get("status") == "unmapped":
        warnings.append(atom_delta["message"])
    balance = mass_charge_balance(parts)
    warnings.extend(balance.get("messages", []))
    classes = infer_reaction_classes(parts)
    conditions = []
    for cls in classes:
        conditions.extend(CONDITION_FAMILIES.get(cls, []))
    svg, draw_warnings = reaction_to_svg(reaction_smiles)
    warnings.extend(draw_warnings)
    analysis: dict[str, Any] = {
        "mode": mode,
        "reaction_classes": classes,
        "mass_charge_balance": balance,
        "element_delta_without_reagents": balance.get("element_delta_without_reagents", {}),
        "tool_facts": {
            "parsed_reactant_count": len(parts.reactants),
            "parsed_reagent_count": len(parts.reagents),
            "parsed_product_count": len(parts.products),
            "rdkit_reaction_svg_available": bool(svg),
            "mass_charge_balance_status": balance.get("status"),
        },
        "rule_inferences": {
            "plausibility_checks": [
                "Check atom economy, product regiochemistry, stereochemistry, and hidden proton/solvent transfers.",
                "Use literature or a dedicated retrosynthesis engine for route selection before laboratory use.",
            ],
        },
        "llm_assumptions": [
            "No paid reaction database was queried.",
            "Rule-based labels are suggestions, not proof of feasibility.",
        ],
    }
    if mode == "retro":
        analysis["retrosynthesis_hints"] = retrosynthesis_hints(parts.products)
    elif mode == "conditions":
        analysis["condition_families"] = conditions
    elif mode == "forward":
        analysis["forward_notes"] = [
            "The MVP does not predict novel products with a learned model; it checks whether the supplied product is chemically plausible.",
            "For side products, inspect leaving groups, acid/base sites, competing nucleophiles, and redox-sensitive groups.",
        ]
    confidence = 0.72
    if "unclassified_rule_based" in classes:
        confidence = 0.45
    if warnings:
        confidence = min(confidence, 0.62)
    if balance.get("status") == "unbalanced":
        # A reaction that does not conserve atoms or charge is not merely
        # under-annotated: as written it is wrong, so cap harder than for a
        # missing atom map.
        confidence = min(confidence, 0.3)
    return ReactionRecord(
        source="reaction_smiles",
        reaction_smiles=reaction_smiles,
        reactants=reactants,
        reagents=reagents,
        products=products,
        atom_mapping=atom_delta,
        conditions=conditions,
        evidence=[
            "RDKit parsing and functional-group heuristics",
            "Local rule library; no commercial reaction database was consulted",
        ],
        confidence=confidence,
        warnings=warnings,
        analysis=analysis,
    )

