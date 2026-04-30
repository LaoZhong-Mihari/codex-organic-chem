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


def _element_counter(smiles_values: list[str]) -> Counter:
    counter: Counter = Counter()
    if not RDKIT_AVAILABLE:
        return counter
    for smi in smiles_values:
        mol = _mol_from_smiles(smi)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            counter[atom.GetSymbol()] += 1
    return counter


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
    reactant_elements = _element_counter(parts.reactants)
    product_elements = _element_counter(parts.products)
    element_delta = {
        key: product_elements.get(key, 0) - reactant_elements.get(key, 0)
        for key in sorted(set(reactant_elements) | set(product_elements))
        if product_elements.get(key, 0) != reactant_elements.get(key, 0)
    }
    if element_delta:
        warnings.append("Reactant/product heavy-atom balance differs; reagents, leaving groups, or stoichiometry may be omitted.")
    classes = infer_reaction_classes(parts)
    conditions = []
    for cls in classes:
        conditions.extend(CONDITION_FAMILIES.get(cls, []))
    svg, draw_warnings = reaction_to_svg(reaction_smiles)
    warnings.extend(draw_warnings)
    analysis: dict[str, Any] = {
        "mode": mode,
        "reaction_classes": classes,
        "element_delta_without_reagents": dict(element_delta),
        "tool_facts": {
            "parsed_reactant_count": len(parts.reactants),
            "parsed_reagent_count": len(parts.reagents),
            "parsed_product_count": len(parts.products),
            "rdkit_reaction_svg_available": bool(svg),
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

