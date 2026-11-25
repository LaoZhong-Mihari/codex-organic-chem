"""Stepwise intermediate-sequence validation for reaction mechanisms.

Unlike ``mechanism_semantics`` (which requires a full atom-mapped curved-arrow
``MechanismSpec`` 2.1 for figure rendering), :func:`mechanism_steps_check`
audits a plain sequence of proposed intermediates—the form chemists and models
write while deriving a mechanism step by step—and flags chemical bookkeeping
failures before the answer is finalized:

* unparseable or valence-invalid SMILES at any step;
* odd-electron radicals appearing without a radical/SET/homolytic step, or
  unexplained spin-parity jumps;
* ``proton_transfer`` steps that alter the heavy-atom connectivity skeleton or
  change heavy-atom counts;
* elements appearing in intermediates that are absent from reactants, declared
  products, and reagents named in the conditions;
* consecutive step heavy-atom/charge deltas that cannot be explained by adding a
  reactant or gaining/losing common small molecules / ions;
* final intermediate failing to match the target product—distinguishing formula
  mismatch, stereoisomer mismatch, tautomers, and **constitutional isomers**
  (the classic sign of an incomplete rearrangement or ring-opening/reclosure
  cascade).
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolHash

from .structure_compare import StructureError, compare_structures, formula_counter, hill_formula

RADICAL_KEYWORDS = (
    "radical",
    "homoly",
    "single_electron",
    "set",
    "photoredox",
    "atom_transfer",
    "hat",
    "xat",
    "birch",
    "minisci",
)

PROTON_TRANSFER_KEYWORDS = (
    "proton_transfer",
    "acid_base",
    "protonation",
    "deprotonation",
    "tautomer",
)

ISOMERIZATION_KEYWORDS = (
    "rearrangement",
    "1,2_shift",
    "1,2-shift",
    "electrocyclic",
    "sigmatropic",
    "ring_opening_electrocyclic",
    "ring_closing_electrocyclic",
    "cope",
    "claisen",
    "pinacol",
    "wagner",
)

# Common neutral and ionic small species frequently omitted as spectators/byproducts between steps.
COMMON_SMALL_SPECIES = (
    "O",           # H2O
    "[OH-]",       # HO-
    "[OH3+]",      # H3O+
    "CO",          # MeOH
    "[O-]C",       # MeO-
    "CCO",         # EtOH
    "[O-]CC",      # EtO-
    "CC(C)O",      # iPrOH
    "CC(C)(C)O",   # tBuOH
    "N#N",         # N2
    "O=C=O",       # CO2
    "[C-]#[O+]",   # CO
    "N",           # NH3
    "[NH4+]",      # NH4+
    "CN(C)C",      # Me3N
    "CCN(CC)CC",   # Et3N
    "Cl",          # HCl
    "[Cl-]",       # Cl-
    "Br",          # HBr
    "[Br-]",       # Br-
    "I",           # HI
    "[I-]",        # I-
    "F",           # HF
    "[F-]",        # F-
    "ClCl",        # Cl2
    "BrBr",        # Br2
    "II",          # I2
    "CC(=O)O",     # AcOH
    "CC(=O)[O-]",  # AcO-
    "OC(=O)C(F)(F)F",     # TFA
    "[O-]C(=O)C(F)(F)F",  # TFA-
    "O=S=O",       # SO2
    "O=S(=O)(O)O", # H2SO4
    "CS(C)=O",     # DMSO
    "CS(C)C",      # Me2S
    "C=O",         # CH2O
    "CC(C)=O",     # acetone
    "S",           # H2S
    "[SH-]",       # HS-
    "C#N",         # HCN
    "[C-]#N",      # CN-
    "O=P(c1ccccc1)(c1ccccc1)c1ccccc1",  # Ph3PO
    "P(c1ccccc1)(c1ccccc1)c1ccccc1",    # Ph3P
)

REAGENT_ELEMENTS: dict[str, str] = {
    "PPh3": "CPO",
    "Ph3P": "CP",
    "P(OEt)3": "CPO",
    "P(OMe)3": "CPO",
    "POCl3": "POCl",
    "PCl5": "PCl",
    "PCl3": "PCl",
    "PBr3": "PBr",
    "SOCl2": "SOCl",
    "SO2Cl2": "SOCl",
    "COCl2": "COCl",
    "(COCl)2": "COCl",
    "oxalyl": "COCl",
    "Na": "Na",
    "NaOH": "NaO",
    "NaOMe": "NaCO",
    "NaOEt": "NaCO",
    "NaH": "Na",
    "NaBH4": "NaB",
    "NaCNBH3": "NaBCN",
    "NaN3": "NaN",
    "K": "K",
    "KOH": "KO",
    "KOtBu": "KCO",
    "K2CO3": "KCO",
    "Cs2CO3": "CsCO",
    "Li": "Li",
    "LiOH": "LiO",
    "LDA": "LiCN",
    "LiHMDS": "LiCSiN",
    "NaHMDS": "NaCSiN",
    "KHMDS": "KCSiN",
    "BuLi": "LiC",
    "LiAlH4": "LiAl",
    "LAH": "LiAl",
    "DIBAL": "AlC",
    "Mg": "Mg",
    "MgBr": "MgBr",
    "Zn": "Zn",
    "Cu": "Cu",
    "Pd": "Pd",
    "Pt": "Pt",
    "Rh": "Rh",
    "Ru": "Ru",
    "Ir": "Ir",
    "Au": "Au",
    "Ag": "Ag",
    "Hg": "Hg",
    "Ti": "Ti",
    "Sn": "Sn",
    "Bu3SnH": "SnC",
    "SmI2": "SmI",
    "Si": "Si",
    "TMS": "CSi",
    "TBS": "CSi",
    "TBDMS": "CSi",
    "TBDPS": "CSi",
    "TIPS": "CSi",
    "TBAF": "CNF",
    "B": "B",
    "BH3": "B",
    "B2H6": "B",
    "9-BBN": "BC",
    "Bpin": "BCO",
    "BF3": "BF",
    "BBr3": "BBr",
    "BCl3": "BCl",
    "AlCl3": "AlCl",
    "FeCl3": "FeCl",
    "FeBr3": "FeBr",
    "Cl2": "Cl",
    "HCl": "Cl",
    "NCS": "CNClO",
    "Br2": "Br",
    "HBr": "Br",
    "NBS": "CNBrO",
    "I2": "I",
    "HI": "I",
    "NIS": "CNIO",
    "Ts": "CSO",
    "TsOH": "CSO",
    "Ms": "CSO",
    "MsCl": "CSOCl",
    "Tf": "CSFO",
    "TFA": "CFO",
    "TFAA": "CFO",
    "Selectfluor": "CNFOClB",
    "NFSI": "CNSFO",
    "H2S": "S",
    "S8": "S",
    "sulfur": "S",
    "Lawesson": "CPSO",
    "Se": "Se",
    "SeO2": "SeO",
    "H2O": "O",
    "water": "O",
    "aq": "O",
    "OH": "O",
    "ROH": "CO",
    "MeOH": "CO",
    "EtOH": "CO",
    "iPrOH": "CO",
    "tBuOH": "CO",
    "alcohol": "CO",
    "H2SO4": "SO",
    "HNO3": "NO",
    "HNO2": "NO",
    "NaNO2": "NaNO",
    "NO2": "NO",
    "AcOH": "CO",
    "Ac2O": "CO",
    "acetic": "CO",
    "formic": "CO",
    "HCO2H": "CO",
    "CH2O": "CO",
    "formaldehyde": "CO",
    "O2": "O",
    "O3": "O",
    "ozone": "O",
    "H2O2": "O",
    "peroxide": "O",
    "TBHP": "CO",
    "mCPBA": "COCl",
    "DMDO": "CO",
    "OsO4": "OsO",
    "KMnO4": "KMnO",
    "CrO3": "CrO",
    "PCC": "CrCNOCl",
    "DMP": "CIO",
    "IBX": "CIO",
    "Swern": "CSOClN",
    "DMF": "CNO",
    "DMSO": "CSO",
    "THF": "CO",
    "dioxane": "CO",
    "CH2Cl2": "CCl",
    "DCM": "CCl",
    "CHCl3": "CCl",
    "MeCN": "CN",
    "CO2": "CO",
    "NH3": "N",
    "NH4": "N",
    "NH2OH": "NO",
    "N2H4": "N",
    "hydrazine": "N",
    "amine": "CN",
    "morpholine": "CNO",
    "piperidine": "CN",
    "pyridine": "CN",
    "Et3N": "CN",
    "DIPEA": "CN",
    "DBU": "CN",
    "DMAP": "CN",
    "ammonia": "N",
    "azide": "N",
    "DIAD": "CNO",
    "DEAD": "CNO",
    "AIBN": "CN",
    "BPO": "CO",
}


def _parse_species_list(raw: str | list[str] | None) -> tuple[list[str], list[Chem.Mol], list[str]]:
    if raw is None:
        return [], [], []
    if isinstance(raw, str):
        # Keep multi-fragment SMILES intact if passed as list, or split dot-separated strings if top-level
        tokens = [part.strip() for part in raw.split(".") if part.strip()]
    else:
        tokens = [str(item).strip() for item in raw if str(item).strip()]
    mols: list[Chem.Mol] = []
    errors: list[str] = []
    for token in tokens:
        mol = Chem.MolFromSmiles(token)
        if mol is None:
            errors.append(token)
        else:
            mols.append(mol)
    return tokens, mols, errors


def _fragments(mol: Chem.Mol | None) -> list[Chem.Mol]:
    return list(Chem.GetMolFrags(mol, asMols=True)) if mol is not None else []


def _heavy_counter(mol: Chem.Mol) -> Counter[str]:
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def _skeleton_hash(mol: Chem.Mol) -> str:
    copy_mol = Chem.Mol(mol)
    for atom in copy_mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)
    return rdMolHash.MolHash(copy_mol, rdMolHash.HashFunction.ElementGraph)


def _nostereo_smiles(mol: Chem.Mol) -> str:
    copy_mol = Chem.Mol(mol)
    for atom in copy_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.RemoveStereochemistry(copy_mol)
    return Chem.MolToSmiles(copy_mol)


def _odd_radical_atoms(mol: Chem.Mol) -> list[dict[str, Any]]:
    """Atoms with an odd number of unpaired electrons (excluding singlet carbenes/nitrenes [C]/[N])."""
    out: list[dict[str, Any]] = []
    for atom in mol.GetAtoms():
        rad = atom.GetNumRadicalElectrons()
        if rad % 2 == 1:
            out.append({"atom_index": atom.GetIdx(), "symbol": atom.GetSymbol(), "radical_electrons": rad})
    return out


def _sub_counter(a: Counter[str], b: Counter[str]) -> dict[str, int]:
    keys = sorted(set(a) | set(b))
    return {k: a.get(k, 0) - b.get(k, 0) for k in keys if a.get(k, 0) != b.get(k, 0)}


def _delta_str(diff: dict[str, int]) -> str:
    if not diff:
        return "0"
    order = (["C"] if "C" in diff else []) + (["H"] if "H" in diff else []) + sorted(
        k for k in diff if k not in {"C", "H"}
    )
    return " ".join(f"{diff[k]:+d}{k}" for k in order)


def _reagent_elements_from_conditions(conditions: str) -> set[str]:
    found: set[str] = set()
    if not conditions:
        return found
    lower_cond = conditions.lower()
    for name, elements in REAGENT_ELEMENTS.items():
        if name.lower() in lower_cond if len(name) >= 3 else re.search(rf"\b{re.escape(name)}\b", conditions):
            found |= set(re.findall(r"[A-Z][a-z]?", elements))
    for token in re.split(r"[,;\s/]+", conditions):
        cleaned = token.strip("()[]{}")
        if not cleaned or cleaned.isalpha() and len(cleaned) == 1 and cleaned not in {"S", "P", "B", "I", "F", "K"}:
            continue
        mol = Chem.MolFromSmiles(cleaned)
        if mol is not None and mol.GetNumHeavyAtoms() > 0:
            found |= set(_heavy_counter(mol))
    return found


_SMALL_HEAVY_UNITS: list[Counter[str]] | None = None


def _small_heavy_units() -> list[Counter[str]]:
    global _SMALL_HEAVY_UNITS
    if _SMALL_HEAVY_UNITS is None:
        units: list[Counter[str]] = []
        seen: set[tuple[tuple[str, int], ...]] = set()
        for smi in COMMON_SMALL_SPECIES:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            hc = _heavy_counter(mol)
            key = tuple(sorted(hc.items()))
            if key and key not in seen:
                seen.add(key)
                units.append(hc)
        # Elemental sulfur oligomers S2..S8 for Willgerodt-Kindler / polysulfide chemistry
        for k in range(2, 9):
            units.append(Counter({"S": k}))
        _SMALL_HEAVY_UNITS = units
    return _SMALL_HEAVY_UNITS


def _heavy_delta_explained(delta: dict[str, int], extra_units: list[Counter[str]]) -> bool:
    if not delta:
        return True
    units = [*extra_units, *_small_heavy_units()]
    signed_units: list[dict[str, int]] = [{}]
    for u in units:
        if not u:
            continue
        pos = dict(u)
        neg = {k: -v for k, v in u.items()}
        if pos not in signed_units:
            signed_units.append(pos)
        if neg not in signed_units:
            signed_units.append(neg)
    for a in signed_units:
        if a == delta:
            return True
    for i, a in enumerate(signed_units):
        if not a:
            continue
        for b in signed_units[i:]:
            if not b:
                continue
            merged: dict[str, int] = {}
            for k in set(a) | set(b):
                val = a.get(k, 0) + b.get(k, 0)
                if val:
                    merged[k] = val
            if merged == delta:
                return True
    return False


def _normalize_step_entry(entry: str | dict[str, Any], idx: int) -> dict[str, Any]:
    if isinstance(entry, str):
        return {
            "step": idx,
            "smiles": entry.strip(),
            "type": "",
            "subtype": "",
            "label": "",
        }
    if not isinstance(entry, dict):
        raise ValueError(f"Step {idx} must be a SMILES string or an object.")
    smiles = (
        entry.get("intermediate_smiles")
        or entry.get("smiles")
        or entry.get("intermediate")
        or entry.get("product_smiles")
        or ""
    )
    return {
        "step": entry.get("step", idx),
        "smiles": str(smiles).strip(),
        "type": str(entry.get("type") or "").strip(),
        "subtype": str(entry.get("subtype") or "").strip(),
        "label": str(entry.get("label") or entry.get("description") or "").strip(),
    }


def mechanism_steps_check(
    steps: list[str | dict[str, Any]],
    reactants: str | list[str] | None = None,
    products: str | list[str] | None = None,
    conditions: str | None = None,
    reaction: str | None = None,
) -> dict[str, Any]:
    """Audit a proposed sequence of reaction intermediates step by step.

    Parameters
    ----------
    steps:
        Ordered list of intermediates, either plain SMILES strings or dicts
        carrying ``intermediate_smiles`` / ``smiles`` plus optional ``type``
        and ``subtype`` (compatible with oMeBench and custom mechanism JSON).
    reactants:
        Starting material SMILES (list or ``.``-separated string).
    products:
        Target product SMILES (list or ``.``-separated string).
    conditions:
        Free-text reagents/solvents over the arrow (e.g. ``"NaOH, H2O, rt"``).
    reaction:
        Optional ``reactants>reagents>products`` reaction SMILES shorthand when
        ``reactants`` / ``products`` are not passed separately.
    """
    if reaction and ">" in reaction:
        parts = reaction.split(">")
        if len(parts) == 3:
            if reactants is None and parts[0].strip():
                reactants = parts[0].strip()
            if not conditions and parts[1].strip():
                conditions = parts[1].strip()
            if products is None and parts[2].strip():
                products = parts[2].strip()

    if not isinstance(steps, list) or not steps:
        return {
            "status": "error",
            "issues": [{"severity": "blocking", "step": None, "code": "empty_steps", "message": "Provide a non-empty list of mechanism steps."}],
            "warnings": ["Provide a non-empty list of mechanism steps."],
        }

    try:
        normalized_steps = [_normalize_step_entry(item, i) for i, item in enumerate(steps, start=1)]
    except ValueError as exc:
        return {
            "status": "error",
            "issues": [{"severity": "blocking", "step": None, "code": "malformed_step", "message": str(exc)}],
            "warnings": [str(exc)],
        }

    r_tokens, r_mols, r_bad = _parse_species_list(reactants)
    p_tokens, p_mols, p_bad = _parse_species_list(products)
    cond_text = str(conditions or "").strip()

    issues: list[dict[str, Any]] = []
    for bad in r_bad:
        issues.append(
            {
                "severity": "blocking",
                "step": 0,
                "code": "unparseable_reactant",
                "message": f"Reactant SMILES could not be parsed by RDKit: {bad!r}",
            }
        )
    for bad in p_bad:
        issues.append(
            {
                "severity": "blocking",
                "step": None,
                "code": "unparseable_product",
                "message": f"Target product SMILES could not be parsed by RDKit: {bad!r}",
            }
        )

    inter_mols: list[Chem.Mol | None] = []
    for entry in normalized_steps:
        smi = entry["smiles"]
        mol = Chem.MolFromSmiles(smi) if smi else None
        inter_mols.append(mol)
        if mol is None:
            issues.append(
                {
                    "severity": "blocking",
                    "step": entry["step"],
                    "code": "unparseable_intermediate",
                    "message": f"Step {entry['step']} intermediate SMILES fails RDKit valence/syntax sanitization: {smi!r}",
                }
            )

    any_type_labels = any(entry["type"] or entry["subtype"] for entry in normalized_steps)
    declared_radical = any(
        any(kw in f"{entry['type']} {entry['subtype']} {cond_text}".lower() for kw in RADICAL_KEYWORDS)
        for entry in normalized_steps
    ) or any(bool(_odd_radical_atoms(m)) for m in r_mols + p_mols)

    allowed_elements: set[str] = set()
    for m in r_mols + p_mols:
        allowed_elements |= set(_heavy_counter(m))
    allowed_elements |= _reagent_elements_from_conditions(cond_text)

    reactant_heavy_units = [_heavy_counter(m) for m in r_mols]
    step_reports: list[dict[str, Any]] = []

    for idx, (entry, mol) in enumerate(zip(normalized_steps, inter_mols)):
        step_no = entry["step"]
        step_tag = f"{entry['type']} / {entry['subtype']}".strip(" /")
        if mol is None:
            step_reports.append({"step": step_no, "smiles": entry["smiles"], "valid_smiles": False})
            continue

        frags = _fragments(mol)
        main_frag = max(frags, key=lambda f: (f.GetNumHeavyAtoms(), f.GetNumAtoms()))
        full_formula = formula_counter(mol)
        heavy_counts = _heavy_counter(mol)
        charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        radicals = _odd_radical_atoms(mol)

        if radicals and not declared_radical:
            severity = "blocking" if any_type_labels else "review"
            issues.append(
                {
                    "severity": severity,
                    "step": step_no,
                    "code": "unexpected_radical_intermediate",
                    "message": (
                        f"Step {step_no} has odd-electron radical atom(s) ({radicals}) in a polar mechanism "
                        "without a declared radical/SET/homolytic step; check formal charges and valences."
                    ),
                }
            )

        if allowed_elements:
            foreign = sorted(set(heavy_counts) - allowed_elements)
            if foreign:
                issues.append(
                    {
                        "severity": "blocking",
                        "step": step_no,
                        "code": "element_from_nowhere",
                        "message": (
                            f"Step {step_no} introduces element(s) {foreign} not present in reactants, "
                            "products, or named reagents in conditions."
                        ),
                    }
                )

        # Compare against previous step (or reactants at step 1)
        prev_mol = inter_mols[idx - 1] if idx > 0 else None
        is_pt = any(kw in step_tag.lower() for kw in PROTON_TRANSFER_KEYWORDS)
        is_iso = any(kw in step_tag.lower() for kw in ISOMERIZATION_KEYWORDS)

        heavy_delta: dict[str, int] | None = None
        formula_delta: dict[str, int] | None = None
        charge_delta: int | None = None
        step_explained = True

        if idx == 0 and r_mols:
            step_explained = any(
                _heavy_delta_explained(_sub_counter(heavy_counts, _heavy_counter(r)), reactant_heavy_units)
                or _heavy_delta_explained(dict(heavy_counts), reactant_heavy_units)
                for r in r_mols
            )
            # Sum all reactants as reference delta for display
            total_r_formula: Counter[str] = Counter()
            total_r_heavy: Counter[str] = Counter()
            for r in r_mols:
                total_r_formula.update(formula_counter(r))
                total_r_heavy.update(_heavy_counter(r))
            formula_delta = _sub_counter(full_formula, total_r_formula)
            heavy_delta = _sub_counter(heavy_counts, total_r_heavy)
            charge_delta = charge - sum(sum(a.GetFormalCharge() for a in r.GetAtoms()) for r in r_mols)

            if is_pt and main_frag.GetNumHeavyAtoms() > 4:
                prev_skeletons = {_skeleton_hash(f) for r in r_mols for f in _fragments(r)}
                if _skeleton_hash(main_frag) not in prev_skeletons:
                    issues.append(
                        {
                            "severity": "blocking",
                            "step": step_no,
                            "code": "proton_transfer_changes_skeleton",
                            "message": (
                                f"Step {step_no} is labelled as a proton transfer ({step_tag}) "
                                "but alters the heavy-atom connectivity skeleton of the starting material."
                            ),
                        }
                    )
        elif idx > 0 and prev_mol is not None:
            prev_frags = _fragments(prev_mol)
            prev_main = max(prev_frags, key=lambda f: (f.GetNumHeavyAtoms(), f.GetNumAtoms()))
            prev_formula = formula_counter(prev_mol)
            prev_heavy = _heavy_counter(prev_mol)
            prev_charge = sum(a.GetFormalCharge() for a in prev_mol.GetAtoms())

            formula_delta = _sub_counter(full_formula, prev_formula)
            heavy_delta = _sub_counter(heavy_counts, prev_heavy)
            charge_delta = charge - prev_charge
            step_explained = _heavy_delta_explained(heavy_delta, reactant_heavy_units)

            if is_pt:
                prev_skeletons = {_skeleton_hash(f) for f in prev_frags}
                if _skeleton_hash(main_frag) not in prev_skeletons:
                    issues.append(
                        {
                            "severity": "blocking",
                            "step": step_no,
                            "code": "proton_transfer_changes_skeleton",
                            "message": (
                                f"Step {step_no} is labelled as a proton transfer ({step_tag}) "
                                "but changes heavy-atom connectivity relative to step "
                                f"{normalized_steps[idx - 1]['step']}."
                            ),
                        }
                    )
                else:
                    main_f_delta = _sub_counter(formula_counter(main_frag), formula_counter(prev_main))
                    non_h = {k: v for k, v in main_f_delta.items() if k != "H"}
                    if non_h and _skeleton_hash(main_frag) == _skeleton_hash(prev_main):
                        issues.append(
                            {
                                "severity": "review",
                                "step": step_no,
                                "code": "proton_transfer_heavy_atom_change",
                                "message": (
                                    f"Step {step_no} is labelled proton transfer but the main species formula "
                                    f"changes by {_delta_str(main_f_delta)}."
                                ),
                            }
                        )
            elif is_iso and heavy_delta:
                issues.append(
                    {
                        "severity": "review",
                        "step": step_no,
                        "code": "rearrangement_changes_heavy_formula",
                        "message": (
                            f"Step {step_no} is labelled as an intramolecular rearrangement/pericyclic step "
                            f"({step_tag}) but changes heavy-atom counts by {_delta_str(heavy_delta)}."
                        ),
                    }
                )

            if not step_explained:
                issues.append(
                    {
                        "severity": "review",
                        "step": step_no,
                        "code": "unexplained_step_delta",
                        "message": (
                            f"Step {step_no} heavy-atom change ({_delta_str(heavy_delta)}) from the previous "
                            "intermediate is not explained by adding a reactant or gaining/losing a common small molecule; "
                            "check for a dropped/extra carbon or missing intermediate."
                        ),
                    }
                )

        step_reports.append(
            {
                "step": step_no,
                "type": step_tag or None,
                "canonical_smiles": Chem.MolToSmiles(mol),
                "formula": hill_formula(full_formula),
                "formal_charge": charge,
                "formula_delta_from_previous": _delta_str(formula_delta) if formula_delta is not None else None,
                "heavy_delta_from_previous": _delta_str(heavy_delta) if heavy_delta is not None else None,
                "charge_delta_from_previous": charge_delta,
                "radical_atoms": len(radicals),
                "delta_explained": step_explained,
            }
        )

    final_vs_product: dict[str, Any] | None = None
    if p_mols and inter_mols and inter_mols[-1] is not None:
        last_mol = inter_mols[-1]
        last_frags = _fragments(last_mol)
        target_frags = [f for p in p_mols for f in _fragments(p)]
        organic_targets = [f for f in target_frags if any(a.GetSymbol() == "C" for a in f.GetAtoms())] or target_frags
        main_target = max(organic_targets, key=lambda f: (f.GetNumHeavyAtoms(), f.GetNumAtoms()))
        main_last = max(last_frags, key=lambda f: (f.GetNumHeavyAtoms(), f.GetNumAtoms()))

        last_nostereo = {_nostereo_smiles(f) for f in last_frags}
        target_nostereo = {_nostereo_smiles(f) for f in organic_targets}

        try:
            comp = compare_structures(Chem.MolToSmiles(main_last), Chem.MolToSmiles(main_target))
        except StructureError as exc:
            comp = {"relation": "error", "warnings": [str(exc)]}

        relation = comp.get("relation")
        final_vs_product = {
            "relation": relation,
            "last_intermediate_smiles": Chem.MolToSmiles(main_last),
            "target_product_smiles": Chem.MolToSmiles(main_target),
            "formula_delta": comp.get("formula_delta"),
            "differing_centres": comp.get("differing_centres", []),
        }

        if relation == "identical" or (last_nostereo & target_nostereo and relation == "same_constitution_stereo_undefined"):
            if relation == "same_constitution_stereo_undefined":
                issues.append(
                    {
                        "severity": "review",
                        "step": normalized_steps[-1]["step"],
                        "code": "final_product_stereo_undefined",
                        "message": "Last intermediate has the target connectivity, but stereochemistry differs in definition.",
                    }
                )
        elif relation == "constitutional_isomers":
            issues.append(
                {
                    "severity": "blocking",
                    "step": normalized_steps[-1]["step"],
                    "code": "incomplete_cascade_constitutional_isomer",
                    "message": (
                        "Last intermediate is a constitutional isomer of the given product (same formula, "
                        "different connectivity). Do NOT assume the target SMILES is a typo: look for a "
                        "subsequent ring-opening/reclosure (retro-aldol, retro-Claisen/Dieckmann), "
                        "1,2-shift, or rearrangement cascade that converts your intermediate into the target."
                    ),
                }
            )
        elif relation == "tautomers":
            issues.append(
                {
                    "severity": "blocking",
                    "step": normalized_steps[-1]["step"],
                    "code": "final_product_tautomer",
                    "message": (
                        "Last intermediate is a tautomer of the given product rather than the exact drawn "
                        "tautomer; add the concluding tautomerization / proton-transfer step(s)."
                    ),
                }
            )
        elif relation in {"enantiomers", "diastereomers"}:
            issues.append(
                {
                    "severity": "blocking",
                    "step": normalized_steps[-1]["step"],
                    "code": "final_product_stereoisomer_mismatch",
                    "message": (
                        f"Last intermediate and given product are {relation} "
                        f"(differing centres: {comp.get('differing_centres', [])})."
                    ),
                }
            )
        else:
            issues.append(
                {
                    "severity": "blocking",
                    "step": normalized_steps[-1]["step"],
                    "code": "final_product_mismatch",
                    "message": (
                        f"Last intermediate does not match the given product (relation={relation}, "
                        f"formula_delta={comp.get('formula_delta')})."
                    ),
                }
            )

    n_blocking = sum(1 for item in issues if item["severity"] == "blocking")
    n_review = sum(1 for item in issues if item["severity"] == "review")
    status = "blocked" if n_blocking else ("review" if n_review else "consistent")

    return {
        "status": status,
        "n_steps": len(normalized_steps),
        "n_blocking": n_blocking,
        "n_review": n_review,
        "issues": issues,
        "final_vs_product": final_vs_product,
        "steps": step_reports,
        "scope": (
            "Checks intermediate valence, radical/spin parity, proton-transfer skeleton invariance, "
            "element provenance, step-to-step atom/charge continuity, and closure to the target product. "
            "Passing ('consistent') confirms chemical bookkeeping, not kinetic barrier or regioselectivity."
        ),
    }
