from __future__ import annotations

import csv
import importlib.util
import json
import re
import tempfile
from pathlib import Path
from typing import Any

from .external import tool_statuses
from .rdkit_tools import RDKIT_AVAILABLE

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through public functions
    from rdkit import Chem
else:  # pragma: no cover
    Chem = None


ABBREVIATION_SMILES = {
    "H": None,
    "OH": "*O",
    "HO": "*O",
    "O": "*O",
    "SH": "*S",
    "HS": "*S",
    "NH2": "*N",
    "H2N": "*N",
    "NMe2": "*N(C)C",
    "NEt2": "*N(CC)CC",
    "OTs": "*OS(=O)(=O)c1ccc(C)cc1",
    "TsO": "*OS(=O)(=O)c1ccc(C)cc1",
    "SO2Tol": "*S(=O)(=O)c1ccc(C)cc1",
    "Ts": "*S(=O)(=O)c1ccc(C)cc1",
    "OMs": "*OS(=O)(=O)C",
    "MsO": "*OS(=O)(=O)C",
    "Ms": "*S(=O)(=O)C",
    "OTf": "*OS(=O)(=O)C(F)(F)F",
    "TfO": "*OS(=O)(=O)C(F)(F)F",
    "Tf": "*S(=O)(=O)C(F)(F)F",
    "NTf": "*NS(=O)(=O)C(F)(F)F",
    "SePh": "*[Se]c1ccccc1",
    "SPh": "*Sc1ccccc1",
    "SOPh": "*S(=O)c1ccccc1",
    "SO2Ph": "*S(=O)(=O)c1ccccc1",
    "MEM": "*COCCOC",
    "MEMO": "*OCOCCOC",
    "OMEM": "*OCOCCOC",
    "MOM": "*COC",
    "MOMO": "*OCOC",
    "OMOM": "*OCOC",
    "THP": "*C1CCCCO1",
    "OTHP": "*OC1CCCCO1",
    "Bn": "*Cc1ccccc1",
    "OBn": "*OCc1ccccc1",
    "BnO": "*OCc1ccccc1",
    "PMB": "*Cc1ccc(OC)cc1",
    "OPMB": "*OCc1ccc(OC)cc1",
    "PMBO": "*OCc1ccc(OC)cc1",
    "Tr": "*C(c1ccccc1)(c1ccccc1)c1ccccc1",
    "OTr": "*OC(c1ccccc1)(c1ccccc1)c1ccccc1",
    "Ac": "*C(=O)C",
    "OAc": "*OC(=O)C",
    "AcO": "*OC(=O)C",
    "NHAc": "*NC(=O)C",
    "Bz": "*C(=O)c1ccccc1",
    "OBz": "*OC(=O)c1ccccc1",
    "BzO": "*OC(=O)c1ccccc1",
    "Piv": "*C(=O)C(C)(C)C",
    "OPiv": "*OC(=O)C(C)(C)C",
    "PivO": "*OC(=O)C(C)(C)C",
    "Boc": "*C(=O)OC(C)(C)C",
    "NHBoc": "*NC(=O)OC(C)(C)C",
    "BocNH": "*NC(=O)OC(C)(C)C",
    "Cbz": "*C(=O)OCc1ccccc1",
    "NHCbz": "*NC(=O)OCc1ccccc1",
    "Fmoc": "*C(=O)OCC1c2ccccc2-c2ccccc12",
    "NHFmoc": "*NC(=O)OCC1c2ccccc2-c2ccccc12",
    "TMS": "*[Si](C)(C)C",
    "OTMS": "*O[Si](C)(C)C",
    "TMSO": "*O[Si](C)(C)C",
    "TES": "*[Si](CC)(CC)CC",
    "OTES": "*O[Si](CC)(CC)CC",
    "TBS": "*[Si](C)(C)C(C)(C)C",
    "TBDMS": "*[Si](C)(C)C(C)(C)C",
    "OTBS": "*O[Si](C)(C)C(C)(C)C",
    "TBSO": "*O[Si](C)(C)C(C)(C)C",
    "OTBDMS": "*O[Si](C)(C)C(C)(C)C",
    "TIPS": "*[Si](C(C)C)(C(C)C)C(C)C",
    "OTIPS": "*O[Si](C(C)C)(C(C)C)C(C)C",
    "TBDPS": "*[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "OTBDPS": "*O[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "CHO": "*C=O",
    "CO2H": "*C(=O)O",
    "COOH": "*C(=O)O",
    "CO2Me": "*C(=O)OC",
    "COOMe": "*C(=O)OC",
    "COOEt": "*C(=O)OCC",
    "CO2Et": "*C(=O)OCC",
    "COOC2H5": "*C(=O)OCC",
    "CO2tBu": "*C(=O)OC(C)(C)C",
    "CN": "*C#N",
    "NO2": "*[N+](=O)[O-]",
    "CF3": "*C(F)(F)F",
    "CCl3": "*C(Cl)(Cl)Cl",
    "Bpin": "*B1OC(C)(C)C(C)(C)O1",
    "Ph": "*c1ccccc1",
    "Me": "*C",
    "Et": "*CC",
    "nPr": "*CCC",
    "iPr": "*C(C)C",
    "nBu": "*CCCC",
    "iBu": "*CC(C)C",
    "sBu": "*C(C)CC",
    "tBu": "*C(C)(C)C",
    "Cy": "*C1CCCCC1",
    "Allyl": "*CC=C",
    "Vinyl": "*C=C",
}


def _mol_from_smiles(smiles: str, sanitize: bool = True) -> Any:
    if not RDKIT_AVAILABLE:
        return None
    try:
        return Chem.MolFromSmiles(smiles, sanitize=sanitize)
    except Exception:
        return None


def _canonical_smiles(smiles: str, *, isomeric: bool = False) -> str | None:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)
    except Exception:
        return None


def _has_low_confidence_label_artifact_atoms(smiles: str | None) -> bool:
    mol = _mol_from_smiles(smiles or "")
    if mol is None:
        return False
    # Elements routinely drawn in organic schemes; a low-confidence candidate
    # containing anything else (B, Si, metals, ...) probably misread an atom
    # label. Halogens must be here: aryl F/Cl/Br/I are everywhere in real
    # schemes, and penalizing them punishes correct reads.
    common_scheme_atoms = {0, 1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}
    return any(atom.GetAtomicNum() not in common_scheme_atoms for atom in mol.GetAtoms())


def _largest_fragment_smiles(smiles: str) -> str | None:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    frag = max(frags, key=lambda item: item.GetNumAtoms())
    return Chem.MolToSmiles(frag, canonical=True, isomericSmiles=False)


def structure_metrics(smiles: str | None) -> dict[str, Any]:
    if not smiles or not RDKIT_AVAILABLE:
        return {
            "sanitize_ok": False,
            "dummy_count": None,
            "fragment_count": None,
            "heavy_atom_count": None,
            "ring_count": None,
            "aromatic_ring_count": None,
            "stereo_warning": None,
            "canonical_smiles": None,
            "main_fragment_smiles": None,
        }
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return {
            "sanitize_ok": False,
            "dummy_count": None,
            "fragment_count": None,
            "heavy_atom_count": None,
            "ring_count": None,
            "aromatic_ring_count": None,
            "stereo_warning": None,
            "canonical_smiles": None,
            "main_fragment_smiles": None,
        }
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    aromatic_ring_count = 0
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    return {
        "sanitize_ok": True,
        "dummy_count": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0),
        "fragment_count": len(Chem.GetMolFrags(mol)),
        "heavy_atom_count": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1),
        "ring_count": mol.GetRingInfo().NumRings(),
        "aromatic_ring_count": aromatic_ring_count,
        "stereo_warning": any(label == "?" for _, label in chiral_centers),
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        "isomeric_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "main_fragment_smiles": _largest_fragment_smiles(smiles),
    }


def compare_smiles(candidate_smiles: str | None, gold_smiles: str | None) -> dict[str, Any]:
    candidate_metrics = structure_metrics(candidate_smiles)
    gold_metrics = structure_metrics(gold_smiles)
    return {
        "exact_graph_match": bool(
            candidate_metrics["canonical_smiles"]
            and candidate_metrics["canonical_smiles"] == gold_metrics["canonical_smiles"]
        ),
        "exact_stereo_match": bool(
            candidate_metrics.get("isomeric_smiles")
            and candidate_metrics.get("isomeric_smiles") == gold_metrics.get("isomeric_smiles")
        ),
        "main_fragment_match": bool(
            candidate_metrics["main_fragment_smiles"]
            and candidate_metrics["main_fragment_smiles"] == gold_metrics["main_fragment_smiles"]
        ),
        "candidate_metrics": candidate_metrics,
        "gold_metrics": gold_metrics,
    }


def candidate_score(record: dict[str, Any]) -> float:
    return _candidate_score(record)


def _candidate_score(record: dict[str, Any], *, max_heavy_atom_count: int | None = None) -> float:
    smiles = record.get("canonical_smiles") or record.get("isomeric_smiles") or record.get("smiles")
    metrics = structure_metrics(smiles)
    score = 0.0
    if metrics["sanitize_ok"]:
        score += 10.0
    else:
        score -= 10.0
    confidence = record.get("confidence")
    adapter_confidence = record.get("metadata", {}).get("adapter_confidence") if isinstance(record.get("metadata"), dict) else None
    selected_confidence: float | None = None
    for value in (adapter_confidence, confidence):
        try:
            if value is not None:
                selected_confidence = float(value)
                score += selected_confidence * 3.0
                break
        except (TypeError, ValueError):
            pass
    dummy_count = metrics["dummy_count"] or 0
    fragment_count = metrics["fragment_count"] or 1
    heavy_atom_count = metrics["heavy_atom_count"] or 0
    if metrics["sanitize_ok"]:
        score += min(heavy_atom_count, 60) * 0.08
    if max_heavy_atom_count and max_heavy_atom_count >= 8 and metrics["sanitize_ok"]:
        complexity_ratio = heavy_atom_count / max_heavy_atom_count if heavy_atom_count else 0.0
        if heavy_atom_count < 4:
            score -= 12.0
        elif complexity_ratio < 0.35:
            score -= 6.0
        elif complexity_ratio < 0.6:
            score -= 2.5
    if metrics["sanitize_ok"] and heavy_atom_count >= 8 and selected_confidence is not None and selected_confidence <= 0.6:
        score -= 2.0
        if _has_low_confidence_label_artifact_atoms(smiles):
            score -= 4.0
    score -= dummy_count * 1.1
    score -= max(0, fragment_count - 1) * 0.8
    if metrics["stereo_warning"]:
        score -= 0.35
    warnings = record.get("warnings") or []
    score -= min(len(warnings), 6) * 0.08
    metadata = record.get("metadata", {}) if isinstance(record.get("metadata"), dict) else {}
    # Consensus: independent engines agreeing on one structure is the strongest
    # signal OCSR has, worth more than any single engine's own confidence.
    # Agreement across preprocessing variants of the same image helps too, but
    # more weakly - variants are correlated evidence, not independent reads.
    try:
        tool_agreement = int(metadata.get("tool_agreement", 1) or 1)
    except (TypeError, ValueError):
        tool_agreement = 1
    try:
        variant_agreement = int(metadata.get("variant_agreement", 1) or 1)
    except (TypeError, ValueError):
        variant_agreement = 1
    if metrics["sanitize_ok"]:
        score += min(tool_agreement - 1, 3) * 2.2
        score += min(variant_agreement - 1, 4) * 0.7
    image_variant = metadata.get("image_variant")
    if image_variant and image_variant != "original":
        # Preprocessed variants often rescue faint skeletons, but can disturb
        # wedge/dash geometry. Let them win on graph quality, not on tiny score
        # ties against the original crop.
        score -= 0.2
    return score


def rank_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    max_heavy_atom_count = 0
    for candidate in candidates:
        smiles = candidate.get("canonical_smiles") or candidate.get("isomeric_smiles") or candidate.get("smiles")
        metrics = structure_metrics(smiles)
        if metrics["sanitize_ok"]:
            max_heavy_atom_count = max(max_heavy_atom_count, int(metrics["heavy_atom_count"] or 0))
    ranked = sorted(candidates, key=lambda item: _candidate_score(item, max_heavy_atom_count=max_heavy_atom_count), reverse=True)
    for idx, candidate in enumerate(ranked, start=1):
        candidate.setdefault("metadata", {})
        if isinstance(candidate["metadata"], dict):
            candidate["metadata"]["rank"] = idx
            candidate["metadata"]["ensemble_score"] = round(
                _candidate_score(candidate, max_heavy_atom_count=max_heavy_atom_count), 4
            )
    return ranked


def _definition_smiles(value: str) -> str | None:
    stripped = value.strip()
    if stripped in ABBREVIATION_SMILES:
        return ABBREVIATION_SMILES[stripped]
    for key, smiles in ABBREVIATION_SMILES.items():
        if key.upper() == stripped.upper():
            return smiles
    return stripped


def _clean_label(value: str) -> str:
    return value.strip().strip("[]").replace("^", "")


def _dummy_labels(atom: Any) -> list[str]:
    labels: list[str] = []
    if atom.HasProp("atomLabel"):
        labels.append(_clean_label(atom.GetProp("atomLabel")))
    if atom.GetAtomMapNum():
        labels.append(str(atom.GetAtomMapNum()))
    if atom.GetIsotope():
        labels.append(f"R{atom.GetIsotope()}")
        labels.append(str(atom.GetIsotope()))
    labels.append("*")
    return [label for label in labels if label]


def _definition_for_atom(atom: Any, definitions: dict[str, str], dummy_count: int) -> tuple[str | None, str | None]:
    key_map = {str(key): str(value) for key, value in definitions.items()}
    upper_map = {str(key).upper(): (str(key), str(value)) for key, value in definitions.items()}
    for label in _dummy_labels(atom):
        if label in key_map:
            return label, key_map[label]
        if label.upper() in upper_map:
            original, value = upper_map[label.upper()]
            return original, value
    if dummy_count == 1 and len(key_map) == 1:
        key, value = next(iter(key_map.items()))
        return key, value
    return None, None


def _fragment_with_attachment(definition_smiles: str) -> tuple[Any, int, list[str]]:
    warnings: list[str] = []
    frag = _mol_from_smiles(definition_smiles, sanitize=True)
    if frag is None:
        return None, 0, [f"Could not parse abbreviation replacement SMILES: {definition_smiles}"]
    dummy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() == 0]
    if not dummy_atoms:
        return frag, 0, warnings
    marker = dummy_atoms[0]
    neighbors = list(marker.GetNeighbors())
    if len(neighbors) != 1:
        return None, 0, [f"Replacement fragment has unsupported dummy attachment: {definition_smiles}"]
    attach_idx = neighbors[0].GetIdx()
    marker_idx = marker.GetIdx()
    rw = Chem.RWMol(frag)
    rw.RemoveAtom(marker_idx)
    if marker_idx < attach_idx:
        attach_idx -= 1
    return rw.GetMol(), attach_idx, warnings


def resolve_placeholders(smiles: str, definitions: dict[str, str]) -> dict[str, Any]:
    warnings: list[str] = []
    if not RDKIT_AVAILABLE:
        return {
            "status": "unavailable",
            "smiles": smiles,
            "warnings": ["RDKit is unavailable; dummy replacement could not be validated."],
            "unresolved_labels": list(definitions),
        }
    mol = _mol_from_smiles(smiles, sanitize=False)
    if mol is None:
        return {"status": "error", "smiles": smiles, "warnings": [f"Invalid placeholder SMILES: {smiles}"]}
    rw = Chem.RWMol(mol)
    used_definitions: dict[str, str] = {}
    unresolved: list[str] = []
    for _ in range(128):
        current = rw.GetMol()
        dummy_atoms = [atom for atom in current.GetAtoms() if atom.GetAtomicNum() == 0]
        if not dummy_atoms:
            break
        target = None
        target_key = None
        target_value = None
        for atom in dummy_atoms:
            key, value = _definition_for_atom(atom, definitions, len(dummy_atoms))
            if value is not None:
                target = atom
                target_key = key
                target_value = value
                break
        if target is None or target_value is None:
            for atom in dummy_atoms:
                unresolved.extend(_dummy_labels(atom)[:1])
            break
        idx = target.GetIdx()
        neighbors = list(target.GetNeighbors())
        if len(neighbors) > 1:
            unresolved.append(target_key or _dummy_labels(target)[0])
            warnings.append(f"Dummy {target_key or '*'} has more than one neighbor; leaving it unresolved.")
            break
        if len(neighbors) == 0:
            rw.RemoveAtom(idx)
            warnings.append(f"Removed disconnected dummy {target_key or '*'} during placeholder resolution.")
            continue
        parent_idx = neighbors[0].GetIdx()
        parent_atom = current.GetAtomWithIdx(parent_idx)
        parent_tag = parent_atom.GetChiralTag()
        old_nbrs = [n.GetIdx() for n in parent_atom.GetNeighbors()]
        bond = current.GetBondBetweenAtoms(parent_idx, idx)
        bond_type = bond.GetBondType() if bond is not None else Chem.BondType.SINGLE
        replacement_smiles = _definition_smiles(target_value)
        used_definitions[str(target_key)] = target_value
        if replacement_smiles is None:
            # Convert dummy to explicit H in place so tetrahedral bond-order parity is preserved, then RemoveHs later.
            h_atom = rw.GetAtomWithIdx(idx)
            h_atom.SetAtomicNum(1)
            h_atom.SetIsotope(0)
            h_atom.SetAtomMapNum(0)
            if h_atom.HasProp("atomLabel"):
                h_atom.ClearProp("atomLabel")
            continue
        fragment, attach_idx, frag_warnings = _fragment_with_attachment(replacement_smiles)
        warnings.extend(frag_warnings)
        if fragment is None:
            unresolved.append(target_key or _dummy_labels(target)[0])
            break
        rw.RemoveAtom(idx)
        shifted_parent_idx = parent_idx - 1 if idx < parent_idx else parent_idx
        base = rw.GetMol()
        offset = base.GetNumAtoms()
        new_attach_idx = offset + attach_idx
        combo = Chem.CombineMols(base, fragment)
        rw = Chem.RWMol(combo)
        rw.AddBond(shifted_parent_idx, new_attach_idx, bond_type)
        if parent_tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            mapped_old = [
                new_attach_idx if u == idx else (u - 1 if u > idx else u)
                for u in old_nbrs
            ]
            new_nbrs = [n.GetIdx() for n in rw.GetAtomWithIdx(shifted_parent_idx).GetNeighbors()]
            if len(new_nbrs) == len(mapped_old) and set(new_nbrs) == set(mapped_old):
                from .structure_compare import _permutation_parity

                rw.GetAtomWithIdx(shifted_parent_idx).SetChiralTag(parent_tag)
                if _permutation_parity(new_nbrs, mapped_old) == 1:
                    rw.GetAtomWithIdx(shifted_parent_idx).InvertChirality()
    resolved_mol = rw.GetMol()
    try:
        Chem.SanitizeMol(resolved_mol)
        resolved_mol = Chem.RemoveHs(resolved_mol)
        resolved_smiles = Chem.MolToSmiles(resolved_mol, canonical=True, isomericSmiles=True)
    except Exception as exc:
        resolved_smiles = Chem.MolToSmiles(resolved_mol, canonical=True, isomericSmiles=True)
        warnings.append(f"Resolved structure did not fully sanitize: {exc}")
    unresolved_labels = sorted(set(label for label in unresolved if label))
    if any(atom.GetAtomicNum() == 0 for atom in resolved_mol.GetAtoms()):
        unresolved_labels.extend(
            label
            for atom in resolved_mol.GetAtoms()
            if atom.GetAtomicNum() == 0
            for label in _dummy_labels(atom)[:1]
            if label not in unresolved_labels
        )
    return {
        "status": "resolved" if not unresolved_labels else "partial",
        "input_smiles": smiles,
        "smiles": resolved_smiles,
        "definitions": used_definitions,
        "unresolved_labels": unresolved_labels,
        "warnings": warnings,
        "metrics": structure_metrics(resolved_smiles),
    }


def _safe_extract_py_literals(source: Path) -> dict[str, Any]:
    import ast

    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    extracted: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {"A_TO_I", "A_TO_I_NUMBERS"}:
                    extracted[target.id] = ast.literal_eval(node.value)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
            if node.target.id in {"A_TO_I", "A_TO_I_NUMBERS"}:
                extracted[node.target.id] = ast.literal_eval(node.value)
    return extracted


def load_gold_smiles(path: str | Path) -> list[dict[str, str]]:
    source = Path(path).expanduser().resolve()
    if source.suffix == ".py":
        extracted = _safe_extract_py_literals(source)
        entries = extracted.get("A_TO_I")
        number_map = extracted.get("A_TO_I_NUMBERS", {})
        if entries is None:
            raise ValueError(f"{source} does not define a literal A_TO_I list/tuple.")
        rows: list[dict[str, str]] = []
        for item in entries:
            if len(item) == 3:
                label, compound, smiles = item
            elif len(item) == 2:
                label, smiles = item
                compound = number_map.get(label, "")
            else:
                raise ValueError(f"Unsupported A_TO_I entry: {item!r}")
            rows.append({"label": str(label), "compound": str(compound), "smiles": str(smiles)})
        return rows
    delimiter = "\t" if source.suffix.lower() in {".tsv", ".tab"} else ","
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = []
        for row in reader:
            label = row.get("label") or row.get("letter") or row.get("id") or ""
            compound = row.get("compound") or row.get("number") or row.get("source") or ""
            smiles = (
                row.get("smiles")
                or row.get("gold_smiles")
                or row.get("canonical_smiles")
                or row.get("current_canonical_smiles")
                or ""
            )
            if smiles:
                rows.append({"label": label, "compound": compound, "smiles": smiles})
        return rows


def load_legend_map(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    if source.suffix.lower() == ".json":
        data = json.loads(source.read_text(encoding="utf-8"))
        result: dict[str, dict[str, str]] = {}
        for compound, definitions in data.items():
            if isinstance(definitions, dict):
                result[str(compound)] = {str(key): str(value) for key, value in definitions.items()}
        return result
    delimiter = "\t" if source.suffix.lower() in {".tsv", ".tab"} else ","
    result: dict[str, dict[str, str]] = {}
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            compound = row.get("compound") or row.get("number") or row.get("source") or ""
            placeholder = row.get("placeholder") or row.get("label") or row.get("symbol") or ""
            value = row.get("value") or row.get("definition") or row.get("group") or ""
            if compound and placeholder and value:
                result.setdefault(str(compound), {})[str(placeholder)] = str(value)
    return result


def _find_crop(image_dir: Path, label: str, compound: str) -> Path | None:
    patterns = [
        f"{label}_source_{compound}.png",
        f"{label}_*.png",
        f"*_{compound}.png",
        f"{compound}.png",
        f"{label}.png",
    ]
    for pattern in patterns:
        matches = sorted(image_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def _candidate_smiles(candidate: dict[str, Any]) -> str | None:
    return (
        candidate.get("canonical_smiles")
        or candidate.get("isomeric_smiles")
        or candidate.get("smiles")
        or candidate.get("reaction_smiles")
    )


def benchmark_ocsr(gold_smiles: str | Path, image_dir: str | Path) -> dict[str, Any]:
    from .ocsr import parse_image

    gold_rows = load_gold_smiles(gold_smiles)
    crop_dir = Path(image_dir).expanduser().resolve()
    rows: list[dict[str, Any]] = []
    tool_summary: dict[str, dict[str, int]] = {}
    exact_matches = 0
    main_fragment_matches = 0
    no_candidates = 0
    for gold in gold_rows:
        image = _find_crop(crop_dir, gold["label"], gold["compound"])
        if image is None:
            rows.append({**gold, "status": "missing_crop", "candidates": [], "warnings": ["No matching crop image found."]})
            no_candidates += 1
            continue
        parsed = parse_image(str(image), kind="molecule")
        candidates = []
        for candidate in parsed.get("candidates", []):
            smiles = _candidate_smiles(candidate)
            comparison = compare_smiles(smiles, gold["smiles"])
            metrics = comparison["candidate_metrics"]
            tool = candidate.get("metadata", {}).get("ocsr_tool", candidate.get("source", "unknown"))
            candidate_row = {
                "tool": tool,
                "smiles": smiles,
                "confidence": candidate.get("confidence"),
                "exact_graph_match": comparison["exact_graph_match"],
                "main_fragment_match": comparison["main_fragment_match"],
                "dummy_count": metrics["dummy_count"],
                "fragment_count": metrics["fragment_count"],
                "stereo_warning": metrics["stereo_warning"],
                "warnings": candidate.get("warnings", []),
            }
            tool_summary.setdefault(tool, {"candidates": 0, "exact_graph_matches": 0, "main_fragment_matches": 0})
            tool_summary[tool]["candidates"] += 1
            if comparison["exact_graph_match"]:
                tool_summary[tool]["exact_graph_matches"] += 1
            if comparison["main_fragment_match"]:
                tool_summary[tool]["main_fragment_matches"] += 1
            candidates.append(candidate_row)
        best = candidates[0] if candidates else None
        if best and best["exact_graph_match"]:
            exact_matches += 1
        if best and best["main_fragment_match"]:
            main_fragment_matches += 1
        if not candidates:
            no_candidates += 1
        rows.append(
            {
                **gold,
                "image": str(image),
                "status": "ok" if candidates else "no_candidates",
                "best": best,
                "candidates": candidates,
                "warnings": parsed.get("warnings", []),
            }
        )
    return {
        "status": "ok",
        "gold_smiles": str(Path(gold_smiles).expanduser().resolve()),
        "image_dir": str(crop_dir),
        "summary": {
            "gold_count": len(gold_rows),
            "evaluated_count": len(rows),
            "exact_graph_matches": exact_matches,
            "main_fragment_matches": main_fragment_matches,
            "no_candidates": no_candidates,
            "tools": tool_summary,
        },
        "rows": rows,
        "tool_status": tool_statuses(),
    }


def benchmark_to_tsv(payload: dict[str, Any]) -> str:
    columns = [
        "label",
        "compound",
        "status",
        "tool",
        "smiles",
        "exact_graph_match",
        "main_fragment_match",
        "dummy_count",
        "fragment_count",
        "stereo_warning",
        "warnings",
    ]
    lines = ["\t".join(columns)]
    for row in payload.get("rows", []):
        candidates = row.get("candidates") or [{}]
        for candidate in candidates:
            values = {
                "label": row.get("label", ""),
                "compound": row.get("compound", ""),
                "status": row.get("status", ""),
                "tool": candidate.get("tool", ""),
                "smiles": candidate.get("smiles", ""),
                "exact_graph_match": candidate.get("exact_graph_match", ""),
                "main_fragment_match": candidate.get("main_fragment_match", ""),
                "dummy_count": candidate.get("dummy_count", ""),
                "fragment_count": candidate.get("fragment_count", ""),
                "stereo_warning": candidate.get("stereo_warning", ""),
                "warnings": "; ".join(str(item) for item in candidate.get("warnings", row.get("warnings", []))),
            }
            lines.append("\t".join(str(values[column]) for column in columns))
    return "\n".join(lines) + "\n"


def _natural_sort_key(path: Path) -> list[int | str]:
    parts = re.split(r"(\d+)", path.stem.lower())
    return [int(part) if part.isdigit() else part for part in parts]


def _infer_crop_identity(path: Path) -> tuple[str, str]:
    match = re.match(r"(?P<label>[A-Za-z]{1,2})_source_(?P<compound>\d+)", path.stem)
    if match:
        return match.group("label").upper(), match.group("compound")
    match = re.match(r"(?:crop|cmpd|compound|mol|item|opt)[_-]?(?P<compound>\d+)$", path.stem, re.IGNORECASE)
    if match:
        return path.stem, match.group("compound")
    match = re.match(r"(?P<label>[A-Za-z])[_-](?P<compound>.*)$", path.stem)
    if match:
        return match.group("label").upper(), match.group("compound") if match.group("compound").isdigit() else ""
    return path.stem, path.stem if path.stem.isdigit() else ""


def parse_scheme(
    image: str | Path,
    crops: str | Path | list[dict[str, Any]],
    gold_map: str | Path | None = None,
) -> dict[str, Any]:
    from .image_tools import crop_image
    from .ocsr import parse_image

    try:
        image_path = Path(str(image)).expanduser().resolve()
    except OSError as exc:
        return {"status": "error", "warnings": [f"Invalid scheme image path: {exc}"], "compounds": []}

    temp_crop_dir: tempfile.TemporaryDirectory[str] | None = None
    crop_dir: Path | None = None
    if isinstance(crops, list) or (isinstance(crops, str) and crops.lstrip().startswith(("[", "{"))):
        temp_crop_dir = tempfile.TemporaryDirectory(prefix="codex-chem-scheme-crops-")
        cropped = crop_image(image_path, boxes=crops, output_dir=temp_crop_dir.name, scale=2.0)
        if cropped.get("status") != "ok":
            temp_crop_dir.cleanup()
            return {"status": "error", "image": str(image_path), "warnings": cropped.get("warnings", []), "compounds": []}
        crop_dir = Path(temp_crop_dir.name)
    else:
        try:
            crop_dir = Path(crops).expanduser().resolve()
        except OSError as exc:
            return {
                "status": "error",
                "image": str(image_path),
                "warnings": [
                    f"Invalid crops directory path ({exc}). Pass either a directory of crop images "
                    "or a JSON list of {'id': 'A', 'box': [x1, y1, x2, y2]} crop specifications."
                ],
                "compounds": [],
            }
        if not crop_dir.is_dir():
            return {
                "status": "error",
                "image": str(image_path),
                "warnings": [f"Crops directory does not exist: {crop_dir}"],
                "compounds": [],
            }

    try:
        legend_definitions = load_legend_map(gold_map)
        compounds: list[dict[str, Any]] = []
        allowed_exts = {".png", ".jpg", ".jpeg", ".webp"}
        crop_files = sorted(
            (p for p in crop_dir.iterdir() if p.is_file() and p.suffix.lower() in allowed_exts and p.stem != "grid_overlay"),
            key=_natural_sort_key,
        )
        for crop in crop_files:
            label, compound = _infer_crop_identity(crop)
            parsed = parse_image(str(crop), kind="molecule")
            raw_candidates = parsed.get("candidates", [])
            definitions = legend_definitions.get(compound, {}) or legend_definitions.get(label, {})
            resolved_candidates = []
            for candidate in raw_candidates:
                smiles = _candidate_smiles(candidate)
                if smiles and definitions:
                    resolved_candidates.append(resolve_placeholders(smiles, definitions))
            top_smiles = None
            if resolved_candidates:
                top_smiles = resolved_candidates[0].get("smiles")
            elif raw_candidates:
                top_smiles = _candidate_smiles(raw_candidates[0])
            compounds.append(
                {
                    "label": label,
                    "compound": compound,
                    "crop": str(crop),
                    "raw_candidates": raw_candidates,
                    "legend_definitions": definitions,
                    "resolved_candidates": resolved_candidates,
                    "selected_smiles_for_consistency": top_smiles,
                    "warnings": parsed.get("warnings", []),
                }
            )
        consistency_input = [
            {"label": item.get("label") or item.get("compound") or str(idx), "smiles": item.get("selected_smiles_for_consistency")}
            for idx, item in enumerate(compounds, start=1)
            if item.get("selected_smiles_for_consistency")
        ]
        return {
            "status": "ok",
            "image": str(image_path),
            "crops": str(crop_dir) if temp_crop_dir is None else "inline_boxes",
            "legend_definitions_source": str(Path(gold_map).expanduser().resolve()) if gold_map else "not_provided",
            "legend_definitions": legend_definitions,
            "compounds": compounds,
            "route_consistency_warnings": route_consistency_warnings(consistency_input),
            "tool_status": tool_statuses(),
        }
    finally:
        if temp_crop_dir is not None:
            temp_crop_dir.cleanup()


def _benzene_relations(smiles: str | None) -> list[str]:
    mol = _mol_from_smiles(smiles or "")
    if mol is None:
        return []
    relations: list[str] = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        substituted_positions = []
        ring_set = set(ring)
        for pos, idx in enumerate(ring):
            atom = mol.GetAtomWithIdx(idx)
            has_external = any(neighbor.GetIdx() not in ring_set and neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
            if has_external:
                substituted_positions.append(pos)
        if len(substituted_positions) == 2:
            distance = abs(substituted_positions[0] - substituted_positions[1])
            distance = min(distance, 6 - distance)
            relations.append({1: "ortho", 2: "meta", 3: "para"}[distance])
    return sorted(relations)


def route_consistency_warnings(rows: list[dict[str, Any]]) -> list[str]:
    warnings: list[str] = []
    previous: dict[str, Any] | None = None
    for row in rows:
        smiles = row.get("smiles")
        metrics = structure_metrics(smiles)
        relations = _benzene_relations(smiles)
        current = {**row, "metrics": metrics, "benzene_relations": relations}
        if previous is not None:
            prev_label = previous.get("label", "?")
            label = current.get("label", "?")
            prev_relations = previous.get("benzene_relations") or []
            if prev_relations and relations and prev_relations != relations:
                warnings.append(
                    f"{prev_label}->{label} aromatic substitution relation changed from "
                    f"{','.join(prev_relations)} to {','.join(relations)}; likely OCSR error unless the route rearranges that arene."
                )
            prev_metrics = previous.get("metrics", {})
            if metrics["sanitize_ok"] and prev_metrics.get("sanitize_ok"):
                if abs((metrics["ring_count"] or 0) - (prev_metrics.get("ring_count") or 0)) > 1:
                    warnings.append(f"{prev_label}->{label} ring count changed sharply; verify core scaffold recognition.")
                prev_heavy = prev_metrics.get("heavy_atom_count") or 0
                heavy = metrics["heavy_atom_count"] or 0
                if prev_heavy and heavy and (heavy < prev_heavy * 0.7 or heavy > prev_heavy * 1.45):
                    warnings.append(f"{prev_label}->{label} heavy atom count changed sharply; verify side-chain anchors and fragments.")
        previous = current
    return warnings
