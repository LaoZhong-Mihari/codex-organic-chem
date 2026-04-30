from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .models import CalculationRecord, MoleculeRecord

try:  # pragma: no cover - import availability is environment dependent
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D

    RDLogger.DisableLog("rdApp.warning")
    RDKIT_AVAILABLE = True
    RDKIT_IMPORT_ERROR: Exception | None = None
except Exception as exc:  # pragma: no cover
    Chem = None
    AllChem = None
    Crippen = None
    Descriptors = None
    Draw = None
    rdBase = None
    rdMolDescriptors = None
    rdMolDraw2D = None
    RDKIT_AVAILABLE = False
    RDKIT_IMPORT_ERROR = exc


METAL_ATOMIC_NUMBERS = {
    3,
    4,
    11,
    12,
    13,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    55,
    56,
    57,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
}

FUNCTIONAL_GROUP_SMARTS = {
    "alcohol": "[OX2H][CX4]",
    "phenol": "[OX2H]c",
    "amine": "[NX3;H2,H1,H0;!$(NC=O)]",
    "amide": "[NX3][CX3](=[OX1])",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[CX3](=O)[OX2H0][#6]",
    "aldehyde": "[CX3H1](=O)[#6,H]",
    "ketone": "[#6][CX3](=O)[#6]",
    "alkene": "C=C",
    "alkyne": "C#C",
    "aryl_halide": "c[F,Cl,Br,I]",
    "alkyl_halide": "[CX4][F,Cl,Br,I]",
    "boronic_acid": "B(O)O",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
}


@dataclass
class ParsedMol:
    mol: Any
    warnings: list[str]


def rdkit_version() -> str | None:
    if not RDKIT_AVAILABLE:
        return None
    return rdBase.rdkitVersion


def unavailable_record(source: str, message: str | None = None) -> MoleculeRecord:
    warnings = ["RDKit is unavailable; structure could not be validated."]
    if message:
        warnings.append(message)
    elif RDKIT_IMPORT_ERROR:
        warnings.append(str(RDKIT_IMPORT_ERROR))
    return MoleculeRecord(source=source, confidence=0.0, warnings=warnings)


def parse_molecule(smiles: str | None = None, molfile: str | None = None) -> ParsedMol:
    if not RDKIT_AVAILABLE:
        return ParsedMol(None, ["RDKit is unavailable."])
    warnings: list[str] = []
    if smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return ParsedMol(None, [f"Invalid SMILES: {smiles}"])
    elif molfile:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False)
        if mol is None:
            return ParsedMol(None, ["Invalid MolBlock input."])
    else:
        return ParsedMol(None, ["Either smiles or molfile is required."])
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        warnings.append(f"RDKit sanitization warning: {exc}")
    return ParsedMol(mol, warnings)


def molecule_to_svg(mol: Any, width: int = 420, height: int = 300) -> str:
    if not RDKIT_AVAILABLE or mol is None:
        return ""
    draw_mol = Chem.Mol(mol)
    AllChem.Compute2DCoords(draw_mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, draw_mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def reaction_to_svg(reaction_smiles: str, width: int = 900, height: int = 260) -> tuple[str | None, list[str]]:
    if not RDKIT_AVAILABLE:
        return None, ["RDKit is unavailable; reaction drawing could not be generated."]
    warnings: list[str] = []
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if rxn is None:
            return None, ["RDKit could not parse reaction SMILES."]
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText(), warnings
    except Exception as exc:
        return None, [f"Reaction drawing failed: {exc}"]


def molecule_record_from_mol(mol: Any, source: str, extra_warnings: list[str] | None = None) -> MoleculeRecord:
    if not RDKIT_AVAILABLE or mol is None:
        return unavailable_record(source)
    warnings = list(extra_warnings or [])
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    molblock = Chem.MolToMolBlock(mol)
    svg = molecule_to_svg(mol)
    inchi_key = None
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as exc:  # pragma: no cover - optional RDKit InChI build
        warnings.append(f"InChIKey generation unavailable: {exc}")
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if any(label == "?" for _, label in chiral_centers):
        warnings.append("One or more stereocenters are unassigned; verify stereochemistry manually.")
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        warnings.append("Wildcard/dummy atoms detected; polymer/query structures require manual review.")
    if any(atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms()):
        warnings.append("Metal-containing or organometallic structures may be outside organic-rule heuristics.")
    if mol.GetNumHeavyAtoms() > 80:
        warnings.append("Large molecule detected; OCSR and mechanism heuristics need manual review.")
    confidence = 0.95 if not warnings else max(0.5, 0.9 - 0.08 * len(warnings))
    metadata = {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "rdkit_version": rdkit_version(),
        "functional_groups": detect_functional_groups(mol),
    }
    return MoleculeRecord(
        source=source,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        inchi_key=inchi_key,
        molblock=molblock,
        svg=svg,
        confidence=confidence,
        warnings=warnings,
        metadata=metadata,
    )


def normalize_structure(smiles: str | None = None, molfile: str | None = None, source: str = "user") -> MoleculeRecord:
    parsed = parse_molecule(smiles=smiles, molfile=molfile)
    if parsed.mol is None:
        return MoleculeRecord(source=source, confidence=0.0, warnings=parsed.warnings)
    return molecule_record_from_mol(parsed.mol, source=source, extra_warnings=parsed.warnings)


def detect_functional_groups(mol: Any) -> list[str]:
    if not RDKIT_AVAILABLE or mol is None:
        return []
    groups: list[str] = []
    for name, smarts in FUNCTIONAL_GROUP_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            groups.append(name)
    return groups


def descriptors_record(smiles: str) -> CalculationRecord:
    parsed = parse_molecule(smiles=smiles)
    if parsed.mol is None:
        return CalculationRecord(
            method="rdkit_descriptors",
            tool_version=rdkit_version(),
            parameters={"smiles": smiles},
            status="error",
            warnings=parsed.warnings,
        )
    mol = parsed.mol
    results = {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "exact_molecular_weight": Descriptors.ExactMolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "h_bond_donors": rdMolDescriptors.CalcNumHBD(mol),
        "h_bond_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(mol),
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "functional_groups": detect_functional_groups(mol),
    }
    return CalculationRecord(
        method="rdkit_descriptors",
        tool_version=rdkit_version(),
        parameters={"smiles": smiles},
        results=results,
        warnings=parsed.warnings,
    )


def conformer_record(smiles: str, num_confs: int = 8, max_iters: int = 200) -> CalculationRecord:
    parsed = parse_molecule(smiles=smiles)
    if parsed.mol is None:
        return CalculationRecord(
            method="rdkit_conformers",
            tool_version=rdkit_version(),
            parameters={"smiles": smiles, "num_confs": num_confs, "max_iters": max_iters},
            status="error",
            warnings=parsed.warnings,
        )
    mol = Chem.AddHs(parsed.mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 12648430
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(num_confs), params=params))
    warnings = list(parsed.warnings)
    if not conf_ids:
        return CalculationRecord(
            method="rdkit_conformers",
            tool_version=rdkit_version(),
            parameters={"smiles": smiles, "num_confs": num_confs, "max_iters": max_iters},
            status="error",
            warnings=warnings + ["RDKit conformer embedding produced no conformers."],
        )
    energies = []
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    for conf_id in conf_ids:
        energy = None
        if mmff_props is not None:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
                if ff is not None:
                    ff.Minimize(maxIts=int(max_iters))
                    energy = ff.CalcEnergy()
            except Exception as exc:
                warnings.append(f"MMFF optimization failed for conformer {conf_id}: {exc}")
        if energy is None:
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                ff.Minimize(maxIts=int(max_iters))
                energy = ff.CalcEnergy()
            except Exception as exc:
                warnings.append(f"UFF optimization failed for conformer {conf_id}: {exc}")
        energies.append({"conformer_id": int(conf_id), "energy": energy})
    best = min((e for e in energies if e["energy"] is not None), key=lambda x: x["energy"], default=None)
    return CalculationRecord(
        method="rdkit_conformers",
        tool_version=rdkit_version(),
        parameters={"smiles": smiles, "num_confs": num_confs, "max_iters": max_iters},
        results={
            "conformer_count": len(conf_ids),
            "energies": energies,
            "lowest_energy": best,
            "molblock_3d": Chem.MolToMolBlock(mol, confId=best["conformer_id"]) if best else None,
        },
        warnings=warnings,
    )


def charges_record(smiles: str) -> CalculationRecord:
    parsed = parse_molecule(smiles=smiles)
    if parsed.mol is None:
        return CalculationRecord(
            method="rdkit_gasteiger_charges",
            tool_version=rdkit_version(),
            parameters={"smiles": smiles},
            status="error",
            warnings=parsed.warnings,
        )
    mol = Chem.AddHs(parsed.mol)
    warnings = list(parsed.warnings)
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception as exc:
        return CalculationRecord(
            method="rdkit_gasteiger_charges",
            tool_version=rdkit_version(),
            parameters={"smiles": smiles},
            status="error",
            warnings=warnings + [f"Gasteiger charge calculation failed: {exc}"],
        )
    atom_charges = []
    for atom in mol.GetAtoms():
        raw = atom.GetProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else "nan"
        try:
            charge = float(raw)
        except ValueError:
            charge = None
        atom_charges.append(
            {
                "atom_index": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
                "gasteiger_charge": charge,
            }
        )
    return CalculationRecord(
        method="rdkit_gasteiger_charges",
        tool_version=rdkit_version(),
        parameters={"smiles": smiles},
        results={"atom_charges": atom_charges},
        warnings=warnings,
    )

