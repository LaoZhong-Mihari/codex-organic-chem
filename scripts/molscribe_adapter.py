#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any


def _checkpoint_path(explicit: str | None) -> str:
    if explicit:
        return explicit
    env_path = os.environ.get("CODEX_CHEM_MOLSCRIBE_CHECKPOINT")
    if env_path:
        return env_path
    os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
    try:
        from huggingface_hub import hf_hub_download
    except Exception as exc:
        raise RuntimeError(
            "MolScribe checkpoint was not provided and huggingface_hub is unavailable. "
            "Install huggingface_hub or set CODEX_CHEM_MOLSCRIBE_CHECKPOINT."
        ) from exc
    return hf_hub_download("yujieq/MolScribe", "swin_base_char_aux_1m.pth")


def _bond_type(value: str):
    from rdkit import Chem

    normalized = str(value).lower()
    if "triple" in normalized:
        return Chem.BondType.TRIPLE
    if "double" in normalized:
        return Chem.BondType.DOUBLE
    if "aromatic" in normalized:
        return Chem.BondType.AROMATIC
    return Chem.BondType.SINGLE


def _normalized_label(label: str) -> str:
    return label.strip().strip("[]").replace(" ", "").replace("-", "")


def _add_atom(mol, symbol: str, *, charge: int = 0, aromatic: bool = False) -> int:
    from rdkit import Chem

    atom = Chem.Atom(symbol)
    atom.SetFormalCharge(charge)
    atom.SetIsAromatic(aromatic)
    return mol.AddAtom(atom)


def _add_methanesulfonyl(mol, root_symbol: str) -> int:
    from rdkit import Chem

    root = _add_atom(mol, root_symbol)
    sulfur = _add_atom(mol, "S")
    oxygen_a = _add_atom(mol, "O")
    oxygen_b = _add_atom(mol, "O")
    methyl = _add_atom(mol, "C")
    mol.AddBond(root, sulfur, Chem.BondType.SINGLE)
    mol.AddBond(sulfur, oxygen_a, Chem.BondType.DOUBLE)
    mol.AddBond(sulfur, oxygen_b, Chem.BondType.DOUBLE)
    mol.AddBond(sulfur, methyl, Chem.BondType.SINGLE)
    return root


def _add_tosyl(mol, root_symbol: str) -> int:
    from rdkit import Chem

    root = _add_atom(mol, root_symbol)
    sulfur = _add_atom(mol, "S")
    oxygen_a = _add_atom(mol, "O")
    oxygen_b = _add_atom(mol, "O")
    aryl = [_add_atom(mol, "C", aromatic=True) for _ in range(6)]
    methyl = _add_atom(mol, "C")
    mol.AddBond(root, sulfur, Chem.BondType.SINGLE)
    mol.AddBond(sulfur, oxygen_a, Chem.BondType.DOUBLE)
    mol.AddBond(sulfur, oxygen_b, Chem.BondType.DOUBLE)
    mol.AddBond(sulfur, aryl[0], Chem.BondType.SINGLE)
    for idx in range(6):
        mol.AddBond(aryl[idx], aryl[(idx + 1) % 6], Chem.BondType.AROMATIC)
    mol.AddBond(aryl[3], methyl, Chem.BondType.SINGLE)
    return root


def _add_abbreviation(mol, label: str) -> int | None:
    from rdkit import Chem

    normalized = _normalized_label(label)
    if normalized in {"C@", "C@@", "CH@", "CH@@", "C@H", "C@@H"}:
        atom = Chem.Atom("C")
        atom.SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CW if "@@" not in normalized else Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        )
        return mol.AddAtom(atom)
    if normalized in {"Me", "CH3"}:
        return _add_atom(mol, "C")
    if normalized in {"tBu", "tertBu", "tertButyl", "tButyl"}:
        central = _add_atom(mol, "C")
        for _ in range(3):
            methyl = _add_atom(mol, "C")
            mol.AddBond(central, methyl, Chem.BondType.SINGLE)
        return central
    if normalized in {"Et", "C2H5"}:
        first = _add_atom(mol, "C")
        second = _add_atom(mol, "C")
        mol.AddBond(first, second, Chem.BondType.SINGLE)
        return first
    if normalized in {"OMe", "MeO"}:
        oxygen = _add_atom(mol, "O")
        methyl = _add_atom(mol, "C")
        mol.AddBond(oxygen, methyl, Chem.BondType.SINGLE)
        return oxygen
    if normalized in {"OEt", "EtO"}:
        oxygen = _add_atom(mol, "O")
        carbon_a = _add_atom(mol, "C")
        carbon_b = _add_atom(mol, "C")
        mol.AddBond(oxygen, carbon_a, Chem.BondType.SINGLE)
        mol.AddBond(carbon_a, carbon_b, Chem.BondType.SINGLE)
        return oxygen
    if normalized in {"OMs", "MsO"}:
        return _add_methanesulfonyl(mol, "O")
    if normalized in {"OTs", "TsO", "OTos", "TosO"}:
        return _add_tosyl(mol, "O")
    if normalized in {"NMs", "NHMs", "MsN", "MsNH"}:
        return _add_methanesulfonyl(mol, "N")
    if normalized in {"NTs", "NHTs", "TsN", "TsNH", "NTos", "NHTos", "TosN", "TosNH"}:
        return _add_tosyl(mol, "N")
    if normalized in {"NO2", "O2N"}:
        nitrogen = _add_atom(mol, "N", charge=1)
        oxygen_double = _add_atom(mol, "O")
        oxygen_single = _add_atom(mol, "O", charge=-1)
        mol.AddBond(nitrogen, oxygen_double, Chem.BondType.DOUBLE)
        mol.AddBond(nitrogen, oxygen_single, Chem.BondType.SINGLE)
        return nitrogen
    if normalized in {"CO2H", "COOH"}:
        carbonyl = _add_atom(mol, "C")
        oxygen_double = _add_atom(mol, "O")
        oxygen_single = _add_atom(mol, "O")
        mol.AddBond(carbonyl, oxygen_double, Chem.BondType.DOUBLE)
        mol.AddBond(carbonyl, oxygen_single, Chem.BondType.SINGLE)
        return carbonyl
    if normalized in {"CO2Me", "COOMe"}:
        carbonyl = _add_atom(mol, "C")
        oxygen_double = _add_atom(mol, "O")
        oxygen_single = _add_atom(mol, "O")
        methyl = _add_atom(mol, "C")
        mol.AddBond(carbonyl, oxygen_double, Chem.BondType.DOUBLE)
        mol.AddBond(carbonyl, oxygen_single, Chem.BondType.SINGLE)
        mol.AddBond(oxygen_single, methyl, Chem.BondType.SINGLE)
        return carbonyl
    if normalized in {"CO2Et", "COOEt"}:
        carbonyl = _add_atom(mol, "C")
        oxygen_double = _add_atom(mol, "O")
        oxygen_single = _add_atom(mol, "O")
        carbon_a = _add_atom(mol, "C")
        carbon_b = _add_atom(mol, "C")
        mol.AddBond(carbonyl, oxygen_double, Chem.BondType.DOUBLE)
        mol.AddBond(carbonyl, oxygen_single, Chem.BondType.SINGLE)
        mol.AddBond(oxygen_single, carbon_a, Chem.BondType.SINGLE)
        mol.AddBond(carbon_a, carbon_b, Chem.BondType.SINGLE)
        return carbonyl
    if normalized == "CN":
        carbon = _add_atom(mol, "C")
        nitrogen = _add_atom(mol, "N")
        mol.AddBond(carbon, nitrogen, Chem.BondType.TRIPLE)
        return carbon
    if normalized == "Ph":
        atoms = [_add_atom(mol, "C", aromatic=True) for _ in range(6)]
        for idx in range(6):
            mol.AddBond(atoms[idx], atoms[(idx + 1) % 6], Chem.BondType.AROMATIC)
        return atoms[0]
    return None


def _molscribe_graph_to_smiles(output: dict[str, Any]) -> tuple[str | None, str | None, list[str]]:
    try:
        from rdkit import Chem
    except Exception as exc:
        return None, None, [f"RDKit unavailable in MolScribe adapter fallback: {exc}"]

    atoms = output.get("atoms") or []
    bonds = output.get("bonds") or []
    if not atoms or not bonds:
        return None, None, ["MolScribe graph fallback skipped because atoms or bonds are missing."]

    mol = Chem.RWMol()
    index_map: dict[int, int] = {}
    warnings: list[str] = []
    for idx, atom_payload in enumerate(atoms):
        label = str(atom_payload.get("atom_symbol", "")).strip()
        if not label:
            warnings.append(f"MolScribe atom {idx} had an empty label.")
            continue
        expanded = _add_abbreviation(mol, label)
        if expanded is not None:
            index_map[idx] = expanded
            continue
        normalized = _normalized_label(label)
        try:
            index_map[idx] = _add_atom(mol, normalized)
        except Exception:
            warnings.append(f"Unsupported MolScribe atom label in fallback graph: {label!r}.")

    for bond_payload in bonds:
        endpoints = bond_payload.get("endpoint_atoms") or []
        if len(endpoints) != 2:
            warnings.append(f"Skipping malformed MolScribe bond endpoints: {endpoints!r}.")
            continue
        try:
            begin = index_map[int(endpoints[0])]
            end = index_map[int(endpoints[1])]
        except Exception:
            warnings.append(f"Skipping MolScribe bond with unresolved endpoints: {endpoints!r}.")
            continue
        if begin == end or mol.GetBondBetweenAtoms(begin, end):
            continue
        mol.AddBond(begin, end, _bond_type(str(bond_payload.get("bond_type", "single"))))

    try:
        candidate = mol.GetMol()
        Chem.SanitizeMol(candidate)
        return Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True), Chem.MolToMolBlock(candidate), warnings
    except Exception as exc:
        warnings.append(f"MolScribe graph fallback could not sanitize candidate: {exc}")
        return None, None, warnings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run MolScribe on one molecule image and print JSON for codex-chem.")
    parser.add_argument("--input", required=True, help="Path to a cropped molecule image.")
    parser.add_argument("--checkpoint", help="Path to a MolScribe checkpoint. Defaults to CODEX_CHEM_MOLSCRIBE_CHECKPOINT or Hugging Face download.")
    parser.add_argument("--device", default=os.environ.get("CODEX_CHEM_MOLSCRIBE_DEVICE", "cpu"))
    args = parser.parse_args(argv)

    image_path = Path(args.input).expanduser().resolve()
    if not image_path.exists():
        print(f"Input image does not exist: {image_path}", file=sys.stderr)
        return 2

    try:
        import torch
        from molscribe import MolScribe

        model = MolScribe(_checkpoint_path(args.checkpoint), device=torch.device(args.device))
        output = model.predict_image_file(
            str(image_path),
            return_atoms_bonds=True,
            return_confidence=True,
        )
    except Exception as exc:
        print(f"MolScribe adapter failed: {exc}", file=sys.stderr)
        return 1

    smiles = output.get("smiles")
    molfile = output.get("molfile")
    adapter_warnings: list[str] = []
    if not smiles or smiles == "<invalid>":
        fallback_smiles, fallback_molfile, adapter_warnings = _molscribe_graph_to_smiles(output)
        if fallback_smiles:
            smiles = fallback_smiles
            molfile = fallback_molfile or molfile

    payload = {
        "smiles": smiles,
        "canonical_smiles": smiles,
        "molfile": molfile,
        "confidence": output.get("confidence"),
        "atoms": output.get("atoms"),
        "bonds": output.get("bonds"),
        "adapter_warnings": adapter_warnings,
        "source": "molscribe",
        "image_path": str(image_path),
    }
    print(json.dumps(payload, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
