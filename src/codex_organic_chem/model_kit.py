"""A molecular model kit: 3D questions answered from conformer ensembles.

A chemist picks up a model to see which face of a carbonyl is open, whether two
atoms can reach each other, whether a ring system is flat or whether a drawn
ring fusion can exist. These functions answer the same questions from an
ETKDG/MMFF conformer ensemble (optionally refined with GFN2-xTB) and report the
numbers per conformer so the conclusion can be checked.

Scope: sterics and ground-state geometry only. Nothing here locates a
transition state; torsional, electronic, chelation and solvent effects on
selectivity are outside the model and must be argued separately.

Atom indices are those of the input SMILES in reading order (hydrogens written
explicitly in the SMILES count; implicit hydrogens are appended after all
heavy atoms). Call :func:`atom_table` to see them, or pass atom-map numbers
with ``by_map=True``.
"""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .rdkit_tools import RDKIT_AVAILABLE
from .structure_compare import StructureError, canonical, formula_counter, hill_formula

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through tests
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdCIPLabeler, rdMolTransforms
else:  # pragma: no cover
    np = None
    Chem = None
    AllChem = None
    rdCIPLabeler = None
    rdMolTransforms = None


GAS_CONSTANT_KCAL = 1.987204258640832e-3
TEMPERATURE_K = 298.15
HARTREE_TO_KCAL = 627.5094740631
DEFAULT_SEED = 20261115
SCOPE_NOTE = (
    "Steric/ground-state model (ETKDG + MMFF94 conformers). No transition state, "
    "electronic, torsional, chelation or explicit-solvent effect is included."
)


@dataclass
class Ensemble:
    mol: Any  # hydrogens added
    conf_ids: list[int]
    energies: list[float]  # kcal/mol, absolute force-field energies
    weights: list[float]  # Boltzmann weights at 298 K
    force_field: str
    warnings: list[str]

    @property
    def relative(self) -> list[float]:
        low = min(self.energies)
        return [e - low for e in self.energies]

    def positions(self, conf_id: int) -> Any:
        return np.asarray(self.mol.GetConformer(conf_id).GetPositions(), dtype=float)


def _parse_keep_order(smiles: str) -> Any:
    if not RDKIT_AVAILABLE:
        raise StructureError("RDKit is not available.")
    mol = Chem.MolFromSmiles(smiles.strip()) if isinstance(smiles, str) else None
    if mol is None:
        raise StructureError(f"RDKit cannot parse SMILES: {smiles!r}")
    return mol


def atom_table(smiles: str) -> dict[str, Any]:
    """Atom indices, map numbers and neighbours for addressing atoms."""
    mol = _parse_keep_order(smiles)
    rows = []
    for atom in mol.GetAtoms():
        rows.append(
            {
                "index": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "map": atom.GetAtomMapNum() or None,
                "hybridization": str(atom.GetHybridization()).lower(),
                "n_h": atom.GetTotalNumHs(),
                "charge": atom.GetFormalCharge(),
                "neighbors": [n.GetIdx() for n in atom.GetNeighbors()],
                "in_ring": atom.IsInRing(),
            }
        )
    return {"status": "ok", "smiles": smiles, "atoms": rows}


def _resolve_atoms(mol: Any, atoms: Sequence[int], by_map: bool) -> list[int]:
    resolved: list[int] = []
    for ref in atoms:
        ref = int(ref)
        if by_map:
            hits = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == ref]
            if len(hits) != 1:
                raise StructureError(f"Atom-map number {ref} matches {len(hits)} atoms.")
            resolved.append(hits[0])
        else:
            if ref < 0 or ref >= mol.GetNumAtoms():
                raise StructureError(f"Atom index {ref} is out of range (0..{mol.GetNumAtoms() - 1}).")
            resolved.append(ref)
    return resolved


def build_ensemble(
    smiles: str,
    *,
    n_confs: int | None = None,
    seed: int = DEFAULT_SEED,
    max_iters: int = 500,
) -> Ensemble:
    """ETKDGv3 conformers, MMFF94s-optimized, with Boltzmann weights."""
    base = _parse_keep_order(smiles)
    warnings: list[str] = []
    from .structure_compare import parse as _strict_parse, undefined_stereo_count

    open_count, _ = undefined_stereo_count(_strict_parse(smiles))
    if open_count:
        warnings.append(
            f"{open_count} stereo element(s) are undefined; the embedding picks configurations "
            "arbitrarily, so 3D conclusions may describe a different stereoisomer."
        )
    mol = Chem.AddHs(base)
    n_rot = AllChem.CalcNumRotatableBonds(base)
    if n_confs is None:
        n_confs = min(200, 20 + 8 * n_rot)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.35
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, int(n_confs), params))
    if not conf_ids:
        params.useRandomCoords = True
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, int(n_confs), params))
    if not conf_ids:
        raise StructureError(
            "No 3D conformer could be embedded. With defined stereochemistry this usually means the "
            "drawn configuration is geometrically impossible or extremely strained."
        )
    force_field = "MMFF94s"
    energies: list[float] = []
    if AllChem.MMFFHasAllMoleculeParams(mol):
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=max_iters, mmffVariant="MMFF94s")
        energies = [float(energy) for _, energy in results]
    else:
        force_field = "UFF"
        warnings.append("MMFF94 parameters are missing for this structure; UFF was used.")
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iters)
        energies = [float(energy) for _, energy in results]
    low = min(energies)
    beta = 1.0 / (GAS_CONSTANT_KCAL * TEMPERATURE_K)
    raw = [math.exp(-(e - low) * beta) for e in energies]
    total = sum(raw)
    weights = [value / total for value in raw]
    return Ensemble(mol, conf_ids, energies, weights, force_field, warnings)


def _ensemble_summary(ens: Ensemble) -> dict[str, Any]:
    rel = ens.relative
    return {
        "n_conformers": len(ens.conf_ids),
        "force_field": ens.force_field,
        "n_within_3_kcal": sum(1 for value in rel if value <= 3.0),
        "temperature_k": TEMPERATURE_K,
    }


def _ranked(ens: Ensemble, limit: int = 8) -> list[int]:
    """Positions (not conformer ids) ordered by energy."""
    return sorted(range(len(ens.conf_ids)), key=lambda i: ens.energies[i])[:limit]


# --------------------------------------------------------------------------
# distances and dihedrals


def distance_profile(
    smiles: str,
    atom_a: int,
    atom_b: int,
    *,
    threshold: float = 4.0,
    by_map: bool = False,
    n_confs: int | None = None,
) -> dict[str, Any]:
    """Can two atoms reach each other? Distance over the conformer ensemble."""
    ens = build_ensemble(smiles, n_confs=n_confs)
    a, b = _resolve_atoms(ens.mol, [atom_a, atom_b], by_map)
    rel = ens.relative
    distances = []
    for pos, conf_id in enumerate(ens.conf_ids):
        xyz = ens.positions(conf_id)
        distances.append(float(np.linalg.norm(xyz[a] - xyz[b])))
    closest = min(range(len(distances)), key=lambda i: distances[i])
    within = [i for i, d in enumerate(distances) if d <= threshold]
    result = {
        "status": "ok",
        "task": "distance",
        "atoms": [a, b],
        "atom_symbols": [ens.mol.GetAtomWithIdx(a).GetSymbol(), ens.mol.GetAtomWithIdx(b).GetSymbol()],
        "threshold_angstrom": threshold,
        "distance_in_lowest_conformer": round(distances[_ranked(ens, 1)[0]], 2),
        "boltzmann_mean_distance": round(sum(w * d for w, d in zip(ens.weights, distances)), 2),
        "closest_approach": {
            "distance": round(distances[closest], 2),
            "relative_energy_kcal": round(rel[closest], 2),
        },
        "population_within_threshold": round(sum(ens.weights[i] for i in within), 4),
        "lowest_energy_conformer_within_threshold_kcal": (
            round(min(rel[i] for i in within), 2) if within else None
        ),
        "ensemble": _ensemble_summary(ens),
        "warnings": list(ens.warnings),
        "scope": SCOPE_NOTE,
    }
    if not within:
        result["interpretation"] = (
            f"No sampled conformer brings the atoms within {threshold} A; closest approach is "
            f"{distances[closest]:.2f} A. Absence from a finite sample is not proof of impossibility."
        )
    else:
        penalty = min(rel[i] for i in within)
        result["interpretation"] = (
            f"A conformer within {threshold} A exists {penalty:.1f} kcal/mol above the minimum "
            "(force-field energy); reaching it is a pre-organisation cost, not a barrier."
        )
    return result


def dihedral_profile(
    smiles: str,
    atoms: Sequence[int],
    *,
    by_map: bool = False,
    n_confs: int | None = None,
) -> dict[str, Any]:
    """Dihedral a-b-c-d over the ensemble: anti/syn-periplanar populations."""
    if len(atoms) != 4:
        raise StructureError("A dihedral needs exactly four atom references.")
    ens = build_ensemble(smiles, n_confs=n_confs)
    a, b, c, d = _resolve_atoms(ens.mol, atoms, by_map)
    rel = ens.relative
    angles = [
        float(rdMolTransforms.GetDihedralDeg(ens.mol.GetConformer(cid), a, b, c, d)) for cid in ens.conf_ids
    ]

    def population(test: Any) -> tuple[float, float | None]:
        idx = [i for i, angle in enumerate(angles) if test(angle)]
        return (
            round(sum(ens.weights[i] for i in idx), 4),
            round(min(rel[i] for i in idx), 2) if idx else None,
        )

    anti, anti_e = population(lambda x: abs(x) >= 150.0)
    syn, syn_e = population(lambda x: abs(x) <= 30.0)
    gauche, gauche_e = population(lambda x: 30.0 < abs(x) < 90.0)
    lowest = _ranked(ens, 5)
    return {
        "status": "ok",
        "task": "dihedral",
        "atoms": [a, b, c, d],
        "dihedral_in_lowest_conformer": round(angles[lowest[0]], 1),
        "anti_periplanar": {"population": anti, "lowest_relative_energy_kcal": anti_e},
        "syn_periplanar": {"population": syn, "lowest_relative_energy_kcal": syn_e},
        "gauche": {"population": gauche, "lowest_relative_energy_kcal": gauche_e},
        "low_energy_conformers": [
            {"relative_energy_kcal": round(rel[i], 2), "dihedral": round(angles[i], 1)} for i in lowest
        ],
        "ensemble": _ensemble_summary(ens),
        "warnings": list(ens.warnings),
        "scope": SCOPE_NOTE,
    }


# --------------------------------------------------------------------------
# facial shielding of a trigonal centre


def _grid(radius: float, spacing: float) -> Any:
    axis = np.arange(-radius, radius + spacing / 2, spacing)
    gx, gy, gz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    norms = np.linalg.norm(points, axis=1)
    return points[(norms <= radius) & (norms >= 1.2)]


def _vdw_radii(mol: Any) -> Any:
    table = Chem.GetPeriodicTable()
    return np.array([table.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()], dtype=float)


def _heavy_parent(mol: Any, idx: int) -> int:
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() == 1 and atom.GetDegree():
        return atom.GetNeighbors()[0].GetIdx()
    return idx


def _pairwise_trigonal_cip_gt(base_mol: Any, centre: int, nb_a: int, nb_b: int, nb_drop: int) -> int:
    """Compare two trigonal substituents ``nb_a`` vs ``nb_b`` around ``centre`` using IUPAC 2013 CIP rules.

    Returns ``+1`` if ``nb_a > nb_b``, ``-1`` if ``nb_b > nb_a``, or ``0`` if CIP-indistinguishable.
    """
    from .structure_compare import _permutation_parity

    rw = Chem.RWMol(base_mol)
    rw.RemoveAllConformers()
    if rw.GetAtomWithIdx(centre).GetIsAromatic():
        try:
            Chem.Kekulize(rw, clearAromaticFlags=True)
        except Exception:
            pass
    centre_z = rw.GetAtomWithIdx(centre).GetAtomicNum() or 6
    # Open any multiple bond at centre into a single bond + zero-H CIP duplicate atom on the partner.
    for nb in (nb_a, nb_b):
        bond = rw.GetBondBetweenAtoms(centre, nb)
        if bond is None:
            return 0
        mult = int(round(bond.GetBondTypeAsDouble()))
        if mult > 1:
            bond.SetBondType(Chem.BondType.SINGLE)
            for _ in range(mult - 1):
                dup = Chem.Atom(centre_z)
                dup.SetNoImplicit(True)
                dup.SetNumExplicitHs(0)
                dup_idx = rw.AddAtom(dup)
                rw.AddBond(nb, dup_idx, Chem.BondType.SINGLE)
    # Detach the third substituent and attach high-Z (Og=118) and zero-Z (*) reference ligands.
    if rw.GetBondBetweenAtoms(centre, nb_drop) is not None:
        rw.RemoveBond(centre, nb_drop)
    drop_atom = rw.GetAtomWithIdx(nb_drop)
    drop_atom.SetNoImplicit(True)
    c_atom = rw.GetAtomWithIdx(centre)
    c_atom.SetAtomicNum(6)
    c_atom.SetFormalCharge(0)
    c_atom.SetNoImplicit(True)
    c_atom.SetNumExplicitHs(0)
    hi = Chem.Atom(118)
    hi.SetNoImplicit(True)
    hi.SetNumExplicitHs(0)
    lo = Chem.Atom(0)
    lo.SetNoImplicit(True)
    lo.SetNumExplicitHs(0)
    hi_idx = rw.AddAtom(hi)
    lo_idx = rw.AddAtom(lo)
    rw.AddBond(centre, hi_idx, Chem.BondType.SINGLE)
    rw.AddBond(centre, lo_idx, Chem.BondType.SINGLE)

    # Set CW tag on centre; in RDKit, looking from first neighbor toward remaining three, CW means clockwise.
    rw.GetAtomWithIdx(centre).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    probe = rw.GetMol()
    try:
        Chem.SanitizeMol(probe)
        rdCIPLabeler.AssignCIPLabels(probe)
    except Exception:
        return 0
    p_centre = probe.GetAtomWithIdx(centre)
    if not p_centre.HasProp("_CIPCode"):
        return 0
    cip = p_centre.GetProp("_CIPCode")
    bond_neighbors = [b.GetOtherAtomIdx(centre) for b in p_centre.GetBonds()]
    # Reference order [hi_idx, nb_a, nb_b, lo_idx]: looking down centre->lo_idx,
    # if nb_a > nb_b then priorities 1(hi) -> 2(nb_a) -> 3(nb_b) match CW ("R") when parity == 0.
    # Let's check against RDKit's bond-order convention: in [hi_idx, nb_a, nb_b, lo_idx],
    # CHI_TETRAHEDRAL_CW means looking from bond_neighbors[0], [1]->[2]->[3] is clockwise.
    # For target [lo_idx, hi_idx, nb_a, nb_b] (looking from lo_idx, hi->nb_a->nb_b is CW -> viewing from front, hi->nb_a->nb_b is CCW = S!):
    parity = _permutation_parity(bond_neighbors, [lo_idx, hi_idx, nb_a, nb_b])
    # When parity == 0, CW tag on bond_neighbors means looking from lo_idx, hi->nb_a->nb_b is CW,
    # which viewed from the front (away from lo_idx) is CCW = S if nb_a > nb_b, or R if nb_b > nb_a!
    is_s = cip in {"S", "s"}
    if parity == 1:
        is_s = not is_s
    return 1 if is_s else -1


def _trigonal_cip_order(mol: Any, centre: int, neighbors: list[int]) -> tuple[list[int], bool]:
    """Return ``(ordered_neighbors_descending_cip, all_three_distinct)`` using IUPAC 2013 CIP rules."""
    n0, n1, n2 = neighbors
    c01 = _pairwise_trigonal_cip_gt(mol, centre, n0, n1, n2)
    c12 = _pairwise_trigonal_cip_gt(mol, centre, n1, n2, n0)
    c02 = _pairwise_trigonal_cip_gt(mol, centre, n0, n2, n1)
    if 0 in (c01, c12, c02):
        return sorted(neighbors), False
    wins = {
        n0: (1 if c01 > 0 else 0) + (1 if c02 > 0 else 0),
        n1: (1 if c01 < 0 else 0) + (1 if c12 > 0 else 0),
        n2: (1 if c02 < 0 else 0) + (1 if c12 < 0 else 0),
    }
    if set(wins.values()) != {0, 1, 2}:
        return sorted(neighbors), False
    ordered = sorted(neighbors, key=lambda idx: wins[idx], reverse=True)
    return ordered, True


def face_accessibility(
    smiles: str,
    atom: int,
    *,
    by_map: bool = False,
    n_confs: int | None = None,
    radius: float = 3.5,
    cone_half_angle: float = 60.0,
) -> dict[str, Any]:
    """Which face of a trigonal (sp2) atom is sterically more open?

    For each conformer the space within ``radius`` of the atom and inside a
    cone around each face normal is sampled on a grid; ``shielding`` is the
    fraction of that volume occupied by other atoms' van der Waals spheres.
    Faces are named Re/Si when the three substituents have distinct CIP ranks.
    """
    ens = build_ensemble(smiles, n_confs=n_confs)
    (centre,) = _resolve_atoms(ens.mol, [atom], by_map)
    mol = ens.mol
    centre_atom = mol.GetAtomWithIdx(centre)
    neighbors = [n.GetIdx() for n in centre_atom.GetNeighbors()]
    if len(neighbors) != 3:
        raise StructureError(
            f"Atom {centre} ({centre_atom.GetSymbol()}) has {len(neighbors)} neighbours; "
            "face analysis needs a trigonal atom with exactly three (e.g. a carbonyl or alkene carbon)."
        )
    ordered, distinct = _trigonal_cip_order(mol, centre, neighbors)
    # From the +normal side the sequence ordered[0] -> [1] -> [2] runs
    # counter-clockwise; with CIP-descending order that face is Si.
    plus_name, minus_name = ("Si", "Re") if distinct else ("A", "B")
    radii = _vdw_radii(mol)
    max_r = float(radii.max()) if radii.size else 2.0
    grid = _grid(radius, 0.25)
    grid_dirs = grid / np.linalg.norm(grid, axis=1)[:, None]
    cos_limit = math.cos(math.radians(cone_half_angle))
    all_others = np.array([i for i in range(mol.GetNumAtoms()) if i != centre], dtype=int)
    per_conf = []
    shield_atoms: dict[str, dict[int, float]] = {plus_name: {}, minus_name: {}}
    for pos, conf_id in enumerate(ens.conf_ids):
        xyz = ens.positions(conf_id)
        origin = xyz[centre]
        v0, v1, v2 = (xyz[i] - origin for i in ordered)
        normal = np.cross(v0, v1) + np.cross(v1, v2) + np.cross(v2, v0)
        length = float(np.linalg.norm(normal))
        if length < 1e-6:
            continue
        normal /= length
        cosines = grid_dirs @ normal
        plus_mask = cosines >= cos_limit
        minus_mask = cosines <= -cos_limit
        cone_mask = plus_mask | minus_mask
        cone_points = origin + grid[cone_mask]
        # Cull atoms too far to intersect the sampling sphere around origin
        other_xyz = xyz[all_others]
        near_filter = np.linalg.norm(other_xyz - origin[None, :], axis=1) <= (radius + max_r)
        others = all_others[near_filter]
        if others.size and cone_points.shape[0]:
            delta = cone_points[:, None, :] - xyz[others][None, :, :]
            dist = np.linalg.norm(delta, axis=2)
            inside = dist < radii[others][None, :]
            buried = inside.any(axis=1)
        else:
            dist = np.empty((cone_points.shape[0], 0), dtype=float)
            inside = np.zeros((cone_points.shape[0], 0), dtype=bool)
            buried = np.zeros((cone_points.shape[0],), dtype=bool)
        sub_plus = plus_mask[cone_mask]
        sub_minus = minus_mask[cone_mask]
        values = {}
        for name, mask in ((plus_name, sub_plus), (minus_name, sub_minus)):
            total = int(mask.sum())
            values[name] = float(buried[mask].sum()) / total if total else 0.0
            sub = inside[mask]
            if sub.size:
                masked_dist = np.where(sub, dist[mask], np.inf)
                nearest = masked_dist.argmin(axis=1)
                hit = sub.any(axis=1)
                for column in nearest[hit]:
                    atom_idx = int(others[int(column)])
                    parent = _heavy_parent(mol, atom_idx)
                    if parent in neighbors or parent == centre:
                        continue
                    bucket = shield_atoms[name]
                    bucket[parent] = bucket.get(parent, 0.0) + ens.weights[pos] / max(total, 1)
        per_conf.append({"pos": pos, **values})
    if not per_conf:
        raise StructureError("The centre is not planar enough to define two faces.")
    rel = ens.relative
    mean = {
        name: sum(ens.weights[row["pos"]] * row[name] for row in per_conf) for name in (plus_name, minus_name)
    }
    weight_total = sum(ens.weights[row["pos"]] for row in per_conf)
    mean = {name: value / weight_total for name, value in mean.items()}
    more_open = min(mean, key=lambda k: mean[k])
    agree = sum(
        ens.weights[row["pos"]] for row in per_conf if row[more_open] < row[plus_name if more_open == minus_name else minus_name]
    ) / weight_total
    difference = abs(mean[plus_name] - mean[minus_name])
    if difference < 0.04:
        verdict = "no clear steric facial bias"
    elif difference < 0.10:
        verdict = f"modest steric preference for the {more_open} face"
    else:
        verdict = f"clear steric preference for the {more_open} face"

    def describe(face: str) -> list[dict[str, Any]]:
        items = sorted(shield_atoms[face].items(), key=lambda kv: kv[1], reverse=True)[:3]
        return [
            {
                "heavy_atom_index": idx,
                "symbol": mol.GetAtomWithIdx(idx).GetSymbol(),
                "share_of_cone_volume": round(share, 3),
            }
            for idx, share in items
            if share >= 0.005
        ]

    ordered_conf = sorted(per_conf, key=lambda row: rel[row["pos"]])[:6]
    result = {
        "status": "ok",
        "task": "face_accessibility",
        "atom": centre,
        "atom_symbol": centre_atom.GetSymbol(),
        "substituents_by_priority" if distinct else "substituents": [
            {"index": idx, "symbol": mol.GetAtomWithIdx(idx).GetSymbol()} for idx in ordered
        ],
        "face_definition": (
            "Re/Si from CIP priorities of the three substituents."
            if distinct
            else (
                "Two substituents have equal CIP rank, so Re/Si is undefined. Face A is the side from "
                f"which atoms {ordered[0]} -> {ordered[1]} -> {ordered[2]} run counter-clockwise."
            )
        ),
        "shielding": {name: round(value, 3) for name, value in mean.items()},
        "more_open_face": more_open if difference >= 0.04 else None,
        "shielding_difference": round(difference, 3),
        "population_agreeing": round(agree, 3),
        "verdict": verdict,
        "shielding_atoms": {name: describe(name) for name in (plus_name, minus_name)},
        "low_energy_conformers": [
            {
                "relative_energy_kcal": round(rel[row["pos"]], 2),
                plus_name: round(row[plus_name], 3),
                minus_name: round(row[minus_name], 3),
            }
            for row in ordered_conf
        ],
        "ensemble": _ensemble_summary(ens),
        "warnings": list(ens.warnings),
        "scope": SCOPE_NOTE
        + " Addition can still favour the more shielded face (torsional/Felkin-Anh, chelation, "
        "small nucleophiles on cyclohexanones).",
    }
    return result


# --------------------------------------------------------------------------
# planarity, embedding feasibility and isomer energies


def _plane_rmsd(xyz: Any) -> tuple[float, float]:
    centred = xyz - xyz.mean(axis=0)
    _, singular, vt = np.linalg.svd(centred)
    normal = vt[-1]
    deviations = centred @ normal
    return float(np.sqrt((deviations**2).mean())), float(np.abs(deviations).max())


def _planarity_label(rmsd: float) -> str:
    if rmsd < 0.05:
        return "planar"
    if rmsd < 0.20:
        return "nearly planar"
    return "nonplanar"


def _xtb_optimize(mol: Any, conf_id: int, solvent: str | None, timeout_s: int) -> tuple[Any, float] | None:
    from .semiempirical import _base_args, _read_xtb_json
    from .utils import executable_status, run_command

    if executable_status("xtb", ("--version",)).status != "available":
        return None
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    unpaired = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    with tempfile.TemporaryDirectory(prefix="codex-chem-kit-xtb-") as tmp:
        work = Path(tmp)
        (work / "input.xyz").write_text(Chem.MolToXYZBlock(mol, confId=conf_id), encoding="utf-8")
        try:
            code, _, _ = run_command(
                ["xtb", "input.xyz", "--opt", "--gfn", "2", "--json", *_base_args(charge, unpaired, solvent)],
                timeout_s=timeout_s,
                cwd=str(work),
            )
        except Exception:
            return None
        out = work / "xtbopt.xyz"
        if code != 0 or not out.exists():
            return None
        lines = out.read_text(encoding="utf-8", errors="replace").splitlines()[2:]
        coords = np.array([[float(x) for x in line.split()[1:4]] for line in lines if len(line.split()) >= 4])
        energy = _read_xtb_json(work).get("total energy")
    if coords.shape[0] != mol.GetNumAtoms() or energy is None:
        return None
    return coords, float(energy)


def planarity(
    smiles: str,
    atoms: Sequence[int] | None = None,
    *,
    by_map: bool = False,
    refine: str = "auto",
    solvent: str | None = None,
    n_confs: int | None = None,
    timeout_s: int = 240,
) -> dict[str, Any]:
    """Is a set of atoms (default: every ring system and all heavy atoms) flat?

    ``refine='xtb'`` re-optimizes the lowest conformer with GFN2-xTB before
    measuring; ``'auto'`` does so when xTB is installed, because force fields
    are unreliable for the flat/twisted balance of conjugated systems.
    """
    ens = build_ensemble(smiles, n_confs=n_confs)
    mol = ens.mol
    lowest = ens.conf_ids[_ranked(ens, 1)[0]]
    xyz = ens.positions(lowest)
    geometry = ens.force_field
    warnings = list(ens.warnings)
    if refine in {"auto", "xtb"}:
        refined = _xtb_optimize(mol, lowest, solvent, timeout_s)
        if refined is not None:
            xyz, _ = refined
            geometry = "GFN2-xTB (optimized from the lowest force-field conformer)"
        elif refine == "xtb":
            warnings.append("xTB refinement was requested but is unavailable or failed; force-field geometry used.")
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    selections: list[tuple[str, list[int]]] = []
    if atoms:
        selections.append(("selected_atoms", _resolve_atoms(mol, atoms, by_map)))
    else:
        selections.append(("all_heavy_atoms", heavy))
        ring_info = mol.GetRingInfo()
        systems: list[set[int]] = []
        for ring in ring_info.AtomRings():
            ring_set = set(ring)
            merged = [s for s in systems if s & ring_set]
            for s in merged:
                ring_set |= s
                systems.remove(s)
            systems.append(ring_set)
        for i, system in enumerate(sorted(systems, key=lambda s: min(s)), start=1):
            selections.append((f"ring_system_{i}", sorted(system)))
        for i, ring in enumerate(ring_info.AtomRings(), start=1):
            if len(systems) != len(ring_info.AtomRings()):
                selections.append((f"ring_{i}_size_{len(ring)}", list(ring)))
    report = []
    for label, idx in selections:
        if len(idx) < 4:
            continue
        rmsd, max_dev = _plane_rmsd(xyz[idx])
        report.append(
            {
                "selection": label,
                "atoms": idx,
                "rms_deviation_from_plane": round(rmsd, 3),
                "max_deviation": round(max_dev, 3),
                "label": _planarity_label(rmsd),
            }
        )
    return {
        "status": "ok",
        "task": "planarity",
        "geometry": geometry,
        "selections": report,
        "thresholds": "rms deviation < 0.05 A planar; < 0.20 A nearly planar; otherwise nonplanar",
        "ensemble": _ensemble_summary(ens),
        "warnings": warnings,
        "scope": "Geometry of one optimized minimum in the gas phase (or ALPB solvent if given).",
    }


def embed_check(smiles: str) -> dict[str, Any]:
    """Can the structure exist in 3D with the stereochemistry as drawn?"""
    mol = _parse_keep_order(smiles)
    result: dict[str, Any] = {"status": "ok", "task": "embed_check", "smiles": canonical(Chem.RemoveHs(mol))}
    try:
        ens = build_ensemble(smiles, n_confs=10)
    except StructureError as exc:
        flat = Chem.Mol(mol)
        Chem.RemoveStereochemistry(flat)
        embeds_without = True
        try:
            build_ensemble(Chem.MolToSmiles(flat), n_confs=5)
        except StructureError:
            embeds_without = False
        result.update(
            {
                "embeds": False,
                "embeds_without_stereo": embeds_without,
                "interpretation": (
                    "The constitution embeds but the drawn configuration does not: check the ring "
                    "fusion/bridgehead wedges."
                    if embeds_without
                    else "Neither the drawn nor the stereo-free structure embeds; the ring system is "
                    "probably impossible or far too strained."
                ),
                "detail": str(exc),
            }
        )
        return result
    result.update(
        {
            "embeds": True,
            "ensemble": _ensemble_summary(ens),
            "warnings": list(ens.warnings),
            "interpretation": (
                "A 3D model exists. This does not rank strain; use isomer_energies to compare "
                "stereoisomers or ring-fusion alternatives."
            ),
        }
    )
    return result


def isomer_energies(
    smiles_list: Sequence[str],
    *,
    method: str = "mmff",
    solvent: str | None = None,
    n_confs: int | None = None,
    timeout_s: int = 240,
) -> dict[str, Any]:
    """Relative energies of isomers (same formula and charge) at their best conformer.

    Use for thermodynamic questions: which diastereomer/ring fusion/alkene
    geometry is more stable, i.e. where an equilibrating step ends up. Not a
    kinetic prediction.
    """
    if method not in {"mmff", "xtb"}:
        raise StructureError("method must be 'mmff' or 'xtb'.")
    if len(smiles_list) < 2:
        raise StructureError("Provide at least two isomers.")
    rows = []
    reference: tuple[str, int] | None = None
    warnings: list[str] = []
    for smi in smiles_list:
        mol = _parse_keep_order(smi)
        signature = (
            hill_formula(formula_counter(Chem.RemoveHs(mol))),
            sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        )
        if reference is None:
            reference = signature
        elif signature != reference:
            raise StructureError(
                f"{smi!r} has formula/charge {signature}, not {reference}: energies of different "
                "compositions are not comparable."
            )
        try:
            ens = build_ensemble(smi, n_confs=n_confs)
        except StructureError as exc:
            rows.append({"smiles": smi, "status": "no_3d_model", "detail": str(exc)})
            continue
        warnings.extend(f"{smi}: {w}" for w in ens.warnings)
        row: dict[str, Any] = {
            "smiles": smi,
            "status": "ok",
            "n_conformers": len(ens.conf_ids),
            "energy_kcal": min(ens.energies),
            "method": ens.force_field,
        }
        if method == "xtb":
            best: float | None = None
            for pos in _ranked(ens, 3):
                refined = _xtb_optimize(ens.mol, ens.conf_ids[pos], solvent, timeout_s)
                if refined is not None:
                    energy = refined[1] * HARTREE_TO_KCAL
                    best = energy if best is None else min(best, energy)
            if best is None:
                row["status"] = "xtb_failed"
            else:
                row["energy_kcal"] = best
                row["method"] = "GFN2-xTB" + (f"/ALPB({solvent})" if solvent else "")
        rows.append(row)
    usable = [r for r in rows if r["status"] == "ok"]
    methods = {r["method"] for r in usable}
    if len(methods) > 1:
        warnings.append("Isomers were evaluated with different methods; relative energies are not comparable.")
    if usable:
        low = min(r["energy_kcal"] for r in usable)
        for r in usable:
            r["relative_energy_kcal"] = round(r["energy_kcal"] - low, 2)
            r.pop("energy_kcal")
    for r in rows:
        r.pop("energy_kcal", None)
    ranked = sorted(usable, key=lambda r: r["relative_energy_kcal"])
    gap = ranked[1]["relative_energy_kcal"] if len(ranked) > 1 else None
    resolution = 2.0 if method == "mmff" else 1.5
    return {
        "status": "ok",
        "task": "isomer_energies",
        "formula": reference[0] if reference else None,
        "isomers": rows,
        "most_stable": ranked[0]["smiles"] if ranked else None,
        "gap_to_next_kcal": gap,
        "decisive": bool(gap is not None and gap >= resolution),
        "warnings": warnings,
        "scope": (
            f"Lowest-conformer {'GFN2-xTB electronic' if method == 'xtb' else 'force-field'} energies; "
            f"gaps below ~{resolution} kcal/mol are within method error. No entropy, no explicit solvent, "
            "thermodynamic preference only."
        ),
    }


def ring_relationship(
    smiles: str,
    atom_a: int,
    atom_b: int,
    *,
    by_map: bool = False,
    n_confs: int | None = None,
) -> dict[str, Any]:
    """cis or trans: are the substituents on two ring atoms on the same face?

    The question behind every wedge on a ring. For each ring atom the
    exocyclic substituent of highest CIP rank is used (so H is used only when
    nothing else is there), and "same face" is judged against the local ring
    plane at each atom, which stays meaningful in puckered and large rings.
    Reported per conformer population so a flexible macrocycle shows up as
    such instead of as a coin toss.
    """
    ens = build_ensemble(smiles, n_confs=n_confs)
    mol = ens.mol
    a, b = _resolve_atoms(mol, [atom_a, atom_b], by_map)
    rings = [list(ring) for ring in mol.GetRingInfo().AtomRings() if a in ring and b in ring]
    if a == b or not rings:
        raise StructureError("cis/trans needs two different atoms of the same ring.")
    ring = min(rings, key=len)
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)

    def substituent(idx: int) -> int:
        outside = [n for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetIdx() not in ring]
        if not outside:
            raise StructureError(f"Ring atom {idx} has no substituent outside the ring.")
        ranked = sorted(
            outside,
            key=lambda n: (n.GetAtomicNum() > 1, int(n.GetProp("_CIPRank")) if n.HasProp("_CIPRank") else 0),
        )
        if len(outside) > 1 and all(n.GetAtomicNum() == 1 for n in outside):
            raise StructureError(f"Ring atom {idx} carries two hydrogens; it has no face to compare.")
        return ranked[-1].GetIdx()

    sub_a, sub_b = substituent(a), substituent(b)

    def side(xyz: Any, idx: int, sub: int) -> float:
        position = ring.index(idx)
        before, after = ring[position - 1], ring[(position + 1) % len(ring)]
        normal = np.cross(xyz[after] - xyz[idx], xyz[before] - xyz[idx])
        return float(np.dot(normal, xyz[sub] - xyz[idx]))

    population_cis = 0.0
    n_cis = 0
    for conf_id, weight in zip(ens.conf_ids, ens.weights):
        xyz = mol.GetConformer(conf_id).GetPositions()
        if side(xyz, a, sub_a) * side(xyz, b, sub_b) > 0:
            population_cis += weight
            n_cis += 1
    if population_cis >= 0.9:
        relationship = "cis"
    elif population_cis <= 0.1:
        relationship = "trans"
    else:
        relationship = "conformer_dependent"

    def name(idx: int) -> str:
        return f"{mol.GetAtomWithIdx(idx).GetSymbol()}{idx}"

    return {
        "status": "ok",
        "task": "cis_trans",
        "smiles": smiles,
        "ring_atoms": [a, b],
        "ring_size": len(ring),
        "substituents_compared": {name(a): name(sub_a), name(b): name(sub_b)},
        "relationship": relationship,
        "population_cis": round(population_cis, 3),
        "conformers_cis": f"{n_cis}/{len(ens.conf_ids)}",
        "note": (
            "Substituents are the highest-priority exocyclic neighbours. In rings larger than about seven "
            "atoms 'face' is a conformational property; trust this only when every conformer agrees."
        ),
        "ensemble": _ensemble_summary(ens),
        "warnings": list(ens.warnings),
    }


TASKS = ("atoms", "distance", "dihedral", "face", "cis_trans", "planarity", "embed", "isomer_energies")


def model_kit(task: str, smiles: str | Sequence[str], **options: Any) -> dict[str, Any]:
    """Dispatch used by the CLI and MCP tool; errors come back as data."""
    try:
        if task == "atoms":
            return atom_table(str(smiles))
        if task == "distance":
            atoms = options.pop("atoms")
            return distance_profile(str(smiles), atoms[0], atoms[1], **options)
        if task == "dihedral":
            return dihedral_profile(str(smiles), options.pop("atoms"), **options)
        if task == "face":
            atoms = options.pop("atoms")
            return face_accessibility(str(smiles), atoms[0], **options)
        if task == "cis_trans":
            atoms = options.pop("atoms")
            return ring_relationship(str(smiles), atoms[0], atoms[1], **options)
        if task == "planarity":
            return planarity(str(smiles), options.pop("atoms", None), **options)
        if task == "embed":
            return embed_check(str(smiles))
        if task == "isomer_energies":
            items = [smiles] if isinstance(smiles, str) else list(smiles)
            return isomer_energies(items, **options)
    except StructureError as exc:
        return {"status": "error", "task": task, "warnings": [str(exc)]}
    except (KeyError, IndexError, TypeError) as exc:
        return {"status": "error", "task": task, "warnings": [f"Missing or malformed option: {exc}"]}
    return {"status": "error", "task": task, "warnings": [f"Unknown task {task!r}; use one of {TASKS}."]}
