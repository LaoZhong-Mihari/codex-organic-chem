"""Structure identity, stereochemistry and enumeration checks.

These functions exist to give *discriminating* feedback: an answer that is
wrong in a way a chemist would care about (wrong diastereomer, miscounted
isomers, a stereocentre left undefined) must produce a different result from a
right one. Parsing and canonicalising a structure the model already wrote
correctly carries no such information.

Everything here is graph-based RDKit. Stereo comparisons are done under an
explicit atom mapping and compare spatial arrangement, not CIP letters: a CIP
label can change when a substituent changes although nothing moved in space.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Iterable

from .rdkit_tools import RDKIT_AVAILABLE

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through tests
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdCIPLabeler, rdMolDescriptors
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
else:  # pragma: no cover
    Chem = None
    AllChem = None
    rdCIPLabeler = None
    rdMolDescriptors = None
    EnumerateStereoisomers = None
    StereoEnumerationOptions = None


class StructureError(ValueError):
    """Raised when an input structure cannot be parsed."""


MAX_MAPPINGS = 512
MAX_STEREOISOMERS = 1024


# --------------------------------------------------------------------------
# parsing helpers


def parse(smiles: str, *, keep_maps: bool = False) -> Any:
    """Parse SMILES into a sanitized mol with explicit H removed.

    Atom-map numbers are stripped unless ``keep_maps`` is set, because they
    change canonical SMILES and therefore identity comparisons.
    """
    if not RDKIT_AVAILABLE:
        raise StructureError("RDKit is not available.")
    if not isinstance(smiles, str) or not smiles.strip():
        raise StructureError("Empty SMILES.")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise StructureError(f"RDKit cannot parse SMILES: {smiles!r}")
    mol = Chem.RemoveHs(mol)
    if not keep_maps:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def formula_counter(mol: Any) -> Counter:
    counter: Counter = Counter()
    for atom in mol.GetAtoms():
        counter[atom.GetSymbol()] += 1
        counter["H"] += atom.GetTotalNumHs()
    return +counter


def hill_formula(counter: Counter) -> str:
    """Hill-order formula; negative counts are kept so deltas stay readable."""
    pieces = []
    keys = [k for k in ("C", "H") if counter.get(k)] + sorted(
        k for k in counter if k not in {"C", "H"} and counter[k]
    )
    for key in keys:
        n = counter[key]
        pieces.append(f"{key}{n if n != 1 else ''}")
    return "".join(pieces)


def delta_text(delta: Counter | dict[str, int]) -> str:
    items = [(k, v) for k, v in delta.items() if v]
    if not items:
        return "0"
    order = {"C": 0, "H": 1}
    items.sort(key=lambda kv: (order.get(kv[0], 2), kv[0]))
    return " ".join(f"{v:+d}{k}" for k, v in items)


def canonical(mol: Any, *, stereo: bool = True) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereo)


def mirror(mol: Any) -> Any:
    """Mirror image: invert every tetrahedral centre; E/Z is unchanged."""
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        if atom.GetChiralTag() in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        ):
            atom.InvertChirality()
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return copy


def _neighbor_text(atom: Any) -> str:
    symbols = sorted(n.GetSymbol() for n in atom.GetNeighbors())
    symbols += ["H"] * atom.GetTotalNumHs()
    return ",".join(symbols)


def _cip_labels(mol: Any) -> None:
    try:
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:  # pragma: no cover - very large or unusual structures
        pass


# --------------------------------------------------------------------------
# stereo inventory


def _cumulated(atom: Any) -> bool:
    """sp centre: two double bonds or a triple bond (azide, allene, ketene)."""
    orders = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
    return orders.count(2.0) >= 2 or 3.0 in orders


def potential_stereo(mol: Any) -> list[Any]:
    """RDKit's potential stereo elements, minus double bonds that cannot have
    E/Z because one end is linear. RDKit reports the N=N of some alkyl azides
    as an open double bond, which would make every azide look undefined."""
    kept = []
    for info in Chem.FindPotentialStereo(mol, cleanIt=False):
        if info.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(info.centeredOn)
            if _cumulated(bond.GetBeginAtom()) or _cumulated(bond.GetEndAtom()):
                continue
        kept.append(info)
    return kept


def stereo_elements(mol: Any) -> list[dict[str, Any]]:
    """Every potential stereo element, specified or not."""
    _cip_labels(mol)
    elements: list[dict[str, Any]] = []
    for info in potential_stereo(mol):
        specified = info.specified == Chem.StereoSpecified.Specified
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            atom = mol.GetAtomWithIdx(info.centeredOn)
            elements.append(
                {
                    "kind": "tetrahedral",
                    "atom_index": atom.GetIdx(),
                    "element": atom.GetSymbol(),
                    "substituents": _neighbor_text(atom),
                    "specified": specified,
                    "cip": atom.GetProp("_CIPCode") if specified and atom.HasProp("_CIPCode") else None,
                    "in_ring": atom.IsInRing(),
                }
            )
        elif info.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(info.centeredOn)
            label = None
            if specified:
                if bond.HasProp("_CIPCode"):
                    label = bond.GetProp("_CIPCode")
                elif bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):
                    label = "E"
                elif bond.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):
                    label = "Z"
            elements.append(
                {
                    "kind": "double_bond",
                    "atom_indices": [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
                    "specified": specified,
                    "cip": label,
                    "in_ring": bond.IsInRing(),
                }
            )
    return elements


def consistent_stereoisomers(mol: Any, *, max_elements: int = 8) -> list[Any] | None:
    """Distinct stereoisomers compatible with what the drawing defines.

    A bridgehead of norbornane is a *potential* stereocentre that the ring
    system fixes; leaving it unmarked does not make the structure ambiguous.
    Enumerating only the undefined elements (and discarding isomers that cannot
    be embedded in 3D) separates that case from a genuinely open centre.
    Returns ``None`` when there are too many open elements to enumerate.
    """
    open_elements = [
        info
        for info in potential_stereo(mol)
        if info.specified != Chem.StereoSpecified.Specified
    ]
    if not open_elements:
        return [mol]
    if len(open_elements) > max_elements:
        return None
    in_ring = any(
        mol.GetAtomWithIdx(info.centeredOn).IsInRing()
        for info in open_elements
        if info.type == Chem.StereoType.Atom_Tetrahedral
    )
    options = StereoEnumerationOptions(
        unique=True,
        onlyUnassigned=True,
        maxIsomers=MAX_STEREOISOMERS,
        tryEmbedding=in_ring and len(open_elements) <= 6,
    )
    unique: dict[str, Any] = {}
    for iso in EnumerateStereoisomers(mol, options=options):
        Chem.AssignStereochemistry(iso, cleanIt=True, force=True)
        unique.setdefault(canonical(iso), iso)
    return list(unique.values()) or [mol]


def resolve_constrained_stereo(mol: Any) -> tuple[Any, int | None]:
    """Complete ring-fixed centres; return (mol, number of consistent isomers)."""
    isomers = consistent_stereoisomers(mol)
    if isomers is None:
        return mol, None
    if len(isomers) == 1:
        return isomers[0], 1
    return mol, len(isomers)


def undefined_stereo_count(mol: Any) -> tuple[int, int | None]:
    """(number of genuinely open stereo elements, isomers consistent with the drawing)."""
    resolved, n_consistent = resolve_constrained_stereo(mol)
    if n_consistent == 1:
        return 0, 1
    open_elements = sum(1 for e in stereo_elements(resolved) if not e["specified"])
    return open_elements, n_consistent


def _is_chiral(mol: Any) -> bool:
    return canonical(mol) != canonical(mirror(mol))


def stereo_audit(smiles: str, *, enumerate_isomers: bool = True) -> dict[str, Any]:
    """Inventory of stereo elements plus what the drawn structure commits to.

    Use this after reading a structure from an image: every wedge/hash you saw
    must appear here as a specified centre with the configuration you intended,
    and ``unspecified`` must be empty unless the source really leaves it open.
    """
    drawn = parse(smiles)
    mol, n_consistent = resolve_constrained_stereo(drawn)
    elements = stereo_elements(mol)
    drawn_open = sum(1 for e in stereo_elements(drawn) if not e["specified"])
    if n_consistent == 1:
        # Whatever is still unmarked cannot be varied (cubane, adamantane):
        # it is a symmetry artefact of perception, not an open stereocentre.
        elements = [e for e in elements if e["specified"]]
    unspecified = [e for e in elements if not e["specified"]]
    centres = [e for e in elements if e["kind"] == "tetrahedral"]
    bonds = [e for e in elements if e["kind"] == "double_bond"]
    result: dict[str, Any] = {
        "status": "ok",
        "input_smiles": smiles,
        "canonical_smiles": canonical(mol),
        "formula": hill_formula(formula_counter(mol)),
        "n_stereocenters": len(centres),
        "n_stereo_double_bonds": len(bonds),
        "n_undefined": len(unspecified),
        "n_isomers_consistent_with_drawing": n_consistent,
        "fully_defined": not unspecified,
        "elements": elements,
        "warnings": [],
    }
    if drawn_open and not unspecified and elements:
        result["note"] = (
            f"{drawn_open} unmarked centre(s) are fixed by the ring system; the drawing denotes one "
            "stereoisomer and they were completed in canonical_smiles."
        )
    if unspecified:
        count = f"{n_consistent} stereoisomers" if n_consistent else "many stereoisomers"
        result["warnings"].append(
            f"{len(unspecified)} stereo element(s) are undefined: the structure is consistent with "
            f"{count}, not a single compound."
        )
    else:
        chiral = _is_chiral(mol)
        result["chiral"] = chiral
        result["meso"] = bool(centres) and not chiral
    if enumerate_isomers and elements:
        result["stereoisomers"] = _stereoisomer_summary(mol)
    return result


def _stereoisomer_summary(mol: Any, *, limit: int = 32) -> dict[str, Any]:
    base = Chem.Mol(mol)
    Chem.RemoveStereochemistry(base)
    n_elements = len(list(Chem.FindPotentialStereo(base)))
    if n_elements > 10:
        return {
            "status": "skipped",
            "reason": f"{n_elements} stereo elements (up to {2 ** n_elements} isomers); not enumerated.",
        }
    has_ring_stereo = any(
        mol.GetAtomWithIdx(i.centeredOn).IsInRing()
        for i in Chem.FindPotentialStereo(base)
        if i.type == Chem.StereoType.Atom_Tetrahedral
    )
    # Ring-fused/bridged systems admit graph stereoisomers that cannot exist in
    # 3D (e.g. trans-bridgeheads of norbornane); embedding filters those out.
    try_embedding = has_ring_stereo and n_elements <= 6
    options = StereoEnumerationOptions(
        unique=True, onlyUnassigned=False, maxIsomers=MAX_STEREOISOMERS, tryEmbedding=try_embedding
    )
    isomers = list(EnumerateStereoisomers(base, options=options))
    smiles_set = {}
    for iso in isomers:
        Chem.AssignStereochemistry(iso, cleanIt=True, force=True)
        smiles_set[canonical(iso)] = iso
    achiral, pairs, seen = [], [], set()
    for smi, iso in smiles_set.items():
        if smi in seen:
            continue
        partner = canonical(mirror(iso))
        if partner == smi:
            achiral.append(smi)
            seen.add(smi)
        else:
            pairs.append([smi, partner])
            seen.update({smi, partner})
    summary: dict[str, Any] = {
        "status": "ok",
        "n_stereoisomers": len(smiles_set),
        "n_enantiomer_pairs": len(pairs),
        "n_achiral_isomers": len(achiral),
        "geometric_filter": "3D-embeddable isomers only" if try_embedding else "graph enumeration",
    }
    if len(smiles_set) <= limit:
        summary["enantiomer_pairs"] = pairs
        summary["achiral_isomers"] = achiral
    return summary


def enumerate_stereoisomers(smiles: str, *, only_unassigned: bool = False) -> dict[str, Any]:
    """Count and list distinct stereoisomers, pairing enantiomers and meso forms."""
    mol = parse(smiles)
    if only_unassigned:
        options = StereoEnumerationOptions(unique=True, onlyUnassigned=True, maxIsomers=MAX_STEREOISOMERS)
        isomers = sorted({canonical(m) for m in EnumerateStereoisomers(mol, options=options)})
        return {
            "status": "ok",
            "input_smiles": smiles,
            "mode": "assigned centres kept, unassigned enumerated",
            "n_stereoisomers": len(isomers),
            "stereoisomers": isomers[:64],
        }
    return {
        "status": "ok",
        "input_smiles": smiles,
        "mode": "all stereo elements enumerated (drawn configuration ignored)",
        **_stereoisomer_summary(mol, limit=64),
    }


# --------------------------------------------------------------------------
# mapping-based configuration comparison


def _permutation_parity(reference: list[int], other: list[int]) -> int:
    """0 if ``other`` is an even permutation of ``reference``, else 1."""
    position = {value: idx for idx, value in enumerate(reference)}
    perm = [position[value] for value in other]
    seen = [False] * len(perm)
    parity = 0
    for start in range(len(perm)):
        if seen[start]:
            continue
        length = 0
        node = start
        while not seen[node]:
            seen[node] = True
            node = perm[node]
            length += 1
        parity ^= (length - 1) & 1
    return parity


_TETRAHEDRAL = None


def _tetrahedral_tags() -> tuple[Any, Any]:
    return (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


def tetrahedral_relation(
    atom_a: Any,
    atom_b: Any,
    mapping: dict[int, int],
    *,
    pair_unmapped: bool = False,
) -> tuple[str, dict[str, Any]]:
    """Compare the spatial configuration of two mapped tetrahedral atoms.

    ``mapping`` sends atom indices of molecule A to molecule B. Returns one of
    ``same``, ``inverted``, ``unspecified_a``, ``unspecified_b``,
    ``unspecified_both`` or ``ambiguous``. With ``pair_unmapped`` a single
    replaced substituent on each side is paired (leaving group <-> incoming
    group), so ``same`` then means retention and ``inverted`` inversion.
    """
    tags = _tetrahedral_tags()
    tag_a, tag_b = atom_a.GetChiralTag(), atom_b.GetChiralTag()
    detail: dict[str, Any] = {}
    if tag_a not in tags and tag_b not in tags:
        return "unspecified_both", detail
    if tag_a not in tags:
        return "unspecified_a", detail
    if tag_b not in tags:
        return "unspecified_b", detail
    nbrs_a = [n.GetIdx() for n in atom_a.GetNeighbors()]
    nbrs_b = [n.GetIdx() for n in atom_b.GetNeighbors()]
    if len(nbrs_a) != len(nbrs_b):
        return "ambiguous", {"reason": "different coordination number"}
    image: list[int | None] = []
    for idx in nbrs_a:
        target = mapping.get(idx)
        image.append(target if target in nbrs_b else None)
    missing_a = [nbrs_a[i] for i, value in enumerate(image) if value is None]
    used = {value for value in image if value is not None}
    missing_b = [idx for idx in nbrs_b if idx not in used]
    if len(used) != len(nbrs_a) - len(missing_a):
        return "ambiguous", {"reason": "mapping is not one-to-one around the centre"}
    if missing_a:
        if not pair_unmapped or len(missing_a) != len(missing_b):
            return "ambiguous", {"reason": "neighbours changed and could not be paired"}
        mol_a, mol_b = atom_a.GetOwningMol(), atom_b.GetOwningMol()
        hs_a = [i for i in missing_a if mol_a.GetAtomWithIdx(i).GetAtomicNum() == 1]
        hs_b = [i for i in missing_b if mol_b.GetAtomWithIdx(i).GetAtomicNum() == 1]
        rest_a = [i for i in missing_a if i not in hs_a]
        rest_b = [i for i in missing_b if i not in hs_b]
        pairs: list[tuple[int, int]] = []
        if len(hs_a) == len(hs_b) == 1:
            pairs.append((hs_a[0], hs_b[0]))
        else:
            rest_a, rest_b = missing_a, missing_b
        if len(rest_a) != len(rest_b) or len(rest_a) > 1:
            return "ambiguous", {"reason": "more than one substituent replaced at this centre"}
        if rest_a:
            pairs.append((rest_a[0], rest_b[0]))
            detail["replaced"] = (
                f"{mol_a.GetAtomWithIdx(rest_a[0]).GetSymbol()}"
                f"->{mol_b.GetAtomWithIdx(rest_b[0]).GetSymbol()}"
            )
        lookup = dict(pairs)
        image = [value if value is not None else lookup[nbrs_a[i]] for i, value in enumerate(image)]
    parity = _permutation_parity(nbrs_b, [int(v) for v in image])
    same_tag = tag_a == tag_b
    return ("same" if same_tag == (parity == 0) else "inverted"), detail


def _stereo_bond_cis(bond: Any) -> bool | None:
    stereo = bond.GetStereo()
    if stereo in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):
        return True
    if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):
        return False
    return None


def double_bond_relation(bond_a: Any, mol_b: Any, mapping: dict[int, int]) -> str:
    """``same``, ``inverted``, ``unspecified_*``, ``absent`` or ``ambiguous``."""
    a1, a2 = bond_a.GetBeginAtomIdx(), bond_a.GetEndAtomIdx()
    b1, b2 = mapping.get(a1), mapping.get(a2)
    if b1 is None or b2 is None:
        return "ambiguous"
    bond_b = mol_b.GetBondBetweenAtoms(b1, b2)
    if bond_b is None or bond_b.GetBondType() != Chem.BondType.DOUBLE:
        return "absent"
    cis_a, cis_b = _stereo_bond_cis(bond_a), _stereo_bond_cis(bond_b)
    if cis_a is None and cis_b is None:
        return "unspecified_both"
    if cis_a is None:
        return "unspecified_a"
    if cis_b is None:
        return "unspecified_b"
    sa = list(bond_a.GetStereoAtoms())
    sb = list(bond_b.GetStereoAtoms())
    if len(sa) != 2 or len(sb) != 2:
        return "ambiguous"
    # stereo atoms are (neighbour of begin, neighbour of end)
    ref_b1, ref_b2 = (sb[0], sb[1]) if bond_b.GetBeginAtomIdx() == b1 else (sb[1], sb[0])
    m1, m2 = mapping.get(sa[0]), mapping.get(sa[1])
    if m1 is None or m2 is None:
        return "ambiguous"
    cis_b_mapped = cis_b ^ (m1 != ref_b1) ^ (m2 != ref_b2)
    return "same" if cis_a == cis_b_mapped else "inverted"


def _full_mappings(mol_a: Any, mol_b: Any) -> list[dict[int, int]]:
    query = Chem.Mol(mol_a)
    Chem.RemoveStereochemistry(query)
    matches = mol_b.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=MAX_MAPPINGS)
    n = mol_a.GetNumAtoms()
    return [dict(enumerate(match)) for match in matches if len(match) == n and len(set(match)) == n]


def configuration_differences(mol_a: Any, mol_b: Any) -> dict[str, Any] | None:
    """Fewest stereo differences over all constitution-preserving mappings."""
    mappings = _full_mappings(mol_a, mol_b)
    if not mappings:
        return None
    _cip_labels(mol_a)
    _cip_labels(mol_b)
    centres = [
        e["atom_index"] for e in stereo_elements(mol_a) if e["kind"] == "tetrahedral"
    ]
    bonds = [
        mol_a.GetBondBetweenAtoms(*e["atom_indices"]).GetIdx()
        for e in stereo_elements(mol_a)
        if e["kind"] == "double_bond"
    ]
    best: dict[str, Any] | None = None
    for mapping in mappings:
        differing, open_items = [], []
        for idx in centres:
            atom_a = mol_a.GetAtomWithIdx(idx)
            atom_b = mol_b.GetAtomWithIdx(mapping[idx])
            relation, _ = tetrahedral_relation(atom_a, atom_b, mapping)
            entry = {
                "kind": "tetrahedral",
                "atom_index_a": idx,
                "atom_index_b": mapping[idx],
                "element": atom_a.GetSymbol(),
                "substituents": _neighbor_text(atom_a),
                "cip_a": atom_a.GetProp("_CIPCode") if atom_a.HasProp("_CIPCode") else None,
                "cip_b": atom_b.GetProp("_CIPCode") if atom_b.HasProp("_CIPCode") else None,
            }
            if relation == "inverted":
                differing.append(entry)
            elif relation.startswith("unspecified") and relation != "unspecified_both":
                open_items.append({**entry, "state": relation})
            elif relation == "ambiguous":
                open_items.append({**entry, "state": relation})
        for bond_idx in bonds:
            bond_a = mol_a.GetBondWithIdx(bond_idx)
            relation = double_bond_relation(bond_a, mol_b, mapping)
            entry = {
                "kind": "double_bond",
                "atom_indices_a": [bond_a.GetBeginAtomIdx(), bond_a.GetEndAtomIdx()],
                "atom_indices_b": [mapping[bond_a.GetBeginAtomIdx()], mapping[bond_a.GetEndAtomIdx()]],
            }
            if relation == "inverted":
                differing.append(entry)
            elif relation not in {"same", "unspecified_both"}:
                open_items.append({**entry, "state": relation})
        score = (len(differing), len(open_items))
        if best is None or score < best["score"]:
            best = {"score": score, "differing": differing, "open": open_items, "mapping": mapping}
            if score == (0, 0):
                break
    return best


def compare_structures(smiles_a: str, smiles_b: str) -> dict[str, Any]:
    """Classify the relationship between two structures.

    ``relation`` is one of ``identical``, ``enantiomers``, ``diastereomers``,
    ``same_constitution_stereo_undefined``, ``tautomers``,
    ``constitutional_isomers`` or ``different_formula``. For diastereomers the
    differing centres are located on both molecules (atom indices of the input
    SMILES order), so the error can be traced to a specific wedge.
    """
    mol_a, _ = resolve_constrained_stereo(parse(smiles_a))
    mol_b, _ = resolve_constrained_stereo(parse(smiles_b))
    fa, fb = formula_counter(mol_a), formula_counter(mol_b)
    result: dict[str, Any] = {
        "status": "ok",
        "a": {"smiles": canonical(mol_a), "formula": hill_formula(fa)},
        "b": {"smiles": canonical(mol_b), "formula": hill_formula(fb)},
        "warnings": [],
    }
    charge_a = sum(a.GetFormalCharge() for a in mol_a.GetAtoms())
    charge_b = sum(a.GetFormalCharge() for a in mol_b.GetAtoms())
    if fa != fb or charge_a != charge_b:
        delta = Counter(fb)
        delta.subtract(fa)
        result["relation"] = "different_formula"
        result["formula_delta_b_minus_a"] = delta_text({k: v for k, v in delta.items() if v})
        if charge_a != charge_b:
            result["charge_delta_b_minus_a"] = charge_b - charge_a
        return result
    if canonical(mol_a, stereo=False) != canonical(mol_b, stereo=False):
        relation = "constitutional_isomers"
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize

            enumerator = rdMolStandardize.TautomerEnumerator()
            if canonical(enumerator.Canonicalize(mol_a), stereo=False) == canonical(
                enumerator.Canonicalize(mol_b), stereo=False
            ):
                relation = "tautomers"
        except Exception:  # pragma: no cover - standardizer unavailable
            pass
        result["relation"] = relation
        return result
    if canonical(mol_a) == canonical(mol_b):
        result["relation"] = "identical"
        open_a = [e for e in stereo_elements(mol_a) if not e["specified"]]
        if open_a:
            result["warnings"].append(
                f"Both structures leave {len(open_a)} stereo element(s) undefined; identity holds "
                "only at that level of definition."
            )
        return result
    open_a = [e for e in stereo_elements(mol_a) if not e["specified"]]
    open_b = [e for e in stereo_elements(mol_b) if not e["specified"]]
    diff = configuration_differences(mol_a, mol_b)
    if open_a or open_b:
        result["relation"] = "same_constitution_stereo_undefined"
        result["n_undefined_a"] = len(open_a)
        result["n_undefined_b"] = len(open_b)
        if diff is not None:
            result["conflicting_defined_elements"] = diff["differing"]
            result["defined_in_only_one"] = diff["open"]
            result["compatible"] = not diff["differing"]
        return result
    if canonical(mirror(mol_a)) == canonical(mol_b):
        result["relation"] = "enantiomers"
        result["note"] = "Every stereocentre is inverted; double-bond geometry is unchanged."
        return result
    result["relation"] = "diastereomers"
    if diff is not None:
        result["differing_elements"] = diff["differing"]
        result["n_differing"] = len(diff["differing"])
        mirrored = configuration_differences(mirror(mol_a), mol_b)
        if mirrored is not None and len(mirrored["differing"]) < len(diff["differing"]):
            result["note"] = (
                f"Closer to the enantiomer of A: {len(mirrored['differing'])} element(s) differ from "
                "mirror(A). Check whether the whole drawing was mirrored before one centre was changed."
            )
    return result


# --------------------------------------------------------------------------
# symmetry, environments, substitution counting


def _embedded_with_hs(mol: Any, seed: int = 61453) -> Any | None:
    with_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(with_h, params) < 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(with_h, params) < 0:
            return None
    return with_h


def _substituted_smiles(with_h: Any, h_idx: int, atomic_num: int, isotope: int = 0) -> str:
    """Replace one explicit H and read the product's stereo from 3D."""
    rw = Chem.RWMol(with_h)
    atom = rw.GetAtomWithIdx(h_idx)
    atom.SetAtomicNum(atomic_num)
    atom.SetIsotope(isotope)
    product = rw.GetMol()
    Chem.SanitizeMol(product)
    if product.GetNumConformers():
        Chem.AssignStereochemistryFrom3D(product)
    product = Chem.RemoveHs(product)
    return canonical(product)


def symmetry_environments(smiles: str) -> dict[str, Any]:
    """Symmetry-distinct atoms: expected 13C and 1H signal counts.

    Hydrogens are classified by the substitution test (replace each H by a
    label and compare products): identical products -> homotopic, mirror-image
    products -> enantiotopic (one signal in an achiral medium), otherwise
    diastereotopic/constitutionally distinct (separate signals in principle).
    Fast conformational averaging is implied because the test is graph-based.
    Accidental overlap and slow exchange are not modelled.
    """
    mol = parse(smiles)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
    heavy_classes: dict[str, dict[int, list[int]]] = {}
    for atom in mol.GetAtoms():
        heavy_classes.setdefault(atom.GetSymbol(), {}).setdefault(ranks[atom.GetIdx()], []).append(atom.GetIdx())
    result: dict[str, Any] = {
        "status": "ok",
        "input_smiles": smiles,
        "canonical_smiles": canonical(mol),
        "heavy_atom_environments": {
            symbol: {"count": len(groups), "classes": sorted(groups.values())}
            for symbol, groups in sorted(heavy_classes.items())
        },
        "n_carbon_environments": len(heavy_classes.get("C", {})),
        "warnings": [],
    }
    with_h = _embedded_with_hs(mol)
    if with_h is None:
        result["warnings"].append("3D embedding failed; hydrogen environments were not classified.")
        return result
    exact: dict[str, list[int]] = {}
    for atom in with_h.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        parent = atom.GetNeighbors()[0].GetIdx()
        try:
            key = _substituted_smiles(with_h, atom.GetIdx(), 1, isotope=3)
        except Exception:
            key = f"unresolved-{atom.GetIdx()}"
        exact.setdefault(key, []).append(parent)
    merged: dict[str, dict[str, Any]] = {}
    for key, parents in exact.items():
        try:
            partner = canonical(mirror(Chem.MolFromSmiles(key)))
        except Exception:
            partner = key
        group = min(key, partner)
        entry = merged.setdefault(group, {"n_h": 0, "on_atoms": set(), "enantiotopic": False})
        entry["n_h"] += len(parents)
        entry["on_atoms"].update(parents)
        if partner != key and partner in exact:
            entry["enantiotopic"] = True
    environments = [
        {
            "n_h": entry["n_h"],
            "on_heavy_atoms": sorted(entry["on_atoms"]),
            "contains_enantiotopic_h": entry["enantiotopic"],
        }
        for entry in merged.values()
    ]
    environments.sort(key=lambda e: (e["on_heavy_atoms"], e["n_h"]))
    appearances: Counter = Counter(a for e in environments for a in e["on_heavy_atoms"])
    diastereotopic_atoms = sorted(atom for atom, count in appearances.items() if count > 1)
    result["hydrogen_environments"] = environments
    result["n_hydrogen_environments"] = len(environments)
    result["n_hydrogen_environments_chiral_medium"] = len(exact)
    result["atoms_with_diastereotopic_h"] = diastereotopic_atoms
    result["note"] = (
        "n_hydrogen_environments counts enantiotopic H as equivalent (achiral solvent); "
        "diastereotopic H on one carbon count separately and may or may not resolve."
    )
    return result


_SUBSTITUENT_ATOMS = {"F": 9, "Cl": 17, "Br": 35, "I": 53, "D": 1, "OH": 8, "NH2": 7, "CH3": 6, "Me": 6}


def substitution_products(
    smiles: str,
    substituent: str = "Cl",
    *,
    n_substitutions: int = 1,
    on_elements: Iterable[str] = ("C",),
) -> dict[str, Any]:
    """Count distinct products of replacing H by a substituent.

    Answers "how many monochloro (or dichloro) derivatives" both ways:
    constitutional isomers only, and counting stereoisomers separately.
    """
    if substituent not in _SUBSTITUENT_ATOMS:
        raise StructureError(f"Unsupported substituent {substituent!r}; use one of {sorted(_SUBSTITUENT_ATOMS)}.")
    if n_substitutions not in (1, 2):
        raise StructureError("n_substitutions must be 1 or 2.")
    mol, _ = resolve_constrained_stereo(parse(smiles))
    atomic_num = _SUBSTITUENT_ATOMS[substituent]
    isotope = 2 if substituent == "D" else 0
    allowed = set(on_elements)
    undefined = [None] * undefined_stereo_count(mol)[0]
    with_h = _embedded_with_hs(mol)
    method = "3D substitution (configuration read from an embedded model)"
    if with_h is None:
        with_h = Chem.AddHs(mol)
        method = "graph substitution (3D embedding failed; new stereocentres left undefined)"
    sites = [
        atom.GetIdx()
        for atom in with_h.GetAtoms()
        if atom.GetAtomicNum() == 1 and not atom.GetIsotope() and atom.GetNeighbors()[0].GetSymbol() in allowed
    ]
    if n_substitutions == 1:
        combos: list[tuple[int, ...]] = [(i,) for i in sites]
    else:
        combos = [(i, j) for pos, i in enumerate(sites) for j in sites[pos + 1 :]]
    stereo_all: set[str] = set()
    for combo in combos:
        rw = Chem.RWMol(with_h)
        for idx in combo:
            target = rw.GetAtomWithIdx(idx)
            target.SetAtomicNum(atomic_num)
            target.SetIsotope(isotope)
        product = rw.GetMol()
        try:
            Chem.SanitizeMol(product)
            if product.GetNumConformers():
                Chem.AssignStereochemistryFrom3D(product)
            stereo_all.add(canonical(Chem.RemoveHs(product)))
        except Exception:
            continue
    constitutional = sorted({canonical(Chem.MolFromSmiles(s), stereo=False) for s in stereo_all})
    warnings = []
    if undefined:
        warnings.append(
            f"The parent has {len(undefined)} undefined stereo element(s); one arbitrary configuration "
            "was embedded, so the stereoisomer count applies to a single parent stereoisomer."
        )
    return {
        "status": "ok",
        "input_smiles": smiles,
        "substituent": substituent,
        "n_substitutions": n_substitutions,
        "n_constitutional_isomers": len(constitutional),
        "n_including_stereoisomers": len(stereo_all),
        "constitutional_isomers": constitutional[:64],
        "stereoisomers": sorted(stereo_all)[:64],
        "method": method,
        "warnings": warnings,
        "note": (
            "Each product keeps the parent skeleton's configuration; a pair of enantiomeric products "
            "counts as two in n_including_stereoisomers."
        ),
    }
