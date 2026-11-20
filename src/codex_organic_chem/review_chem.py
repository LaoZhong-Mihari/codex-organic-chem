"""Chemistry behind the review workbench.

Two jobs: describe a structure so a person can check it at a glance (drawing
with CIP labels and clickable atoms, formula, every stereo element, how it
differs from a reference), and apply the one-click corrections people need most
often when a drawing was misread (flip one centre, set an open centre, flip a
double bond, mirror the molecule). Pure functions; the HTTP layer only calls
these.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from . import structure_compare as sc
from .structure_compare import StructureError

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdCIPLabeler, rdDepictor, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D

    RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover - RDKit is a hard dependency in practice
    RDKIT_AVAILABLE = False

EDIT_OPS = ("invert", "set", "clear", "mirror", "clear_all")

# highlight colours (r, g, b): what changed against the reference, what the
# caller asked to point at, and what is still undefined
_CHANGED = (1.0, 0.72, 0.30)
_POINTED = (0.55, 0.78, 1.0)
_OPEN = (1.0, 0.86, 0.86)


def _mol(smiles: str | None = None, molfile: str | None = None) -> Any:
    if not RDKIT_AVAILABLE:
        raise StructureError("RDKit is not available.")
    if molfile and molfile.strip():
        mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=True)
        if mol is None:
            raise StructureError("RDKit cannot parse the Molfile.")
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.AssignStereochemistryFrom3D(mol) if mol.GetConformer().Is3D() else Chem.AssignChiralTypesFromBondDirs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        return mol
    return sc.parse(smiles or "")


def working_smiles(smiles: str | None = None, molfile: str | None = None) -> str:
    """Isomeric SMILES that keeps the atom order of the input.

    Its atom order defines the indices used by :func:`describe` and
    :func:`apply_edit`, and it is written without canonical reordering so that
    "C5" stays "C5" after an edit and still means the fifth atom of the SMILES
    the assistant sent.
    """
    return Chem.MolToSmiles(_mol(smiles, molfile), canonical=False)


def _cip(mol: Any) -> None:
    try:
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:  # pragma: no cover
        pass


def _element_id(element: dict[str, Any]) -> str:
    if element["kind"] == "tetrahedral":
        return f"a{element['atom_index']}"
    i, j = element["atom_indices"]
    return f"b{i}-{j}"


def _atom_name(mol: Any, idx: int) -> str:
    return f"{mol.GetAtomWithIdx(idx).GetSymbol()}{idx + 1}"


def stereo_inventory(mol: Any) -> list[dict[str, Any]]:
    """Stereo elements with stable ids and human-readable names.

    Centres whose configuration is fixed by a ring system once the rest is
    drawn (norbornane bridgeheads) are reported as ``constrained`` rather than
    open, so nobody is asked to assign something geometry already decides.
    """
    resolved, n_consistent = sc.resolve_constrained_stereo(Chem.Mol(mol))
    fixed = n_consistent == 1
    out = []
    for element in sc.stereo_elements(mol):
        item = dict(element)
        item["id"] = _element_id(element)
        if element["kind"] == "tetrahedral":
            item["name"] = _atom_name(mol, element["atom_index"])
            item["atoms"] = [element["atom_index"]]
        else:
            i, j = element["atom_indices"]
            item["name"] = f"{_atom_name(mol, i)}={_atom_name(mol, j)}"
            item["atoms"] = [i, j]
        item["constrained"] = bool(fixed and not element["specified"])
        out.append(item)
    return out


def _bond_drawn_right(mol: Any, bond: Any) -> bool:
    """Do the 2D coordinates show this double bond with its E/Z?

    A wedge is derived from the coordinates when the picture is drawn, so any
    layout can show any configuration at a centre. The geometry of a double
    bond *is* the coordinates: a layout copied from the other isomer shows the
    other isomer, whatever the label next to it says.
    """
    same_side = (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)
    opposite = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
    stereo = bond.GetStereo()
    ends = list(bond.GetStereoAtoms())
    if stereo not in same_side + opposite or len(ends) != 2:
        return True
    conf = mol.GetConformer()
    begin = conf.GetAtomPosition(bond.GetBeginAtomIdx())
    finish = conf.GetAtomPosition(bond.GetEndAtomIdx())
    ax, ay = finish.x - begin.x, finish.y - begin.y

    def side(idx: int) -> float:
        point = conf.GetAtomPosition(idx)
        return ax * (point.y - begin.y) - ay * (point.x - begin.x)

    first, second = side(ends[0]), side(ends[1])
    if abs(first) < 1e-3 or abs(second) < 1e-3:
        return True  # collinear in this drawing: nothing to contradict
    return (first * second > 0) == (stereo in same_side)


def _mismatched_double_bonds(mol: Any) -> list[Any]:
    return [bond for bond in mol.GetBonds() if not _bond_drawn_right(mol, bond)]


def _double_bonds_drawn_right(mol: Any) -> bool:
    return not _mismatched_double_bonds(mol)


def _smaller_side(mol: Any, bond: Any) -> set[int] | None:
    """Atoms on the smaller side of an acyclic bond (None for a ring bond)."""
    if bond.IsInRing():
        return None
    begin, finish = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

    def reach(start: int, blocked: int) -> set[int]:
        seen, stack = {start}, [start]
        while stack:
            for neighbour in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                idx = neighbour.GetIdx()
                if idx not in seen and idx != blocked:
                    seen.add(idx)
                    stack.append(idx)
        return seen

    left, right = reach(begin, finish), reach(finish, begin)
    return left if len(left) <= len(right) else right


def _superpose(mol: Any, template: Any, match: tuple[int, ...]) -> None:
    """Turn (and if it fits better, mirror) a fresh drawing onto the old one."""
    import numpy as np

    conf, reference = mol.GetConformer(), template.GetConformer()
    moving = np.array([[conf.GetAtomPosition(match[i]).x, conf.GetAtomPosition(match[i]).y] for i in range(len(match))])
    fixed = np.array([[reference.GetAtomPosition(i).x, reference.GetAtomPosition(i).y] for i in range(len(match))])
    centre_moving, centre_fixed = moving.mean(axis=0), fixed.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - centre_moving).T @ (fixed - centre_fixed))
    rotation = u @ vt  # best orthogonal map; a reflection is a legitimate redrawing
    for idx in range(mol.GetNumAtoms()):
        point = conf.GetAtomPosition(idx)
        x, y = (np.array([point.x, point.y]) - centre_moving) @ rotation + centre_fixed
        conf.SetAtomPosition(idx, (float(x), float(y), 0.0))


def lay_out(mol: Any, template: Any | None = None) -> Any:
    """Copy of ``mol`` with 2D coordinates, kept as close to ``template`` (an
    earlier drawing of the same constitution) as the structure allows.

    A changed centre needs no new coordinates, so the old layout is reused
    exactly and only the wedge moves. A changed double bond does: the smaller
    side of that bond is redrawn and the rest stays put. If that cannot be done
    (a ring double bond), the molecule is redrawn and turned onto the old
    drawing.
    """
    from rdkit.Geometry import Point2D

    draw = Chem.Mol(mol)
    match = None
    if template is not None and template.GetNumConformers():
        try:
            query = Chem.Mol(template)
            Chem.RemoveStereochemistry(query)
            found = draw.GetSubstructMatch(query, useChirality=False)
            if found and len(found) == draw.GetNumAtoms():
                match = found
        except Exception:
            match = None
    if match is None:
        rdDepictor.Compute2DCoords(draw)
        return draw
    conf = template.GetConformer()
    pinned = {match[i]: Point2D(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(len(match))}
    try:
        rdDepictor.Compute2DCoords(draw, coordMap=pinned)
        wrong = _mismatched_double_bonds(draw)
        if not wrong:
            return draw
        free: set[int] = set()
        for bond in wrong:
            side = _smaller_side(draw, bond)
            if side is None:
                free = set(range(draw.GetNumAtoms()))
                break
            free |= side
        kept = {idx: point for idx, point in pinned.items() if idx not in free}
        if len(kept) >= 2:
            rdDepictor.Compute2DCoords(draw, coordMap=kept)
            if _double_bonds_drawn_right(draw):
                return draw
    except Exception:
        pass
    rdDepictor.Compute2DCoords(draw)
    try:
        _superpose(draw, template, match)
    except Exception:
        pass
    return draw


def depict(
    mol: Any,
    *,
    highlight: dict[int, tuple[float, float, float]] | None = None,
    show_indices: bool = False,
    width: int = 460,
    height: int = 340,
    template: Any | None = None,
) -> str:
    """SVG with CIP labels and one invisible ``atom-<idx>`` circle per atom, so
    the page can make atoms clickable. ``template`` keeps the layout of an
    earlier drawing of the same constitution (see :func:`lay_out`)."""
    draw = lay_out(mol, template)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.clearBackground = False
    opts.padding = 0.12
    opts.bondLineWidth = 1.6
    if width >= 400:
        # keep small molecules at a normal drawing scale instead of blowing
        # them up to fill the pane; large ones are still shrunk to fit
        opts.fixedBondLength = 34
    if show_indices:
        for atom in draw.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx() + 1))
    atoms = sorted(highlight or {})
    colours = {idx: (highlight or {})[idx] for idx in atoms}
    radii = {idx: 0.42 for idx in atoms}
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, draw, highlightAtoms=atoms, highlightAtomColors=colours, highlightAtomRadii=radii, highlightBonds=[]
    )
    drawer.TagAtoms(draw, radius=0.45)  # generous click targets, still short of the next atom
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _layout(mol: Any) -> Any:
    laid = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(laid)
    return laid


def _versus(reference: Any, mol: Any) -> dict[str, Any]:
    """How ``mol`` differs from ``reference``, with the atoms of ``mol`` to mark."""
    # written without reordering, so compare_structures reports atom indices
    # in the order of ``mol`` and they line up with the drawing
    result = sc.compare_structures(
        Chem.MolToSmiles(reference, canonical=False), Chem.MolToSmiles(mol, canonical=False)
    )
    relation = result["relation"]
    out: dict[str, Any] = {"relation": relation, "changed_atoms": [], "details": []}
    if relation == "diastereomers":
        for element in result.get("differing_elements", []):
            if element["kind"] == "tetrahedral":
                idx = element["atom_index_b"]
                out["changed_atoms"].append(idx)
                out["details"].append(
                    f"{_atom_name(mol, idx)}: {element.get('cip_a') or '?'} -> {element.get('cip_b') or '?'}"
                )
            else:
                i, j = element["atom_indices_b"]
                out["changed_atoms"] += [i, j]
                out["details"].append(f"{_atom_name(mol, i)}={_atom_name(mol, j)}: E/Z flipped")
        if result.get("note"):
            out["note"] = result["note"]
    elif relation == "enantiomers":
        out["changed_atoms"] = [e["atom_index"] for e in sc.stereo_elements(mol) if e["kind"] == "tetrahedral"]
    elif relation == "same_constitution_stereo_undefined":
        for element in result.get("conflicting_defined_elements", []):
            if element["kind"] == "tetrahedral":
                out["changed_atoms"].append(element["atom_index_b"])
                out["details"].append(f"{_atom_name(mol, element['atom_index_b'])}: inverted")
            else:
                i, j = element["atom_indices_b"]
                out["changed_atoms"] += [i, j]
                out["details"].append(f"{_atom_name(mol, i)}={_atom_name(mol, j)}: E/Z flipped")
        for element in result.get("defined_in_only_one", []):
            state = {"unspecified_a": "now defined", "unspecified_b": "no longer defined"}.get(
                element.get("state"), "changed"
            )
            if element["kind"] == "tetrahedral":
                idx = element["atom_index_b"]
                out["changed_atoms"].append(idx)
                out["details"].append(f"{_atom_name(mol, idx)}: {state}")
            else:
                i, j = element["atom_indices_b"]
                out["changed_atoms"] += [i, j]
                out["details"].append(f"{_atom_name(mol, i)}={_atom_name(mol, j)}: {state}")
    elif relation == "different_formula":
        out["formula_delta"] = result.get("formula_delta_b_minus_a")
        if result.get("charge_delta_b_minus_a"):
            out["charge_delta"] = result["charge_delta_b_minus_a"]
    return out


def describe(
    smiles: str | None = None,
    molfile: str | None = None,
    *,
    reference: str | None = None,
    highlight_atoms: list[int] | None = None,
    show_indices: bool = False,
    width: int = 460,
    height: int = 340,
) -> dict[str, Any]:
    """Everything the page shows for one structure. Errors come back as data."""
    try:
        working = working_smiles(smiles, molfile)
        mol = _mol(working)
    except StructureError as exc:
        return {"status": "invalid", "warnings": [str(exc)], "input": smiles or ""}
    except Exception as exc:  # RDKit raises plain exceptions on odd Molfiles
        return {"status": "invalid", "warnings": [f"Could not read the structure: {exc}"], "input": smiles or ""}
    _cip(mol)
    counter = sc.formula_counter(mol)
    elements = stereo_inventory(mol)
    open_elements = [e for e in elements if not e["specified"] and not e["constrained"]]
    highlight: dict[int, tuple[float, float, float]] = {}
    for element in open_elements:
        for idx in element["atoms"]:
            highlight[idx] = _OPEN
    versus = None
    template = None
    if reference:
        try:
            ref_mol = _mol(working_smiles(reference))
            versus = _versus(ref_mol, mol)
            if versus["relation"] not in {"different_formula", "constitutional_isomers", "tautomers"}:
                template = _layout(ref_mol)
            for idx in versus["changed_atoms"]:
                highlight[idx] = _CHANGED
        except StructureError:
            versus = None
    for idx in highlight_atoms or []:
        if 0 <= idx < mol.GetNumAtoms():
            highlight[idx] = _POINTED
    warnings = []
    if open_elements:
        names = ", ".join(e["name"] for e in open_elements)
        warnings.append(f"{len(open_elements)} stereo element(s) are not defined: {names}.")
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("More than one fragment (salt, solvate or two molecules).")
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if radicals:
        warnings.append(f"{radicals} unpaired electron(s) on the structure.")
    return {
        "status": "ok",
        "smiles": working,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": sc.hill_formula(counter),
        "mol_weight": round(Descriptors.MolWt(mol), 2),
        "exact_mass": round(rdMolDescriptors.CalcExactMolWt(mol), 4),
        "charge": charge,
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "stereo": elements,
        "n_stereo_elements": len(elements),
        "n_undefined": len(open_elements),
        "versus_reference": versus,
        "svg": depict(mol, highlight=highlight, show_indices=show_indices, width=width, height=height, template=template),
        "warnings": warnings,
    }


def _target_atom(mol: Any, target: Any) -> Any:
    try:
        idx = int(str(target).lstrip("a"))
    except (TypeError, ValueError) as exc:
        raise StructureError(f"Not an atom reference: {target!r}") from exc
    if not 0 <= idx < mol.GetNumAtoms():
        raise StructureError(f"Atom index {idx} is outside this structure.")
    return mol.GetAtomWithIdx(idx)


def _target_bond(mol: Any, target: Any) -> Any:
    try:
        i, j = (int(part) for part in str(target).lstrip("b").split("-"))
    except (TypeError, ValueError) as exc:
        raise StructureError(f"Not a bond reference: {target!r}") from exc
    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
        raise StructureError(f"No double bond between atoms {i} and {j}.")
    return bond


def _bond_label(mol: Any, bond: Any) -> str | None:
    # no AssignStereochemistry here: it would re-derive E/Z from bond
    # directions and discard a geometry that was just set on the bond
    probe = Chem.Mol(mol)
    _cip(probe)
    b = probe.GetBondWithIdx(bond.GetIdx())
    if b.HasProp("_CIPCode"):
        return b.GetProp("_CIPCode").upper()
    if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):
        return "E"
    if b.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):
        return "Z"
    return None


def _ensure_stereo_atoms(bond: Any) -> None:
    if len(list(bond.GetStereoAtoms())) == 2:
        return
    begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
    left = sorted(n.GetIdx() for n in begin.GetNeighbors() if n.GetIdx() != end.GetIdx())
    right = sorted(n.GetIdx() for n in end.GetNeighbors() if n.GetIdx() != begin.GetIdx())
    if not left or not right:
        raise StructureError("This double bond has no substituents to define E/Z.")
    bond.SetStereoAtoms(left[0], right[0])


def _flip_bond(bond: Any) -> None:
    _ensure_stereo_atoms(bond)
    stereo = bond.GetStereo()
    cis = (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)
    bond.SetStereo(Chem.BondStereo.STEREOTRANS if stereo in cis else Chem.BondStereo.STEREOCIS)


def _atom_label(mol: Any, idx: int) -> str | None:
    probe = Chem.Mol(mol)
    Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
    _cip(probe)
    atom = probe.GetAtomWithIdx(idx)
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def apply_edit(smiles: str, op: str, target: Any = None, value: str | None = None) -> str:
    """Apply one stereo edit to ``smiles`` (atom indices as in :func:`describe`)
    and return the new canonical SMILES.

    ``invert``  flip a centre (``a<idx>``) or a double bond (``b<i>-<j>``); an
                open element gets an arbitrary first assignment
    ``set``     set a centre to R/S or a double bond to E/Z
    ``clear``   remove the definition of one element
    ``mirror``  invert every tetrahedral centre (the enantiomer)
    ``clear_all`` remove all stereo definitions
    """
    if op not in EDIT_OPS:
        raise StructureError(f"Unknown edit {op!r}; use one of {', '.join(EDIT_OPS)}.")
    mol = Chem.RWMol(_mol(working_smiles(smiles)))
    if op == "mirror":
        return Chem.MolToSmiles(sc.mirror(mol.GetMol()), canonical=False)
    if op == "clear_all":
        plain = mol.GetMol()
        Chem.RemoveStereochemistry(plain)
        return Chem.MolToSmiles(plain, canonical=False)
    is_bond = str(target).startswith("b")
    if is_bond:
        bond = _target_bond(mol, target)
        if op == "clear":
            bond.SetStereo(Chem.BondStereo.STEREONONE)
        elif op == "invert":
            _flip_bond(bond)
        else:
            wanted = str(value or "").upper()
            if wanted not in {"E", "Z"}:
                raise StructureError("Set a double bond to E or Z.")
            _ensure_stereo_atoms(bond)
            if bond.GetStereo() == Chem.BondStereo.STEREONONE:
                bond.SetStereo(Chem.BondStereo.STEREOTRANS)
            if _bond_label(mol, bond) != wanted:
                _flip_bond(bond)
            if _bond_label(mol, bond) != wanted:
                raise StructureError(f"Could not set this double bond to {wanted}.")
    else:
        atom = _target_atom(mol, target)
        tags = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        if op == "clear":
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        elif op == "invert":
            if atom.GetChiralTag() in tags:
                atom.InvertChirality()
            else:
                atom.SetChiralTag(tags[0])
        else:
            wanted = str(value or "").upper()
            if wanted not in {"R", "S"}:
                raise StructureError("Set a centre to R or S.")
            if atom.GetChiralTag() not in tags:
                atom.SetChiralTag(tags[0])
            label = _atom_label(mol, atom.GetIdx())
            if label is None:
                raise StructureError("This atom is not a stereocentre.")
            if label.upper() != wanted:
                atom.InvertChirality()
    out = mol.GetMol()
    if not is_bond:
        Chem.SanitizeMol(out)
        Chem.AssignStereochemistry(out, cleanIt=True, force=True)
    return Chem.MolToSmiles(out, canonical=False)


def formula_delta(before: str, after: str) -> str:
    """Formula change between two SMILES, e.g. ``+1C +2H``."""
    delta = Counter(sc.formula_counter(_mol(after)))
    delta.subtract(sc.formula_counter(_mol(before)))
    return sc.delta_text({k: v for k, v in delta.items() if v})


__all__ = [
    "EDIT_OPS",
    "StructureError",
    "apply_edit",
    "depict",
    "describe",
    "formula_delta",
    "stereo_inventory",
    "working_smiles",
]
