"""Canonical molecule depiction with fixed geometry.

Every chemistry figure in this package draws molecules through
:func:`render_molecule`, so a molecule keeps the same bond length, line width,
and font size no matter which figure it appears in or how big its own cell is.

RDKit's default depiction auto-scales each molecule to fill its canvas, which
makes ethanol draw ~145 px bonds and a Boc-protected piperidine ~23 px bonds in
the same 345x300 cell. That single behaviour is what makes multi-panel figures
look inconsistent. Here the canvas is sized to the molecule instead: the caller
picks a bond length in pixels, and the canvas grows to fit.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from .rdkit_tools import RDKIT_AVAILABLE

if RDKIT_AVAILABLE:  # pragma: no cover - import availability is environment dependent
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
else:  # pragma: no cover
    Chem = None
    Draw = None
    rdDepictor = None
    rdMolDraw2D = None


# prepare_mol rescales every conformer so the median bond is exactly 1.0 mol
# unit; combined with RDKit's fixedBondLength (px per mol unit) this makes
# bond_length_px land exactly, whichever 2D embedder produced the coordinates.
MOL_UNITS_PER_BOND = 1.0

#: Named geometry presets. ``bond_length_px`` is the single most important value:
#: it fixes the drawn bond length, and every other dimension is derived from it
#: so a figure scales as one coherent unit.
FIGURE_PRESETS: dict[str, dict[str, Any]] = {
    # ACS-style: 0.508 cm bonds at ~150 dpi, restrained lines, black atoms.
    "acs": {
        "bond_length_px": 30.0,
        "bond_line_width": 1.6,
        "base_font_size": 0.62,
        "min_font_size": 9,
        "max_font_size": 15,
        "cell_padding": 16.0,
        "monochrome": True,
    },
    "rsc": {
        "bond_length_px": 28.0,
        "bond_line_width": 1.5,
        "base_font_size": 0.60,
        "min_font_size": 8,
        "max_font_size": 14,
        "cell_padding": 15.0,
        "monochrome": True,
    },
    # Slightly larger and coloured; for slides, teaching, and quick previews.
    "presentation": {
        "bond_length_px": 36.0,
        "bond_line_width": 2.0,
        "base_font_size": 0.66,
        "min_font_size": 10,
        "max_font_size": 18,
        "cell_padding": 18.0,
        "monochrome": False,
    },
}
DEFAULT_PRESET = "acs"


def resolve_preset(preset: str | None = None, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    """Merge a named preset with caller overrides, keeping unknown names safe."""
    name = str(preset or DEFAULT_PRESET).lower()
    resolved = dict(FIGURE_PRESETS.get(name, FIGURE_PRESETS[DEFAULT_PRESET]))
    resolved["name"] = name if name in FIGURE_PRESETS else DEFAULT_PRESET
    for key, value in (overrides or {}).items():
        if key in resolved and value is not None:
            resolved[key] = value
    return resolved


@dataclass
class MoleculeDepiction:
    """A drawn molecule plus the geometry an audit needs to check it.

    Coordinates are in the molecule's own SVG space; :meth:`translated_geometry`
    maps them into a composed figure's canvas space.
    """

    svg_body: str
    width: int
    height: int
    bond_lengths_px: list[float] = field(default_factory=list)
    atom_points: list[tuple[float, float]] = field(default_factory=list)
    label_boxes: list[tuple[float, float, float, float]] = field(default_factory=list)
    atom_label_boxes: dict[int, list[tuple[float, float, float, float]]] = field(default_factory=dict)
    ink_box: tuple[float, float, float, float] | None = None
    warnings: list[str] = field(default_factory=list)

    def translated_geometry(self, dx: float, dy: float) -> dict[str, Any]:
        return {
            "bond_lengths_px": list(self.bond_lengths_px),
            "atom_points": [(x + dx, y + dy) for x, y in self.atom_points],
            "label_boxes": [(x0 + dx, y0 + dy, x1 + dx, y1 + dy) for x0, y0, x1, y1 in self.label_boxes],
            "ink_box": (
                (self.ink_box[0] + dx, self.ink_box[1] + dy, self.ink_box[2] + dx, self.ink_box[3] + dy)
                if self.ink_box
                else None
            ),
        }


_ATOM_PATH = re.compile(r"<path\s+class='atom-\d+'\s+d='([^']*)'")
_NUMBER = re.compile(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?")


def _bond_lengths_from_atoms(mol: Any, atom_points: list[tuple[float, float]]) -> list[float]:
    """Measure drawn bond lengths from atom draw coordinates.

    Reading the SVG paths instead would be wrong: RDKit splits wedge and hashed
    stereo bonds into many short sub-paths, which drags the median far below the
    real bond length.
    """
    lengths: list[float] = []
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin >= len(atom_points) or end >= len(atom_points):
            continue
        (x1, y1), (x2, y2) = atom_points[begin], atom_points[end]
        lengths.append(((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5)
    return lengths


def _label_boxes_from_svg(svg: str) -> list[tuple[float, float, float, float]]:
    """Atom labels are emitted as glyph outlines, so their bbox is exact."""
    boxes: list[tuple[float, float, float, float]] = []
    for path_data in _ATOM_PATH.findall(svg):
        numbers = [float(value) for value in _NUMBER.findall(path_data)]
        xs = numbers[0::2]
        ys = numbers[1::2]
        if xs and ys:
            boxes.append((min(xs), min(ys), max(xs), max(ys)))
    return boxes


def _strip_to_body(svg: str) -> str:
    body = re.sub(r"^.*?<!-- END OF HEADER -->", "", svg, flags=re.S)
    body = re.sub(r"<\?xml[^>]*>\s*", "", body).strip()
    body = re.sub(r"<svg[^>]*>", "", body, count=1).strip()
    body = re.sub(r"<rect[^>]*>\s*</rect>", "", body)
    body = re.sub(r"</svg>\s*$", "", body).strip()
    return body


def prepare_mol(smiles: str | None = None, molfile: str | None = None, *, molecule: Any = None) -> Any:
    """Build a depiction-ready mol with 2D coordinates, or raise ValueError."""
    if not RDKIT_AVAILABLE:
        raise ValueError("RDKit is unavailable; molecule depiction requires RDKit.")
    mol = Chem.Mol(molecule) if molecule is not None else None
    if mol is not None:
        pass  # Keep mapped explicit H and any supplied 2D conformer.
    elif molfile:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False)
    elif smiles:
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Could not build a valid molecule for depiction.")
    if mol.GetNumAtoms() == 0:
        raise ValueError("Molecule has no atoms to depict.")
    has_2d = mol.GetNumConformers() > 0 and not mol.GetConformer().Is3D()
    if not has_2d:
        # rdCoordGen gives the cleaner, more journal-like layouts. Use it
        # directly rather than via SetPreferCoordGen, which flips process-global
        # RDKit state and silently changes every other module's depictions.
        try:
            from rdkit.Chem import rdCoordGen

            rdCoordGen.AddCoords(mol)
        except Exception:
            rdDepictor.Compute2DCoords(mol)
    Chem.rdDepictor.StraightenDepiction(mol)
    _normalize_bond_units(mol)
    return mol


def _normalize_bond_units(mol: Any) -> None:
    """Rescale the conformer so the median bond is exactly 1.0 mol unit.

    CoordGen embeds at 1.0 units per bond and Compute2DCoords at 1.5; without
    this, the same preset would draw different absolute bond lengths depending
    on which embedder supplied the coordinates.
    """
    if mol.GetNumBonds() == 0 or mol.GetNumConformers() == 0:
        return
    conf = mol.GetConformer()
    lengths = []
    for bond in mol.GetBonds():
        p = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        q = conf.GetAtomPosition(bond.GetEndAtomIdx())
        lengths.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
    lengths.sort()
    median = lengths[len(lengths) // 2]
    if median <= 1e-6 or abs(median - 1.0) < 1e-6:
        return
    for idx in range(mol.GetNumAtoms()):
        point = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, (point.x / median, point.y / median, 0.0))


def molecule_extent(mol: Any) -> tuple[float, float]:
    """Molecule width/height in RDKit mol units."""
    conf = mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    for idx in range(mol.GetNumAtoms()):
        point = conf.GetAtomPosition(idx)
        xs.append(point.x)
        ys.append(point.y)
    return max(max(xs) - min(xs), 0.0), max(max(ys) - min(ys), 0.0)


def canvas_for(mol: Any, preset: dict[str, Any]) -> tuple[int, int]:
    """Size a canvas so the molecule fits at the preset's exact bond length.

    Undersizing is the one thing that defeats ``fixedBondLength``: RDKit falls
    back to shrink-to-fit and the figure loses its uniform geometry.
    """
    bond_px = float(preset["bond_length_px"])
    pad = float(preset["cell_padding"])
    px_per_unit = bond_px / MOL_UNITS_PER_BOND
    extent_w, extent_h = molecule_extent(mol)
    # Atom labels extend past the outermost atom centres, so reserve a label
    # allowance on every side in addition to the requested padding.
    label_allowance = bond_px * 1.1
    width = extent_w * px_per_unit + 2 * pad + label_allowance
    height = extent_h * px_per_unit + 2 * pad + label_allowance
    return max(int(round(width)), 56), max(int(round(height)), 56)


def render_molecule(
    smiles: str | None = None,
    molfile: str | None = None,
    preset: str | dict[str, Any] | None = None,
    mol: Any = None,
    highlight_atoms: list[int] | None = None,
) -> MoleculeDepiction:
    """Draw one molecule at the preset's fixed geometry.

    The canvas estimate cannot know how far atom labels will extend, so the
    depiction is measured and the canvas grown until the drawn bond length
    matches the preset. Without this, tall or label-heavy molecules quietly come
    back shrunk and the figure loses its uniform geometry.

    Raises ``ValueError`` when the structure cannot be depicted, so callers fail
    loudly instead of emitting a blank cell.
    """
    resolved = preset if isinstance(preset, dict) else resolve_preset(preset)
    if mol is None:
        mol = prepare_mol(smiles=smiles, molfile=molfile)
    width, height = canvas_for(mol, resolved)
    target = float(resolved["bond_length_px"])

    depiction = _draw_at_size(mol, resolved, width, height, highlight_atoms)
    for _ in range(4):
        drawn = max(depiction.bond_lengths_px) if depiction.bond_lengths_px else target
        if drawn >= target * 0.99:
            break
        # Grow by the exact measured shortfall plus a small margin.
        growth = min(target / drawn, 2.0) * 1.02
        width = int(round(width * growth))
        height = int(round(height * growth))
        depiction = _draw_at_size(mol, resolved, width, height, highlight_atoms)

    if depiction.bond_lengths_px:
        drawn = max(depiction.bond_lengths_px)
        if drawn < target * 0.9:
            depiction.warnings.append(
                f"Depiction was scaled below the target bond length "
                f"({drawn:.1f}px drawn vs {target:.1f}px requested)."
            )
    return _trim_to_ink(depiction, float(resolved["cell_padding"]))


def _trim_to_ink(depiction: MoleculeDepiction, padding: float) -> MoleculeDepiction:
    """Crop a depiction to its real ink bounds plus ``padding``.

    The grow loop deliberately overshoots the canvas, so trimming afterwards is
    what lets a composed figure pack cells tightly without clipping. Geometry is
    shifted with the artwork so audit coordinates stay valid.
    """
    if depiction.ink_box is None:
        return depiction
    x0, y0, x1, y1 = depiction.ink_box
    dx = max(0.0, x0 - padding)
    dy = max(0.0, y0 - padding)
    width = int(round(min(depiction.width - dx, x1 - x0 + 2 * padding)))
    height = int(round(min(depiction.height - dy, y1 - y0 + 2 * padding)))
    if width < 8 or height < 8:
        return depiction
    body = depiction.svg_body
    if dx or dy:
        body = f"<g transform='translate({-dx:.2f},{-dy:.2f})'>\n{body}\n</g>"
    return MoleculeDepiction(
        svg_body=body,
        width=width,
        height=height,
        bond_lengths_px=depiction.bond_lengths_px,
        atom_points=[(x - dx, y - dy) for x, y in depiction.atom_points],
        label_boxes=[(a - dx, b - dy, c - dx, d - dy) for a, b, c, d in depiction.label_boxes],
        atom_label_boxes={
            idx: [(a - dx, b - dy, c - dx, d - dy) for a, b, c, d in boxes]
            for idx, boxes in depiction.atom_label_boxes.items()
        },
        ink_box=(x0 - dx, y0 - dy, x1 - dx, y1 - dy),
        warnings=depiction.warnings,
    )


def _draw_at_size(
    mol: Any,
    resolved: dict[str, Any],
    width: int,
    height: int,
    highlight_atoms: list[int] | None,
) -> MoleculeDepiction:
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.fixedBondLength = float(resolved["bond_length_px"])
    options.bondLineWidth = float(resolved["bond_line_width"])
    options.baseFontSize = float(resolved["base_font_size"])
    options.minFontSize = int(resolved["min_font_size"])
    options.maxFontSize = int(resolved["max_font_size"])
    options.padding = 0.0
    options.centreMoleculesBeforeDrawing = True
    options.clearBackground = False
    options.scaleBondWidth = False
    options.explicitMethyl = bool(resolved.get("explicit_methyl", False))
    if resolved.get("monochrome", True):
        options.useBWAtomPalette()
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol, highlightAtoms=highlight_atoms or []
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    atom_points: list[tuple[float, float]] = []
    for idx in range(mol.GetNumAtoms()):
        try:
            point = drawer.GetDrawCoords(idx)
            atom_points.append((float(point.x), float(point.y)))
        except Exception:  # pragma: no cover - defensive
            break

    bond_lengths = _bond_lengths_from_atoms(mol, atom_points)
    label_boxes = _label_boxes_from_svg(svg)
    atom_label_boxes: dict[int, list[tuple[float, float, float, float]]] = {}
    for atom_index, path_data in re.findall(r"<path\s+class='atom-(\d+)'\s+d='([^']*)'", svg):
        numbers = [float(value) for value in _NUMBER.findall(path_data)]
        xs, ys = numbers[0::2], numbers[1::2]
        if xs and ys:
            atom_label_boxes.setdefault(int(atom_index), []).append((min(xs), min(ys), max(xs), max(ys)))
    warnings: list[str] = []

    ink_box: tuple[float, float, float, float] | None = None
    xs = [x for x, _ in atom_points] + [box[0] for box in label_boxes] + [box[2] for box in label_boxes]
    ys = [y for _, y in atom_points] + [box[1] for box in label_boxes] + [box[3] for box in label_boxes]
    if xs and ys:
        ink_box = (min(xs), min(ys), max(xs), max(ys))

    target = float(resolved["bond_length_px"])
    if bond_lengths:
        longest = max(bond_lengths)
        # A shortfall means the canvas clamped the molecule, so the caller's
        # figure would silently lose its uniform bond length.
        if longest < target * 0.9:
            warnings.append(
                f"Depiction was scaled below the target bond length "
                f"({longest:.1f}px drawn vs {target:.1f}px requested); the cell is too small."
            )

    return MoleculeDepiction(
        svg_body=_strip_to_body(svg),
        width=width,
        height=height,
        bond_lengths_px=bond_lengths,
        atom_points=atom_points,
        label_boxes=label_boxes,
        atom_label_boxes=atom_label_boxes,
        ink_box=ink_box,
        warnings=warnings,
    )
