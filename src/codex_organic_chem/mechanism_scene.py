"""Compact, object-anchored mechanism composition using the shared RDKit painter.

Molecular typography is never reimplemented here.  This module places prepared
depictions, routes local electron arrows around their actual ink, and measures
the result.  It is a layout aid, not a guarantee of manuscript acceptance.
All distances in the authoring controls are in bond lengths, not canvas pixels.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from html import escape
from itertools import pairwise
from typing import Any

from .figure_audit import FigureElement, audit_figure
from .figure_style import (
    MoleculeDepiction,
    _normalize_bond_units,
    prepare_mol,
    render_molecule,
    resolve_preset,
)
from .mechanism_semantics import parse_molecule, state_selection

Point = tuple[float, float]
Box = tuple[float, float, float, float]

# Preserve the rounded open V's established proportions. Its alignment, rather
# than a narrower included angle or a filled triangle, controls the curved join.
PAIR_HEAD_LENGTH = 4.8
PAIR_HEAD_HALF_WIDTH = 2.5


def add(p: Point, q: Point) -> Point:
    return p[0] + q[0], p[1] + q[1]


def sub(p: Point, q: Point) -> Point:
    return p[0] - q[0], p[1] - q[1]


def mul(p: Point, scale: float) -> Point:
    return p[0] * scale, p[1] * scale


def norm(p: Point) -> float:
    return math.hypot(*p)


def unit(p: Point) -> Point:
    return mul(p, 1 / (norm(p) or 1))


def union(boxes: list[Box]) -> Box:
    if not boxes:
        return 0, 0, 0, 0
    return (
        min(b[0] for b in boxes),
        min(b[1] for b in boxes),
        max(b[2] for b in boxes),
        max(b[3] for b in boxes),
    )


def shift_box(box: Box, p: Point) -> Box:
    return box[0] + p[0], box[1] + p[1], box[2] + p[0], box[3] + p[1]


def point_box(p: Point, radius: float = 0) -> Box:
    return p[0] - radius, p[1] - radius, p[0] + radius, p[1] + radius


def inside(p: Point, box: Box, pad: float = 0) -> bool:
    return box[0] - pad < p[0] < box[2] + pad and box[1] - pad < p[1] < box[3] + pad


def segment_distance(p: Point, a: Point, b: Point) -> float:
    v, w = sub(b, a), sub(p, a)
    t = max(0.0, min(1.0, (v[0] * w[0] + v[1] * w[1]) / (norm(v) ** 2 or 1)))
    return norm(sub(p, add(a, mul(v, t))))


@dataclass
class Molecule:
    index: int
    spec: dict
    mol: Any
    depiction: MoleculeDepiction
    maps: dict[int, int]
    offset: Point = (0, 0)

    def point(self, atom_map: int) -> Point:
        return add(self.depiction.atom_points[self.maps[atom_map]], self.offset)

    @property
    def ink(self) -> Box:
        return shift_box(
            self.depiction.ink_box
            or (0, 0, self.depiction.width, self.depiction.height),
            self.offset,
        )

    def label_boxes(self, atom_map: int | None = None) -> list[Box]:
        boxes = (
            self.depiction.label_boxes
            if atom_map is None
            else self.depiction.atom_label_boxes.get(self.maps[atom_map], [])
        )
        return [shift_box(b, self.offset) for b in boxes]

    def entry(self) -> dict:
        positions = {
            idx: add(p, self.offset) for idx, p in enumerate(self.depiction.atom_points)
        }
        return {
            "mol": self.mol,
            "spec": self.spec,
            "mol_index": self.index - 1,
            "positions": positions,
            "positions_by_map": {m: positions[i] for m, i in self.maps.items()},
            "map_to_idx": self.maps,
            "mol_box": (
                self.ink[0],
                self.ink[1],
                self.ink[2] - self.ink[0],
                self.ink[3] - self.ink[1],
            ),
            "px_scale": 30,
        }


@dataclass
class Anchor:
    molecule: Molecule
    kind: str
    maps: tuple[int, ...]
    fraction: float = 0.5

    @property
    def point(self) -> Point:
        if self.kind == "bond":
            a, b = (self.molecule.point(m) for m in self.maps)
            return add(a, mul(sub(b, a), self.fraction))
        return self.molecule.point(self.maps[0])

    def object_id(self, panel: int) -> str:
        prefix = f"p{panel}:m{self.molecule.index}:"
        if self.kind == "bond":
            return prefix + "bond:" + "-".join(map(str, sorted(self.maps)))
        if self.kind == "lone_pair":
            return prefix + f"lp:{self.maps[0]}:0"
        return prefix + f"atom:{self.maps[0]}"


@dataclass
class Curve:
    source: Anchor
    target: Anchor
    points: tuple[Point, Point, Point, Point]
    samples: list[Point]
    lp: Point | None
    head: list[Point]
    length: float
    hits: int
    lift: float
    electrons: int
    index: int
    cost: float = 0


@dataclass
class Panel:
    molecules: list[Molecule]
    operators: list[tuple[str, Box]]
    curves: list[Curve]
    issues: list[str]
    warnings: list[str]
    spec: dict
    width: float = 0
    height: float = 0
    offset: Point = (0, 0)
    boxes: list[Box] = field(default_factory=list)


def _rotate(mol: Any, degrees: float) -> None:
    angle = math.radians(degrees)
    c, s = math.cos(angle), math.sin(angle)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (c * p.x - s * p.y, s * p.x + c * p.y, 0))


def _face(mol: Any, maps: dict[int, int], atom_map: int, degrees: float) -> None:
    """Face the largest free angular sector at the reactive atom toward its partner."""
    idx = maps[atom_map]
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(idx)
    angles = []
    for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
        q = conf.GetAtomPosition(neighbor.GetIdx())
        angles.append(math.atan2(q.y - p.y, q.x - p.x) % (2 * math.pi))
    if not angles:
        return
    angles.sort()
    sectors = [
        (
            angles[(i + 1) % len(angles)]
            + (2 * math.pi if i == len(angles) - 1 else 0)
            - a,
            a,
        )
        for i, a in enumerate(angles)
    ]
    gap, start = max(sectors)
    direction = start + gap / 2
    _rotate(mol, degrees - math.degrees(direction))


def _match_orientation(
    mol: Any, maps: dict, references: list[tuple[Any, dict]]
) -> None:
    """Rotate without reflection to preserve the orientation of a shared scaffold."""
    if not references:
        return
    ref, ref_maps = max(references, key=lambda item: len(set(maps) & set(item[1])))
    common = sorted(set(maps) & set(ref_maps))
    if len(common) < 2:
        return

    # Preserve the coordinates of an unchanged common graph. A protonation
    # should not swap two branches just because the next SMILES embedded anew.
    def common_bonds(molecule, mapping):
        return {
            (a, b): bond.GetBondTypeAsDouble()
            for ai, a in enumerate(common)
            for b in common[ai + 1 :]
            if (bond := molecule.GetBondBetweenAtoms(mapping[a], mapping[b]))
            is not None
        }

    if common_bonds(mol, maps) == common_bonds(ref, ref_maps) and (
        set(common) == set(maps) or set(common) == set(ref_maps)
    ):
        from rdkit import Chem
        from rdkit.Chem import rdDepictor

        template = Chem.Mol(ref)
        conf = template.GetConformer()
        for i in range(template.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (1.5 * p.x, 1.5 * p.y, 0))
        rdDepictor.GenerateDepictionMatching2DStructure(
            mol, template, [(ref_maps[m], maps[m]) for m in common]
        )
        _normalize_bond_units(mol)
        return
    a = [mol.GetConformer().GetAtomPosition(maps[m]) for m in common]
    b = [ref.GetConformer().GetAtomPosition(ref_maps[m]) for m in common]
    ac = sum(p.x for p in a) / len(a), sum(p.y for p in a) / len(a)
    bc = sum(p.x for p in b) / len(b), sum(p.y for p in b) / len(b)
    cross = sum(
        (p.x - ac[0]) * (q.y - bc[1]) - (p.y - ac[1]) * (q.x - bc[0])
        for p, q in zip(a, b)
    )
    dot = sum(
        (p.x - ac[0]) * (q.x - bc[0]) + (p.y - ac[1]) * (q.y - bc[1])
        for p, q in zip(a, b)
    )
    _rotate(mol, math.degrees(math.atan2(cross, dot)))


def _arrow_ref(arrow: dict, prefix: str) -> tuple[str, tuple[int, ...]]:
    if prefix + "_bond" in arrow:
        return "bond", tuple(map(int, arrow[prefix + "_bond"]))
    if prefix + "_lone_pair_atom_map" in arrow:
        return "lone_pair", (int(arrow[prefix + "_lone_pair_atom_map"]),)
    if prefix + "_atom_map" in arrow:
        return "atom", (int(arrow[prefix + "_atom_map"]),)
    raise ValueError("compact layout requires chemical atom/bond anchors")


def _find_molecule(
    arrow: dict, prefix: str, molecules: list[Molecule], before: list[int]
) -> Molecule:
    _, maps = _arrow_ref(arrow, prefix)
    hint = arrow.get(prefix + "_molecule_index")
    matches = [
        m
        for m in molecules
        if set(maps) <= set(m.maps) and (m.index == hint if hint else m.index in before)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{prefix} anchor is missing or ambiguous in the starting state"
        )
    return matches[0]


def _make_molecules(
    panel: dict, preset: dict, before: list[int], warnings: list[str], prior: dict
) -> list[Molecule]:
    from rdkit import Chem

    prepared = []
    for index, spec in enumerate(panel.get("molecules", []), 1):
        raw = parse_molecule(spec)
        if raw is None or not raw.GetNumAtoms():
            warnings.append(f"Molecule {index} could not be drawn.")
            continue
        mol = prepare_mol(molecule=raw)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        maps = {
            a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()
        }
        prepared.append((index, spec, mol, maps))
    # Work on molecules, then render them. Never rotate glyphs or reflect stereo.
    facing: dict[int, tuple[int, float]] = {}
    for arrow in panel.get("arrows", []):
        try:
            _, src = _arrow_ref(arrow, "from")
            _, dst = _arrow_ref(arrow, "to")
        except ValueError:
            continue
        matches = []
        for prefix, refs in (("from", src), ("to", dst)):
            hint = arrow.get(prefix + "_molecule_index")
            found = [
                p
                for p in prepared
                if set(refs) <= set(p[3]) and (p[0] == hint if hint else p[0] in before)
            ]
            matches.append(found[0] if len(found) == 1 else None)
        a, b = matches
        if a and b and a[0] != b[0]:
            facing.setdefault(a[0], (src[0], 0 if a[0] < b[0] else 180))
            facing.setdefault(b[0], (dst[0], 180 if a[0] < b[0] else 0))
    references = []
    molecules = []
    for index, spec, mol, maps in prepared:
        previous = prior.get(Chem.MolToSmiles(mol))
        if previous is not None:
            _match_orientation(mol, maps, [previous])
        elif index in facing:
            _face(mol, maps, *facing[index])
        elif index not in before:
            _match_orientation(mol, maps, references)
        rotation = spec.get("depiction", {}).get("rotate_deg")
        if rotation is not None:
            _rotate(mol, float(rotation))
        axis = spec.get("depiction", {}).get("bond_axis")
        if axis:
            a, b = (
                mol.GetConformer().GetAtomPosition(maps[int(m)])
                for m in axis["atom_maps"]
            )
            _rotate(
                mol,
                float(axis.get("angle_deg", 0))
                - math.degrees(math.atan2(b.y - a.y, b.x - a.x)),
            )
        if index in before:
            references.append((mol, maps))
        drawing = Chem.Mol(mol)
        for atom in drawing.GetAtoms():
            atom.SetAtomMapNum(0)
        # A single-carbon reagent (e.g. CH3Br) benefits from a terminal methyl
        # label; larger molecules use the same skeletal style as reaction figures.
        single_carbon = sum(a.GetAtomicNum() == 6 for a in drawing.GetAtoms()) == 1
        depiction = render_molecule(
            mol=drawing, preset={**preset, "explicit_methyl": single_carbon}
        )
        warnings.extend(depiction.warnings)
        molecules.append(Molecule(index, spec, mol, depiction, maps))
    return molecules


def _place(
    panel: dict, molecules: list[Molecule], before: list[int], bond: float
) -> list[tuple[str, Box]]:
    gap = float(panel.get("composition", {}).get("species_gap_bonds", 1.35)) * bond
    # The baseline passes through the reacting atom where possible, instead of
    # the centre of an arbitrary molecule-sized cell.
    baselines = {}
    for arrow in panel.get("arrows", []):
        for prefix in ("from", "to"):
            try:
                m = _find_molecule(arrow, prefix, molecules, before)
                _, maps = _arrow_ref(arrow, prefix)
                if len(maps) == 1:
                    baselines.setdefault(m.index, m.point(maps[0])[1])
            except ValueError:
                continue
    operators = []
    cursor = 0.0
    previous = None
    for m in molecules:
        if previous:
            transition = previous.index in before and m.index not in before
            if transition:
                operators.append(
                    ("reaction", (cursor + 20, -3, cursor + 20 + 1.6 * bond, 3))
                )
                cursor += 1.6 * bond + 40
            else:
                operators.append(
                    ("plus", (cursor + gap / 2 - 4, -4, cursor + gap / 2 + 4, 4))
                )
                cursor += gap
        b = m.ink
        baseline = baselines.get(m.index, (b[1] + b[3]) / 2)
        m.offset = (cursor - b[0], -baseline)
        # Optional exact art direction is dimensionless and does not alter maps.
        shift = m.spec.get("depiction", {}).get("offset_bonds", [0, 0])
        m.offset = add(m.offset, (float(shift[0]) * bond, float(shift[1]) * bond))
        cursor = max(cursor, m.ink[2])
        previous = m
    return operators


def _boundary(anchor: Anchor, toward: Point, gap: float = 2.8) -> Point:
    """Stop at the element glyph, keeping implicit H and charge glyphs separate."""
    p = anchor.point
    d = unit(sub(toward, p))
    if anchor.kind == "bond":
        return add(p, mul(d, 2.0))
    boxes = anchor.molecule.label_boxes(anchor.maps[0])
    # Extending a ray across all of "HO" can incorrectly place an oxygen lone
    # pair beside the H. Leave those glyphs as obstacles for the route search.
    core = [b for b in boxes if inside(p, b, 0.8)]
    if not core and boxes:
        core = [
            min(
                boxes,
                key=lambda b: norm(sub(p, ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2))),
            )
        ]
    distance = 0.8
    for box in core:
        # Slab intersection with each glyph; disconnected superscripts need not
        # enlarge the entire atom label into a fictitious rectangular obstacle.
        lo, hi = 0.0, float("inf")
        for axis in (0, 1):
            a, b = box[axis] - 0.5, box[axis + 2] + 0.5
            if abs(d[axis]) < 1e-8:
                if not a <= p[axis] <= b:
                    hi = -1
                    break
            else:
                t1, t2 = (a - p[axis]) / d[axis], (b - p[axis]) / d[axis]
                lo, hi = max(lo, min(t1, t2)), min(hi, max(t1, t2))
        if hi >= lo:
            distance = max(distance, hi)
    return add(p, mul(d, distance + gap))


def _sample(points: tuple[Point, Point, Point, Point], count: int = 64) -> list[Point]:
    return [
        (
            sum(
                v[0] * w
                for v, w in zip(
                    points,
                    ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t * t, t**3),
                )
            ),
            sum(
                v[1] * w
                for v, w in zip(
                    points,
                    ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t * t, t**3),
                )
            ),
        )
        for t in (i / count for i in range(count + 1))
    ]


def _head(end: Point, control: Point, electrons: int) -> list[Point]:
    d = unit(sub(end, control))
    n = (-d[1], d[0])
    if electrons == 1:
        back = add(end, mul(d, -4.8))
        return [add(back, mul(n, 2.5)), end]
    back = add(end, mul(d, -PAIR_HEAD_LENGTH))
    half_width = PAIR_HEAD_HALF_WIDTH
    return [add(back, mul(n, half_width)), end, add(back, mul(n, -half_width))]


def _shaft_points(curve: Curve) -> tuple[Point, Point, Point, Point]:
    """Join the curved body tangentially to the open V's visible angle bisector.

    Matching the derivative at the tip alone is insufficient on a tight arc:
    the shaft can visibly crowd one wing. A short straight finish runs along
    the bisector, and the cubic meets it without a kink. Keep the source, its
    tangent and the chemical target fixed; adjust only the terminal approach.
    """
    if curve.electrons == 1:
        return curve.points
    samples = _sample(curve.points)
    total = sum(norm(sub(a, b)) for a, b in pairwise(samples))
    finish = min(PAIR_HEAD_LENGTH, total * 0.25)
    remaining = finish
    for i in range(len(samples) - 1, 0, -1):
        length = norm(sub(samples[i], samples[i - 1]))
        if length >= remaining:
            t = (i - remaining / length) / (len(samples) - 1)
            row = list(curve.points)
            prefix = [row[0]]
            while len(row) > 1:
                row = [add(mul(a, 1 - t), mul(b, t)) for a, b in pairwise(row)]
                prefix.append(row[0])
            tip = curve.points[-1]
            direction = unit(sub(tip, curve.points[-2]))
            join = sub(tip, mul(direction, finish))
            handle = norm(sub(prefix[-1], prefix[-2]))
            return prefix[0], prefix[1], sub(join, mul(direction, handle)), join
        remaining -= length
    return curve.points


def _obstacles(
    molecules: list[Molecule], operators: list[tuple[str, Box]]
) -> tuple[list, list]:
    boxes, segments = [], []
    for m in molecules:
        boxes.extend(m.label_boxes())
        for b in m.mol.GetBonds():
            a, c = b.GetBeginAtom(), b.GetEndAtom()
            # Multiple bonds need a wider clearance corridor.
            pad = 2.4 if b.GetBondTypeAsDouble() == 1 else 4.4
            pa = add(m.depiction.atom_points[a.GetIdx()], m.offset)
            pc = add(m.depiction.atom_points[c.GetIdx()], m.offset)
            segments.append(
                (
                    m.index,
                    tuple(sorted((a.GetAtomMapNum(), c.GetAtomMapNum()))),
                    pa,
                    pc,
                    pad,
                )
            )
    boxes.extend(b for _, b in operators)
    return boxes, segments


def _curve_hits(curve: Curve, boxes: list, segments: list) -> int:
    hits = 0
    points = [*curve.samples, *curve.head]
    if curve.lp:
        points.append(curve.lp)
    for p in points:
        if any(inside(p, b, 1.2) for b in boxes):
            hits += 1
        for mi, maps, a, b, pad in segments:
            # Contact is legitimate only immediately at the chemically anchored
            # endpoint, not throughout the first/last 20% of the arrow.
            if (
                curve.source.kind == "bond"
                and mi == curve.source.molecule.index
                and maps == tuple(sorted(curve.source.maps))
                and norm(sub(p, curve.source.point)) < 5.0
            ):
                continue
            if (
                mi == curve.target.molecule.index
                and set(curve.target.maps) <= set(maps)
                and norm(sub(p, curve.target.point)) < 5.2
            ):
                continue
            if segment_distance(p, a, b) < pad:
                hits += 1
    return hits


def _candidates(
    arrow: dict,
    index: int,
    source: Anchor,
    target: Anchor,
    bond: float,
    obstacles: tuple,
) -> list[Curve]:
    a, d = source.point, target.point
    direction = unit(sub(d, a))
    if norm(sub(d, a)) < 1:
        raise ValueError("arrow endpoints coincide")
    n = (-direction[1], direction[0])
    electrons = int(
        arrow.get("electron_count", 1 if arrow.get("kind") == "single_electron" else 2)
    )
    # Old pixel curvatures affect side only. New authoring controls are in bonds.
    preferred_side = -1 if float(arrow.get("curvature", -1)) < 0 else 1
    controls = arrow.get("routing", {})
    if controls.get("side") in ("left", "right"):
        preferred_side = -1 if controls["side"] == "left" else 1
    sides = (
        [preferred_side]
        if controls.get("lock_side")
        else [preferred_side, -preferred_side]
    )
    lifts = (
        [float(controls["lift_bonds"]) * bond]
        if "lift_bonds" in controls
        else [bond * k for k in (0.55, 0.8, 1.05, 1.35)]
    )
    result = []
    for side in sides:
        for lift in lifts:
            c1 = add(add(a, mul(sub(d, a), 0.20)), mul(n, lift * side))
            c2 = add(add(a, mul(sub(d, a), 0.80)), mul(n, lift * side))
            lp = _boundary(source, c1, 3.0) if source.kind == "lone_pair" else None
            start = (
                add(lp, mul(unit(sub(c1, lp)), 3.0))
                if lp
                else _boundary(source, c1, 1.5)
            )
            end = _boundary(target, c2)
            points = start, c1, c2, end
            samples = _sample(points)
            length = sum(norm(sub(p, q)) for p, q in pairwise(samples))
            curve = Curve(
                source,
                target,
                points,
                samples,
                lp,
                _head(end, c2, electrons),
                length,
                0,
                lift * side,
                electrons,
                index,
            )
            if electrons == 2:
                shaft = _shaft_points(curve)
                curve.samples = _sample(shaft)
                curve.samples.extend(
                    add(mul(shaft[-1], 1 - t), mul(end, t))
                    for t in (i / 8 for i in range(1, 9))
                )
                curve.length = sum(norm(sub(p, q)) for p, q in pairwise(curve.samples))
            curve.hits = _curve_hits(curve, *obstacles)
            curve.cost = (
                curve.hits * 1000
                # Rank the established routes by their original arc length;
                # audit the actual refined shaft and head above.
                + length / bond
                + abs(lift / bond) * 0.15
                + (0.08 if side != preferred_side else 0)
            )
            result.append(curve)
    return sorted(result, key=lambda c: c.cost)


def _arrow_crossings(a: Curve, b: Curve) -> int:
    # Heads and lone pairs count as ink; overlapping tails at a homolysed bond
    # are the sole exception, confined to a 5px neighbourhood of that bond.
    hits = 0
    common_source = (
        a.source.kind == b.source.kind == "bond"
        and a.source.molecule.index == b.source.molecule.index
        and set(a.source.maps) == set(b.source.maps)
    )
    for p in [*a.samples[::2], *a.head, *([a.lp] if a.lp else [])]:
        if common_source and norm(sub(p, a.source.point)) < 5:
            continue
        if any(
            segment_distance(p, x, y) < 2.5
            for x, y in zip(b.samples[::2], b.samples[2::2])
        ):
            hits += 1
    return hits


def _route(
    panel: dict,
    molecules: list[Molecule],
    before: list[int],
    operators: list,
    bond: float,
) -> tuple[list[Curve], list[str]]:
    obstacles = _obstacles(molecules, operators)
    issues = []
    beam: list[tuple[float, list[Curve]]] = [(0, [])]
    for index, arrow in enumerate(panel.get("arrows", []), 1):
        try:
            source_mol = _find_molecule(arrow, "from", molecules, before)
            target_mol = _find_molecule(arrow, "to", molecules, before)
            source = Anchor(source_mol, *_arrow_ref(arrow, "from"))
            target = Anchor(target_mol, *_arrow_ref(arrow, "to"))
            options = _candidates(arrow, index, source, target, bond, obstacles)
        except (ValueError, KeyError) as error:
            issues.append(f"Arrow {index}: {error}.")
            continue
        possibilities = []
        for cost, previous in beam:
            for candidate in options:
                crossings = sum(
                    _arrow_crossings(candidate, other)
                    + _arrow_crossings(other, candidate)
                    for other in previous
                )
                possibilities.append(
                    (cost + candidate.cost + crossings * 1000, [*previous, candidate])
                )
        beam = sorted(possibilities, key=lambda item: item[0])[:12]
    chosen = beam[0][1]
    for curve in chosen:
        if curve.hits:
            issues.append(
                f"Arrow {curve.index}: electron arrow intersects molecule or operator ink ({curve.hits} sampled contacts)."
            )
        if curve.length / bond > 3.5:
            issues.append(
                f"Arrow {curve.index}: electron arrow is {curve.length / bond:.1f} bond lengths long; move the reacting sites closer."
            )
    for i, curve in enumerate(chosen):
        for other in chosen[i + 1 :]:
            if _arrow_crossings(curve, other) or _arrow_crossings(other, curve):
                issues.append(
                    f"Arrows {curve.index} and {other.index} overlap; adjust the molecule orientation or routing side."
                )
    return chosen, issues


def _compact_reactive_pair(
    panel: dict,
    molecules: list[Molecule],
    before: list[int],
    operators: list,
    bond: float,
) -> list:
    """Move a small participating reagent next to the actual reactive site.

    A large substrate's bounding box is not a reason for a long electron arrow.
    Search a small set of positions around its reactive atom, then reflow the
    remaining species. The graph and the scaffold orientation stay unchanged.
    This bounded pass handles one intermolecular pair; more complex compositions
    keep explicit art direction and are still checked by the collision audit.
    """
    if not panel.get("composition", {}).get("auto_compact", True):
        return operators
    pairs = []
    for i, arrow in enumerate(panel.get("arrows", []), 1):
        try:
            a = _find_molecule(arrow, "from", molecules, before)
            b = _find_molecule(arrow, "to", molecules, before)
            if a != b:
                pairs.append((i, arrow, a, b))
        except ValueError:
            continue
    if len(pairs) != 1:
        return operators
    index, arrow, a, b = pairs[0]
    if any(m.spec.get("depiction", {}).get("offset_bonds") for m in (a, b)):
        return operators
    source, target = (
        Anchor(a, *_arrow_ref(arrow, "from")),
        Anchor(b, *_arrow_ref(arrow, "to")),
    )
    original = _candidates(
        arrow, index, source, target, bond, _obstacles(molecules, operators)
    )[0]
    if original.length < 2.35 * bond and original.hits == 0:
        return operators
    moving, _fixed = (
        (a, b) if a.mol.GetNumHeavyAtoms() < b.mol.GetNumHeavyAtoms() else (b, a)
    )
    moving_anchor, fixed_anchor = (source, target) if moving is a else (target, source)
    original_offset = moving.offset
    own = sub(moving_anchor.point, moving.offset)
    pivot = fixed_anchor.point
    candidates = []
    for angle in range(-180, 180, 30):
        direction = (math.cos(math.radians(angle)), math.sin(math.radians(angle)))
        for radius in (1.6, 1.9, 2.2):
            moving.offset = sub(add(pivot, mul(direction, radius * bond)), own)
            if any(
                min(moving.ink[2], m.ink[2]) - max(moving.ink[0], m.ink[0]) > -7
                and min(moving.ink[3], m.ink[3]) - max(moving.ink[1], m.ink[1]) > -7
                for m in molecules
                if m is not moving and m.index in before
            ):
                continue
            options = _candidates(
                arrow, index, source, target, bond, _obstacles([a, b], [])
            )
            candidate = options[0]
            bound = union([a.ink, b.ink])
            cost = (
                candidate.cost
                + (bound[2] - bound[0] + bound[3] - bound[1]) / bond * 0.12
            )
            candidates.append((cost, moving.offset))
    moving.offset = original_offset
    if not candidates:
        return operators
    # Route the entire step for the best few positions, so shortening one arrow
    # cannot force another electron arrow through it.
    evaluated = []
    for cost, offset in sorted(candidates)[:5]:
        moving.offset = offset
        trial, issues = _route(panel, [a, b], [a.index, b.index], [], bond)
        evaluated.append(
            (len(issues) * 10000 + sum(c.cost for c in trial) + cost, offset)
        )
    moving.offset = min(evaluated)[1]
    # This is a spatial encounter complex: omit a redundant '+' inside the
    # interacting pair. Keep '+' between spectators and between product species.
    left = min(a.ink[0], b.ink[0])
    a.offset, b.offset = add(a.offset, (-left, 0)), add(b.offset, (-left, 0))
    cursor = max(a.ink[2], b.ink[2])
    operators = []
    gap = 1.35 * bond
    previous_before = True
    for m in [m for m in molecules if m not in (a, b)]:
        if previous_before and m.index not in before:
            operators.append(
                ("reaction", (cursor + 20, -3, cursor + 20 + 1.6 * bond, 3))
            )
            cursor += 1.6 * bond + 40
        else:
            operators.append(
                ("plus", (cursor + gap / 2 - 4, -4, cursor + gap / 2 + 4, 4))
            )
            cursor += gap
        m.offset = add(m.offset, (cursor - m.ink[0], 0))
        cursor = m.ink[2]
        previous_before = m.index in before
    return operators


def _panel(panel: dict, preset: dict, prior: dict | None = None) -> Panel:
    warnings: list[str] = []
    before = state_selection(panel, "starting_state") or list(
        range(1, len(panel.get("molecules", [])) + 1)
    )
    molecules = _make_molecules(panel, preset, before, warnings, prior or {})
    operators = _place(panel, molecules, before, float(preset["bond_length_px"]))
    operators = _compact_reactive_pair(
        panel, molecules, before, operators, float(preset["bond_length_px"])
    )
    curves, issues = _route(
        panel, molecules, before, operators, float(preset["bond_length_px"])
    )
    boxes = [m.ink for m in molecules] + [b for _, b in operators]
    for curve in curves:
        boxes.extend(
            point_box(p, 2.5)
            for p in [*curve.samples, *curve.head, *([curve.lp] if curve.lp else [])]
        )
    return Panel(molecules, operators, curves, issues, warnings, panel, boxes=boxes)


def _points(points: list[Point]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _panel_svg(panel: Panel, index: int, preset: dict, display: dict) -> str:
    parts = [
        f"<g id='panel-{index}' transform='translate({panel.offset[0]:.2f},{panel.offset[1]:.2f})'>"
    ]
    for molecule in panel.molecules:
        dx, dy = molecule.offset
        parts.append(
            f"<g id='p{index}:m{molecule.index}' class='mechanism-molecule' data-painter='rdkit-shared' transform='translate({dx:.2f},{dy:.2f})'>{molecule.depiction.svg_body}</g>"
        )
        if display.get("show_atom_maps"):
            for atom_map in molecule.maps:
                x, y = molecule.point(atom_map)
                parts.append(
                    f"<text class='atom-map' x='{x + 5:.2f}' y='{y + 16:.2f}'>{atom_map}</text>"
                )
    for kind, (x0, y0, x1, y1) in panel.operators:
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        if kind == "plus":
            parts.append(
                f"<path class='scheme-operator' d='M {x0:.2f} {cy:.2f} H {x1:.2f} M {cx:.2f} {y0:.2f} V {y1:.2f}'/>"
            )
        else:
            parts.append(
                f"<path class='reaction-arrow' d='M {x0:.2f} {cy:.2f} H {x1 - 3.6:.2f}'/>"
            )
            parts.append(
                f"<polygon class='reaction-head' points='{_points([(x1 - 6, cy - 2.7), (x1, cy), (x1 - 6, cy + 2.7)])}'/>"
            )
    for curve in panel.curves:
        start, c1, c2, end = _shaft_points(curve)
        if curve.lp:
            direction = unit(sub(start, curve.lp))
            perp = (-direction[1], direction[0])
            dots = [add(curve.lp, mul(perp, sign * 1.8)) for sign in (-1, 1)]
            parts.append(
                f"<g class='lone-pair-object' id='{curve.source.object_id(index)}'>"
            )
            parts.extend(
                f"<circle cx='{x:.2f}' cy='{y:.2f}' r='0.9'/>" for x, y in dots
            )
            parts.append("</g>")
        path = f"M {start[0]:.2f} {start[1]:.2f} C {c1[0]:.2f} {c1[1]:.2f} {c2[0]:.2f} {c2[1]:.2f} {end[0]:.2f} {end[1]:.2f}"
        if curve.electrons == 2:
            tip = curve.points[-1]
            path += f" L {tip[0]:.2f} {tip[1]:.2f}"
        parts.append(
            f"<path class='mech-arrow' d='{path}' data-source-object='{curve.source.object_id(index)}' data-target-object='{curve.target.object_id(index)}' data-electrons='{curve.electrons}' data-length-bonds='{curve.length / float(preset['bond_length_px']):.3f}' data-overlap-score='{curve.hits}' data-routed-curvature='{curve.lift:.2f}'/>"
        )
        head_class = "arrow-radical" if curve.electrons == 1 else "arrow-pair"
        parts.append(
            f"<polyline class='electron-head {head_class}' points='{_points(curve.head)}'/>"
        )
    parts.append("</g>")
    return "\n".join(parts)


def compose_mechanism(spec: dict, display: dict) -> dict:
    layout = spec.get("layout", {})
    preset = resolve_preset(
        spec.get("journal_style") or layout.get("style_preset") or "acs"
    )
    from rdkit import Chem

    panels = []
    prior = {}
    for p in spec.get("panels", []):
        panel = _panel(p, preset, prior)
        panels.append(panel)
        end = state_selection(p, "expected_state") or []
        prior = {
            Chem.MolToSmiles(m.mol): (m.mol, m.maps)
            for m in panel.molecules
            if m.index in end
        }
    margin, row_gap = 18.0, 30.0
    y = margin
    width = 120.0
    parts = []
    if display.get("show_title"):
        title = str(spec.get("title", "Mechanism"))
        parts.append(
            f"<text class='figure-title' x='{margin}' y='{y + 12}'>{escape(title)}</text>"
        )
        width = max(width, 7.4 * len(title) + 2 * margin)
        y += 36
    for index, panel in enumerate(panels, 1):
        x0, y0, x1, y1 = union(panel.boxes)
        if display.get("show_panel_titles"):
            title = str(panel.spec.get("title", f"Step {index}"))
            parts.append(
                f"<text class='panel-title' x='{margin}' y='{y + 11:.2f}'>{escape(title)}</text>"
            )
            width = max(width, 6.6 * len(title) + 2 * margin)
            y += 32
        panel.offset = (margin - x0, y - y0)
        panel.width, panel.height = x1 - x0, y1 - y0
        width = max(width, panel.width + 2 * margin)
        parts.append(_panel_svg(panel, index, preset, display))
        y += panel.height + row_gap
    height = max(72.0, y - row_gap + margin)
    elements = []
    issues = []
    warnings = []
    arrow_metrics = []
    entries = []
    for index, panel in enumerate(panels, 1):
        issues.extend(f"Panel {index}: {issue}" for issue in panel.issues)
        warnings.extend(f"Panel {index}: {warning}" for warning in panel.warnings)
        panel_entries = []
        for molecule in panel.molecules:
            molecule.offset = add(molecule.offset, panel.offset)
            elements.append(
                FigureElement(
                    "molecule",
                    f"panel {index} molecule {molecule.index}",
                    molecule.ink,
                    molecule.depiction.bond_lengths_px,
                    molecule.label_boxes(),
                )
            )
            panel_entries.append(molecule.entry())
        entries.append(panel_entries)
        for kind, box in panel.operators:
            elements.append(
                FigureElement(
                    kind, f"panel {index} {kind}", shift_box(box, panel.offset)
                )
            )
        for curve in panel.curves:
            arrow_metrics.append(
                {
                    "panel": index,
                    "arrow": curve.index,
                    "length_bonds": round(
                        curve.length / float(preset["bond_length_px"]), 3
                    ),
                    "ink_contacts": curve.hits,
                    "source": curve.source.object_id(index),
                    "target": curve.target.object_id(index),
                }
            )
    audit = audit_figure(
        renderer="rdkit_compact_mechanism",
        width=width,
        height=height,
        elements=elements,
        target_bond_px=float(preset["bond_length_px"]),
    )
    audit["blocking_issues"].extend(issues)
    audit["warnings"].extend(warnings)
    audit["electron_arrows"] = arrow_metrics
    audit["ready_for_visual_review"] = not audit["blocking_issues"]
    audit["visual_review_required"] = True
    audit["limits"] = [
        "Collision sampling and label bounds are conservative, not a human legibility judgement.",
        "Orientation in a 2D scheme does not demonstrate stereoelectronic alignment or a transition-state geometry.",
    ]
    svg = (
        f"<svg xmlns='http://www.w3.org/2000/svg' xmlns:rdkit='http://www.rdkit.org/xml' width='{math.ceil(width)}' height='{math.ceil(height)}' viewBox='0 0 {width:.2f} {height:.2f}'>\n"
        "<style>.mech-arrow,.electron-head{fill:none;stroke:#111;stroke-width:1.15;stroke-linecap:round;stroke-linejoin:round}"
        ".reaction-head{fill:#111;stroke:none}"
        ".scheme-operator,.reaction-arrow{fill:none;stroke:#111;stroke-width:1.2}"
        ".figure-title,.panel-title,.atom-map{font-family:Arial,Helvetica,sans-serif;fill:#111;font-weight:400}"
        ".figure-title{font-size:14px}.panel-title{font-size:12px}.atom-map{font-size:9px;fill:#666}</style>\n"
        f"<rect width='{width:.2f}' height='{height:.2f}' fill='white'/>\n"
        + "\n".join(parts)
        + "\n</svg>"
    )
    return {
        "svg": svg,
        "figure_audit": audit,
        "entries": entries,
        "width": width,
        "height": height,
        "preset": preset["name"],
    }
