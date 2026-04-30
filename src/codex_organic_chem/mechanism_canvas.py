from __future__ import annotations

import json
import math
from html import escape
from pathlib import Path
from typing import Any

from .rdkit_tools import RDKIT_AVAILABLE

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through service tests
    from rdkit import Chem
    from rdkit.Chem import AllChem
else:  # pragma: no cover
    Chem = None
    AllChem = None


ATOM_COLORS = {
    "O": "#d62728",
    "N": "#1f55ff",
    "S": "#c7a500",
    "P": "#d97706",
    "F": "#149954",
    "Cl": "#149954",
    "Br": "#8c4a18",
    "I": "#6f3c97",
}


JOURNAL_STYLES: dict[str, dict[str, Any]] = {
    "default": {
        "font_family": "Arial,Helvetica,sans-serif",
        "bond_width": 2.2,
        "bond_light_width": 1.4,
        "arrow_width": 2.2,
        "atom_font_size": 16,
        "caption_font_size": 13,
        "max_bond_px": 48,
        "double_bond_spacing": 3.5,
        "arrow_color": "#111",
    },
    "acs": {
        "font_family": "Arial,Helvetica,sans-serif",
        "bond_width": 2.0,
        "bond_light_width": 1.2,
        "arrow_width": 1.8,
        "atom_font_size": 15,
        "caption_font_size": 12,
        "max_bond_px": 44,
        "double_bond_spacing": 3.0,
        "arrow_color": "#111",
    },
    "rsc": {
        "font_family": "Arial,Helvetica,sans-serif",
        "bond_width": 1.8,
        "bond_light_width": 1.1,
        "arrow_width": 1.7,
        "atom_font_size": 14,
        "caption_font_size": 12,
        "max_bond_px": 42,
        "double_bond_spacing": 3.0,
        "arrow_color": "#111",
    },
}


def _style_for_spec(spec: dict[str, Any]) -> dict[str, Any]:
    layout = spec.get("layout", {})
    style_name = str(
        spec.get("journal_style")
        or layout.get("journal_style")
        or layout.get("style_preset")
        or layout.get("style")
        or "default"
    ).lower()
    style = {**JOURNAL_STYLES["default"], **JOURNAL_STYLES.get(style_name, {})}
    for key in (
        "font_family",
        "bond_width",
        "bond_light_width",
        "arrow_width",
        "atom_font_size",
        "caption_font_size",
        "max_bond_px",
        "double_bond_spacing",
        "arrow_color",
        "atom_colors",
    ):
        if key in layout:
            style[key] = layout[key]
    if _is_publication_presentation(spec) and "atom_colors" not in layout:
        style["atom_colors"] = {}
    style["name"] = style_name if style_name in JOURNAL_STYLES else "default"
    return style


def _is_publication_presentation(spec: dict[str, Any]) -> bool:
    layout = spec.get("layout", {})
    mode = str(
        spec.get("presentation")
        or spec.get("presentation_mode")
        or layout.get("presentation")
        or layout.get("presentation_mode")
        or ""
    ).lower()
    return mode in {"publication", "manuscript", "journal"}


def _display_options(spec: dict[str, Any]) -> dict[str, Any]:
    layout = spec.get("layout", {})
    publication = _is_publication_presentation(spec)

    def flag(name: str, default: bool) -> bool:
        return bool(layout.get(name, default))

    lone_pair_display = str(layout.get("lone_pair_display", "source_only" if publication else "all")).lower()
    return {
        "publication": publication,
        "show_title": flag("show_title", not publication),
        "show_subtitle": flag("show_subtitle", not publication),
        "show_panel_titles": flag("show_panel_titles", not publication),
        "show_atom_maps": flag("show_atom_maps", not publication),
        "show_molecule_labels": flag("show_molecule_labels", not publication),
        "show_arrow_labels": flag("show_arrow_labels", not publication),
        "show_notes": flag("show_notes", not publication),
        "show_footer": flag("show_footer", not publication),
        "show_partial_charges": flag("show_partial_charges", not publication),
        "show_formal_charges": flag("show_formal_charges", True),
        "show_chemdoodle_toolbar": flag("show_chemdoodle_toolbar", not publication),
        "lone_pair_display": lone_pair_display,
    }


def _canvas_metrics(spec: dict[str, Any]) -> dict[str, int]:
    layout = spec.get("layout", {})
    display = _display_options(spec)
    panels = spec.get("panels", [])
    columns = max(1, int(layout.get("columns", min(3, max(1, len(panels))))))
    panel_w = int(layout.get("panel_width", 560))
    panel_h = int(layout.get("panel_height", 360))
    rows = max(1, math.ceil(max(1, len(panels)) / columns))
    default_top = 82 if display["show_title"] or display["show_subtitle"] else 12
    default_bottom = 44 if display["show_footer"] else 12
    top_margin = int(layout.get("top_margin", default_top))
    bottom_margin = int(layout.get("bottom_margin", default_bottom))
    return {
        "columns": columns,
        "panel_width": panel_w,
        "panel_height": panel_h,
        "rows": rows,
        "top_margin": top_margin,
        "bottom_margin": bottom_margin,
        "width": columns * panel_w,
        "height": top_margin + rows * panel_h + bottom_margin,
    }


def _svg_header(width: int, height: int, style: dict[str, Any] | None = None) -> list[str]:
    style = style or JOURNAL_STYLES["default"]
    font = style["font_family"]
    bond_width = float(style["bond_width"])
    bond_light_width = float(style["bond_light_width"])
    arrow_width = float(style["arrow_width"])
    atom_font_size = int(style["atom_font_size"])
    caption_font_size = int(style["caption_font_size"])
    arrow_color = str(style["arrow_color"])
    return [
        "<?xml version='1.0' encoding='UTF-8'?>",
        (
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}px' height='{height}px' "
            f"viewBox='0 0 {width} {height}'>"
        ),
        "<defs>",
        (
            "<marker id='arrow-pair' markerWidth='7' markerHeight='5' refX='6.4' refY='2.5' orient='auto'>"
            f"<path d='M0,0 L7,2.5 L0,5 Z' fill='{arrow_color}'/></marker>"
        ),
        (
            "<marker id='arrow-radical' markerWidth='6' markerHeight='5' refX='5.4' refY='2.5' orient='auto'>"
            f"<path d='M0,0 L6,2.5 L0,2.5 Z' fill='{arrow_color}'/></marker>"
        ),
        (
            "<marker id='arrow-step' markerWidth='12' markerHeight='8' refX='11' refY='4' orient='auto'>"
            "<path d='M0,0 L12,4 L0,8 Z' fill='#111'/></marker>"
        ),
        (
            "<style>"
            f".bond{{stroke:#111;stroke-width:{bond_width};stroke-linecap:round}}"
            f".bond-light{{stroke:#111;stroke-width:{bond_light_width};stroke-linecap:round}}"
            f".atom{{font:{atom_font_size}px {font};font-weight:600}}"
            f".map{{font:10px {font};fill:#777}}"
            f".title{{font:20px {font};fill:#111}}.subtitle{{font:14px {font};fill:#333}}"
            f".panel-title{{font:16px {font};fill:#111}}.note{{font:{caption_font_size}px {font};fill:#333}}"
            f".atom-hydrogen{{font:{atom_font_size}px {font};font-weight:600}}"
            f".atom,.atom-hydrogen{{stroke:#fff;stroke-width:3px;paint-order:stroke fill;stroke-linejoin:round}}"
            f".atom-subscript{{font:{max(9, int(atom_font_size * 0.7))}px {font};font-weight:600}}"
            f".charge{{font:15px {font};font-weight:700}}.partial{{font:14px {font};font-weight:700}}"
            f".charge-ring{{fill:#fff;stroke-width:1.2}}.charge-sign{{font:12px {font};font-weight:700;text-anchor:middle;dominant-baseline:central}}"
            f".lp{{fill:#111}}.mech-arrow{{stroke:{arrow_color};stroke-width:{arrow_width};fill:none;stroke-linecap:round;stroke-linejoin:round}}"
            f".callout{{font:{caption_font_size}px {font};fill:#9a3412}}"
            "</style>"
        ),
        "</defs>",
        "<rect width='100%' height='100%' fill='white'/>",
    ]


def _bond_order(bond: Any) -> int:
    if bond.GetIsAromatic():
        return 1
    raw = bond.GetBondTypeAsDouble()
    if raw >= 2.5:
        return 3
    if raw >= 1.5:
        return 2
    return 1


def _unit_normal(x1: float, y1: float, x2: float, y2: float) -> tuple[float, float]:
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy) or 1.0
    return -dy / length, dx / length


def _transform(
    coords: dict[int, tuple[float, float]],
    box: tuple[float, float, float, float],
    max_bond_px: float = 48,
) -> dict[int, tuple[float, float]]:
    x, y, width, height = box
    xs = [pt[0] for pt in coords.values()]
    ys = [pt[1] for pt in coords.values()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    mol_w = max(max_x - min_x, 1.0)
    mol_h = max(max_y - min_y, 1.0)
    scale = min((width - 32) / mol_w, (height - 42) / mol_h, max_bond_px)
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    return {
        idx: (x + width / 2 + (px - center_x) * scale, y + height / 2 - (py - center_y) * scale)
        for idx, (px, py) in coords.items()
    }


def _apply_alignment(
    positions: dict[int, tuple[float, float]],
    map_to_idx: dict[int, int],
    box: tuple[float, float, float, float],
    align_targets: dict[int, tuple[float, float]] | None,
) -> dict[int, tuple[float, float]]:
    if not align_targets:
        return positions
    x, y, width, height = box
    shifts: list[tuple[float, float]] = []
    for atom_map, atom_idx in map_to_idx.items():
        if atom_map not in align_targets or atom_idx not in positions:
            continue
        target_x = x + align_targets[atom_map][0] * width
        target_y = y + align_targets[atom_map][1] * height
        px, py = positions[atom_idx]
        shifts.append((target_x - px, target_y - py))
    if not shifts:
        return positions
    dx = sum(item[0] for item in shifts) / len(shifts)
    dy = sum(item[1] for item in shifts) / len(shifts)
    dx = max(min(dx, width * 0.18), -width * 0.18)
    dy = max(min(dy, height * 0.18), -height * 0.18)
    return {idx: (px + dx, py + dy) for idx, (px, py) in positions.items()}


def _relative_map_positions(
    mol_entries: list[dict[str, Any]],
    box: tuple[float, float, float, float],
) -> dict[int, tuple[float, float]]:
    x, y, width, height = box
    relative: dict[int, tuple[float, float]] = {}
    for entry in mol_entries:
        for atom_map, (px, py) in entry["positions_by_map"].items():
            if atom_map not in relative:
                relative[atom_map] = ((px - x) / width, (py - y) / height)
    return relative


def _parse_molecule(spec: dict[str, Any]) -> tuple[Any | None, list[str]]:
    warnings: list[str] = []
    if not RDKIT_AVAILABLE:
        return None, ["RDKit is unavailable; mechanism canvas cannot render molecule coordinates."]
    smiles = spec.get("smiles")
    molfile = spec.get("molfile")
    if smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    elif molfile:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False)
    else:
        return None, ["Molecule spec requires smiles or molfile."]
    if mol is None:
        return None, [f"Could not parse molecule spec: {spec.get('label') or smiles or 'molfile'}"]
    try:
        AllChem.Compute2DCoords(mol)
    except Exception as exc:
        warnings.append(f"2D coordinate generation failed: {exc}")
    return mol, warnings


def _map_index(mol: Any) -> dict[int, int]:
    mapping: dict[int, int] = {}
    for atom in mol.GetAtoms():
        number = atom.GetAtomMapNum()
        if number:
            mapping[number] = atom.GetIdx()
    return mapping


def _coords_for(mol: Any) -> dict[int, tuple[float, float]]:
    conf = mol.GetConformer()
    return {atom.GetIdx(): (conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y) for atom in mol.GetAtoms()}


def _hydrogen_text(h_count: int) -> str:
    if h_count <= 0:
        return ""
    if h_count == 1:
        return "H"
    return f"H{h_count}"


def _atom_label(
    atom: Any,
    show_carbons: bool,
    force_label: bool = False,
    atom_pos: tuple[float, float] | None = None,
    neighbor_positions: list[tuple[float, float]] | None = None,
) -> str:
    symbol = atom.GetSymbol()
    if symbol == "C" and atom.GetFormalCharge() == 0 and not (show_carbons or force_label):
        return ""
    h_count = atom.GetTotalNumHs()
    if h_count == 0:
        return symbol
    if symbol == "C" and neighbor_positions and atom_pos and len(neighbor_positions) == 1:
        neighbor_x, _ = neighbor_positions[0]
        if neighbor_x > atom_pos[0]:
            return f"{_hydrogen_text(h_count)}C"
    if h_count == 1:
        return f"{symbol}H"
    return f"{symbol}H{h_count}"


def _atom_label_override(mol_spec: dict[str, Any], atom_map: int | None) -> str | None:
    if not atom_map:
        return None
    labels = mol_spec.get("atom_labels") or mol_spec.get("display_atom_labels") or {}
    if not isinstance(labels, dict):
        return None
    value = labels.get(atom_map, labels.get(str(atom_map)))
    return str(value) if value is not None else None


def _charge_label(charge: int) -> str:
    if charge == 0:
        return ""
    if charge == 1:
        return "+"
    if charge == -1:
        return "-"
    sign = "+" if charge > 0 else "-"
    return f"{abs(charge)}{sign}"


def _estimated_text_width(text: str, font_size: float = 15.0) -> float:
    if not text:
        return 0.0
    return max(font_size * 0.58, len(text) * font_size * 0.52)


def _neighbor_angles(atom_pos: tuple[float, float], neighbor_positions: list[tuple[float, float]]) -> list[float]:
    ax, ay = atom_pos
    return [math.atan2(ny - ay, nx - ax) for nx, ny in neighbor_positions]


def _choose_annotation_angle(
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]],
    preferred_angle: float = -math.pi / 4,
    avoid_positions: list[tuple[float, float]] | None = None,
) -> float:
    """Choose a nearby annotation direction that avoids bonds around an atom."""
    neighbor_angle_values = _neighbor_angles(atom_pos, neighbor_positions)
    ax, ay = atom_pos
    avoid_angle_values = [
        math.atan2(py - ay, px - ax)
        for px, py in (avoid_positions or [])
        if math.hypot(px - ax, py - ay) > 1e-6
    ]
    slots = [
        -math.pi / 2,
        -math.pi / 3,
        -math.pi / 4,
        -math.pi / 6,
        0.0,
        math.pi / 6,
        math.pi / 3,
        math.pi / 2,
        2 * math.pi / 3,
        5 * math.pi / 6,
        math.pi,
        -5 * math.pi / 6,
        -2 * math.pi / 3,
    ]
    scored: list[tuple[float, float]] = []
    for angle in slots:
        min_neighbor_gap = min((_angle_delta(angle, neighbor) for neighbor in neighbor_angle_values), default=math.pi)
        score = min_neighbor_gap * 3.0
        score += max(0.0, math.pi - _angle_delta(angle, preferred_angle)) * 0.35
        for avoid_angle in avoid_angle_values:
            score -= max(0.0, 1.25 - _angle_delta(angle, avoid_angle)) * 8.0
        if min_neighbor_gap < 0.38:
            score -= 6.0
        scored.append((score, angle))
    return max(scored)[1]


def _hydrogen_label_text(h_count: int) -> str:
    return "H" if h_count <= 1 else f"H{h_count}"


def _atom_label_layout(
    atom: Any,
    label: str,
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]] | None = None,
    label_override: str | None = None,
    font_size: float = 15.0,
) -> dict[str, Any]:
    ax, ay = atom_pos
    width = _estimated_text_width(label, font_size)
    half_height = font_size * 0.65
    default = {
        "x": ax,
        "y": ay + 5,
        "anchor": "middle",
        "full_box": (ax - width / 2 - 2, ay - half_height, ax + width / 2 + 2, ay + half_height),
        "core_box": (ax - width / 2 - 2, ay - half_height, ax + width / 2 + 2, ay + half_height),
    }
    if (
        label_override is not None
        or atom.GetSymbol() != "C"
        or atom.GetTotalNumHs() <= 0
        or len(neighbor_positions or []) != 1
    ):
        return default
    c_width = _estimated_text_width("C", font_size)
    if label.endswith("C") and label.startswith("H"):
        text_x = ax + c_width / 2
        return {
            "x": text_x,
            "y": ay + 5,
            "anchor": "end",
            "full_box": (text_x - width - 2, ay - half_height, text_x + 2, ay + half_height),
            "core_box": (ax - c_width / 2 - 2, ay - half_height, ax + c_width / 2 + 2, ay + half_height),
        }
    if label.startswith("C") and "H" in label:
        text_x = ax - c_width / 2
        return {
            "x": text_x,
            "y": ay + 5,
            "anchor": "start",
            "full_box": (text_x - 2, ay - half_height, text_x + width + 2, ay + half_height),
            "core_box": (ax - c_width / 2 - 2, ay - half_height, ax + c_width / 2 + 2, ay + half_height),
        }
    return default


def _atom_label_boxes(
    atom: Any,
    label: str,
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]],
    label_override: str | None = None,
    font_size: float = 15.0,
    avoid_positions: list[tuple[float, float]] | None = None,
) -> list[tuple[float, float, float, float]]:
    if not label:
        return []
    return [
        _atom_label_layout(
            atom,
            label,
            atom_pos,
            neighbor_positions=neighbor_positions,
            label_override=label_override,
            font_size=font_size,
        )["full_box"]
    ]


def _atom_label_core_boxes(
    atom: Any,
    label: str,
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]] | None = None,
    label_override: str | None = None,
    font_size: float = 15.0,
) -> list[tuple[float, float, float, float]]:
    if not label:
        return []
    return [
        _atom_label_layout(
            atom,
            label,
            atom_pos,
            neighbor_positions=neighbor_positions,
            label_override=label_override,
            font_size=font_size,
        )["core_box"]
    ]


def _draw_atom_label(
    parts: list[str],
    atom: Any,
    label: str,
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]],
    color: str,
    label_override: str | None = None,
    font_size: float = 15.0,
    avoid_positions: list[tuple[float, float]] | None = None,
) -> None:
    if not label:
        return
    ax, ay = atom_pos
    layout = _atom_label_layout(
        atom,
        label,
        atom_pos,
        neighbor_positions=neighbor_positions,
        label_override=label_override,
        font_size=font_size,
    )
    parts.append(
        f"<text class='atom' x='{layout['x']:.1f}' y='{layout['y']:.1f}' "
        f"text-anchor='{layout['anchor']}' fill='{color}'>{escape(label)}</text>"
    )


def _normalize_charge_text(label: Any) -> str:
    text = str(label or "").strip().replace("−", "-")
    if text in {"1+", "+1"}:
        return "+"
    if text in {"1-", "-1"}:
        return "-"
    return text


def _draw_formal_charge_marker(
    parts: list[str],
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]],
    label: Any,
    color: str = "#111",
    avoid_positions: list[tuple[float, float]] | None = None,
) -> None:
    text = _normalize_charge_text(label)
    if not text:
        return
    ax, ay = atom_pos
    angle = _choose_annotation_angle(
        (ax, ay),
        neighbor_positions,
        preferred_angle=-math.pi / 4,
        avoid_positions=avoid_positions,
    )
    radius = 6.8 if len(text) <= 1 else 8.8
    offset = 13.0 + max(0.0, radius - 6.8)
    cx = ax + math.cos(angle) * offset
    cy = ay + math.sin(angle) * offset
    parts.append(
        f"<g class='formal-charge-marker' data-anchor-x='{ax:.1f}' data-anchor-y='{ay:.1f}'>"
        f"<circle class='charge-ring' cx='{cx:.1f}' cy='{cy:.1f}' r='{radius:.1f}' stroke='{color}'/>"
        f"<text class='charge-sign' x='{cx:.1f}' y='{cy + 0.4:.1f}' fill='{color}'>{escape(text)}</text>"
        "</g>"
    )


def _draw_bond(
    parts: list[str],
    p1: tuple[float, float],
    p2: tuple[float, float],
    order: int,
    style: dict[str, Any] | None = None,
) -> None:
    x1, y1 = p1
    x2, y2 = p2
    if order == 1:
        parts.append(f"<line class='bond' x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}'/>")
        return
    nx, ny = _unit_normal(x1, y1, x2, y2)
    spacing = float((style or JOURNAL_STYLES["default"]).get("double_bond_spacing", 3.5))
    offsets = [-spacing, spacing] if order == 2 else [-spacing * 1.4, 0.0, spacing * 1.4]
    for offset in offsets:
        cls = "bond" if offset == 0 else "bond-light"
        parts.append(
            f"<line class='{cls}' x1='{x1 + nx * offset:.1f}' y1='{y1 + ny * offset:.1f}' "
            f"x2='{x2 + nx * offset:.1f}' y2='{y2 + ny * offset:.1f}'/>"
        )


def _label_clearance(label: str) -> float:
    if not label:
        return 0.0
    return min(22.0, max(9.0, 4.2 * len(label) + 2.5))


def _shorten_bond_for_labels(
    p1: tuple[float, float],
    p2: tuple[float, float],
    label1: str,
    label2: str,
) -> tuple[tuple[float, float], tuple[float, float]]:
    return _shorten_curve_endpoints(
        p1,
        p2,
        start_gap=_label_clearance(label1),
        end_gap=_label_clearance(label2),
    )


def _draw_lone_pair_at(
    parts: list[str],
    cx: float,
    cy: float,
    angle: float,
    color: str = "#111",
    object_id: str | None = None,
) -> None:
    dx = math.cos(angle + math.pi / 2) * 2.7
    dy = math.sin(angle + math.pi / 2) * 2.7
    if object_id:
        parts.append(
            f"<g class='lone-pair-object' id='{escape(object_id)}' data-object-kind='lone_pair' "
            f"data-anchor-x='{cx:.1f}' data-anchor-y='{cy:.1f}'>"
        )
    parts.append(f"<circle class='lp' cx='{cx - dx:.1f}' cy='{cy - dy:.1f}' r='2.1' style='fill:{color}'/>")
    parts.append(f"<circle class='lp' cx='{cx + dx:.1f}' cy='{cy + dy:.1f}' r='2.1' style='fill:{color}'/>")
    if object_id:
        parts.append("</g>")


def _angle_delta(a: float, b: float) -> float:
    return abs((a - b + math.pi) % (2 * math.pi) - math.pi)


def _candidate_lone_pair_angles(
    atom_pos: tuple[float, float],
    neighbor_positions: list[tuple[float, float]],
    count: int,
    preferred_target: tuple[float, float] | None = None,
) -> list[float]:
    """Choose lone-pair directions as atom-bound annotations, not free canvas dots."""
    ax, ay = atom_pos
    neighbor_angles = [math.atan2(ny - ay, nx - ax) for nx, ny in neighbor_positions]
    preferred_angle = math.atan2(preferred_target[1] - ay, preferred_target[0] - ax) if preferred_target else None
    slots = [index * math.pi / 3 for index in range(6)]
    scored: list[tuple[float, float]] = []
    for angle in slots:
        min_neighbor_gap = min((_angle_delta(angle, neighbor) for neighbor in neighbor_angles), default=math.pi)
        score = min_neighbor_gap
        if preferred_angle is not None:
            score += max(0.0, math.pi - _angle_delta(angle, preferred_angle)) * 1.25
        scored.append((score, angle))
    selected = [angle for _, angle in sorted(scored, reverse=True)[: max(count, 0)]]
    if preferred_angle is not None and selected:
        selected.sort(key=lambda angle: _angle_delta(angle, preferred_angle))
    else:
        selected.sort()
    return selected


def _lone_pair_positions_for_atom(
    atom_map: int,
    mol_entries: list[dict[str, Any]],
    warnings: list[str],
    count: int = 1,
    preferred_target: tuple[float, float] | None = None,
    molecule_index: int | None = None,
) -> list[tuple[float, float, float]]:
    for entry in mol_entries:
        if molecule_index is not None and entry["mol_index"] + 1 != molecule_index:
            continue
        if atom_map not in entry["positions_by_map"]:
            continue
        ax, ay = entry["positions_by_map"][atom_map]
        neighbors = entry.get("neighbor_positions_by_map", {}).get(atom_map, [])
        angles = _candidate_lone_pair_angles((ax, ay), neighbors, count, preferred_target=preferred_target)
        return [(ax + math.cos(angle) * 16, ay + math.sin(angle) * 16, angle) for angle in angles]
    scope = f" in molecule {molecule_index}" if molecule_index is not None else ""
    warnings.append(f"Could not resolve lone-pair atom-map anchor: {atom_map}{scope}")
    return []


def _wrap_text(text: str, width: int = 46) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    length = 0
    for word in words:
        if current and length + len(word) + 1 > width:
            lines.append(" ".join(current))
            current = [word]
            length = len(word)
        else:
            current.append(word)
            length += len(word) + 1
    if current:
        lines.append(" ".join(current))
    return lines or [""]


def _anchor_position(
    anchor: dict[str, Any],
    mol_entries: list[dict[str, Any]],
    warnings: list[str],
) -> tuple[float, float] | None:
    if "xy" in anchor:
        x, y = anchor["xy"]
        return float(x), float(y)
    lone_pair_atom_map = anchor.get("lone_pair_atom_map")
    atom_map = anchor.get("atom_map") or anchor.get("from_atom_map") or anchor.get("to_atom_map")
    bond_map = anchor.get("bond")
    molecule_index = int(anchor["molecule_index"]) if anchor.get("molecule_index") is not None else None
    for entry in mol_entries:
        if molecule_index is not None and entry["mol_index"] + 1 != molecule_index:
            continue
        pos_by_map = entry["positions_by_map"]
        if lone_pair_atom_map and int(lone_pair_atom_map) in pos_by_map:
            ax, ay = pos_by_map[int(lone_pair_atom_map)]
            angle = float(anchor.get("angle", -math.pi / 2))
            return ax + math.cos(angle) * 16, ay + math.sin(angle) * 16
        if atom_map and int(atom_map) in pos_by_map:
            return pos_by_map[int(atom_map)]
        if bond_map:
            a, b = [int(value) for value in bond_map]
            if a in pos_by_map and b in pos_by_map:
                ax, ay = pos_by_map[a]
                bx, by = pos_by_map[b]
                return (ax + bx) / 2, (ay + by) / 2
    warnings.append(f"Could not resolve arrow anchor: {anchor}")
    return None


def _atom_context_for_anchor(
    anchor: dict[str, Any],
    mol_entries: list[dict[str, Any]],
    warnings: list[str],
) -> tuple[tuple[float, float], list[tuple[float, float]]] | None:
    atom_map = anchor.get("atom_map") or anchor.get("from_atom_map") or anchor.get("to_atom_map")
    if not atom_map:
        return None
    molecule_index = int(anchor["molecule_index"]) if anchor.get("molecule_index") is not None else None
    for entry in mol_entries:
        if molecule_index is not None and entry["mol_index"] + 1 != molecule_index:
            continue
        atom_map_int = int(atom_map)
        if atom_map_int not in entry["positions_by_map"]:
            continue
        return (
            entry["positions_by_map"][atom_map_int],
            entry.get("neighbor_positions_by_map", {}).get(atom_map_int, []),
        )
    warnings.append(f"Could not resolve atom context: {anchor}")
    return None


def _mapped_atom_keys(
    atom_map: int,
    mol_entries: list[dict[str, Any]],
    molecule_index: int | None = None,
) -> list[tuple[int, int]]:
    keys: list[tuple[int, int]] = []
    for entry in mol_entries:
        mol_index = entry["mol_index"] + 1
        if molecule_index is not None and mol_index != molecule_index:
            continue
        if atom_map not in entry["positions_by_map"]:
            continue
        keys.append((mol_index, atom_map))
    return keys


def _mapped_bond_keys(
    atom_maps: list[int] | tuple[int, int],
    mol_entries: list[dict[str, Any]],
    molecule_index: int | None = None,
) -> list[tuple[int, int, int]]:
    a, b = sorted(int(value) for value in atom_maps)
    keys: list[tuple[int, int, int]] = []
    for entry in mol_entries:
        mol_index = entry["mol_index"] + 1
        if molecule_index is not None and mol_index != molecule_index:
            continue
        mol = entry["mol"]
        idx_to_map = {idx: atom_map for atom_map, idx in entry["map_to_idx"].items()}
        for bond in mol.GetBonds():
            begin_map = idx_to_map.get(bond.GetBeginAtomIdx())
            end_map = idx_to_map.get(bond.GetEndAtomIdx())
            if begin_map and end_map and sorted((int(begin_map), int(end_map))) == [a, b]:
                keys.append((mol_index, a, b))
                break
    return keys


def _arrow_source_molecule_hint(arrow: dict[str, Any], mol_entries: list[dict[str, Any]]) -> int | None:
    if arrow.get("from_molecule_index") is not None:
        return int(arrow["from_molecule_index"])
    if "from_lone_pair_atom_map" in arrow:
        keys = _mapped_atom_keys(int(arrow["from_lone_pair_atom_map"]), mol_entries)
        return keys[0][0] if len(keys) == 1 else None
    if "from_atom_map" in arrow:
        keys = _mapped_atom_keys(int(arrow["from_atom_map"]), mol_entries)
        return keys[0][0] if len(keys) == 1 else None
    if "from_bond" in arrow:
        keys = _mapped_bond_keys(arrow["from_bond"], mol_entries)
        return keys[0][0] if len(keys) == 1 else None
    return None


def _arrow_target_molecule_hint(arrow: dict[str, Any], mol_entries: list[dict[str, Any]]) -> int | None:
    if arrow.get("to_molecule_index") is not None:
        return int(arrow["to_molecule_index"])
    source_hint = _arrow_source_molecule_hint(arrow, mol_entries)
    if source_hint is not None:
        if "to_atom_map" in arrow and _mapped_atom_keys(int(arrow["to_atom_map"]), mol_entries, source_hint):
            return source_hint
        if "to_bond" in arrow and _mapped_bond_keys(arrow["to_bond"], mol_entries, source_hint):
            return source_hint
    if "to_atom_map" in arrow:
        keys = _mapped_atom_keys(int(arrow["to_atom_map"]), mol_entries)
        return keys[0][0] if len(keys) == 1 else None
    if "to_bond" in arrow:
        keys = _mapped_bond_keys(arrow["to_bond"], mol_entries)
        return keys[0][0] if len(keys) == 1 else None
    return None


def _arrow_atom_keys_for_annotation(arrow: dict[str, Any], mol_entries: list[dict[str, Any]]) -> list[tuple[str, tuple[int, int]]]:
    keys: list[tuple[str, tuple[int, int]]] = []
    if "from_atom_map" in arrow:
        molecule_index = _arrow_source_molecule_hint(arrow, mol_entries)
        for key in _mapped_atom_keys(int(arrow["from_atom_map"]), mol_entries, molecule_index):
            keys.append(("source", key))
    if "from_lone_pair_atom_map" in arrow:
        molecule_index = _arrow_source_molecule_hint(arrow, mol_entries)
        for key in _mapped_atom_keys(int(arrow["from_lone_pair_atom_map"]), mol_entries, molecule_index):
            keys.append(("source", key))
    if "to_atom_map" in arrow:
        molecule_index = _arrow_target_molecule_hint(arrow, mol_entries)
        for key in _mapped_atom_keys(int(arrow["to_atom_map"]), mol_entries, molecule_index):
            keys.append(("target", key))
    return keys


def _anchor_with_molecule_hint(
    anchor: dict[str, Any],
    molecule_hint: int | None,
) -> dict[str, Any]:
    if molecule_hint is None or anchor.get("molecule_index") is not None:
        return anchor
    hinted = dict(anchor)
    hinted["molecule_index"] = molecule_hint
    return hinted


def _ray_exits_box_distance(
    center: tuple[float, float],
    unit: tuple[float, float],
    box: tuple[float, float, float, float],
    padding: float = 3.0,
) -> float | None:
    cx, cy = center
    ux, uy = unit
    x1, y1, x2, y2 = box
    x1 -= padding
    y1 -= padding
    x2 += padding
    y2 += padding
    if not (x1 <= cx <= x2 and y1 <= cy <= y2):
        return None
    distances: list[float] = []
    if ux > 1e-6:
        distances.append((x2 - cx) / ux)
    elif ux < -1e-6:
        distances.append((x1 - cx) / ux)
    if uy > 1e-6:
        distances.append((y2 - cy) / uy)
    elif uy < -1e-6:
        distances.append((y1 - cy) / uy)
    positive = [distance for distance in distances if distance >= 0]
    return min(positive) if positive else None


def _atom_boundary_point(
    atom_context: dict[str, Any],
    toward: tuple[float, float],
    default_radius: float = 8.5,
) -> tuple[float, float]:
    center = atom_context["position"]
    cx, cy = center
    tx, ty = toward
    dx = tx - cx
    dy = ty - cy
    length = math.hypot(dx, dy)
    if length <= 1e-6:
        return center
    unit = (dx / length, dy / length)
    distances = [default_radius]
    for box in atom_context.get("core_label_boxes", atom_context.get("label_boxes", [])):
        distance = _ray_exits_box_distance(center, unit, box, padding=2.5)
        if distance is not None:
            distances.append(distance)
    distance = min(max(max(distances), default_radius), 30.0)
    return cx + unit[0] * distance, cy + unit[1] * distance


def _reaction_split_index(panel: dict[str, Any]) -> int | None:
    molecules = panel.get("molecules", [])
    if not molecules:
        return None
    if "reaction_arrow_after" in panel:
        split = int(panel["reaction_arrow_after"])
        if 0 < split < len(molecules):
            return split
    for index, mol_spec in enumerate(molecules):
        side = str(mol_spec.get("side") or mol_spec.get("role") or "").lower()
        if side in {"product", "products", "expected", "expected_state"} and 0 < index < len(molecules):
            return index
    return None


def _panel_repeated_atom_maps(panel: dict[str, Any]) -> set[int]:
    counts: dict[int, int] = {}
    for mol_spec in panel.get("molecules", []):
        mol, _ = _parse_molecule(mol_spec)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                counts[int(atom_map)] = counts.get(int(atom_map), 0) + 1
    return {atom_map for atom_map, count in counts.items() if count > 1}


def _scoped_svg_object_id(
    panel_index: int,
    mol_index: int,
    kind: str,
    key: str,
    repeated_maps: set[int] | None = None,
) -> str:
    repeated_maps = repeated_maps or set()
    if kind == "atom" and int(key) not in repeated_maps:
        return f"p{panel_index}:atom:{key}"
    if kind == "bond":
        bond_maps = {int(value) for value in key.split("-") if value}
        if not (bond_maps & repeated_maps):
            return f"p{panel_index}:bond:{key}"
    return f"p{panel_index}:m{mol_index + 1}:{kind}:{key}"


def _lone_pair_count_by_map(panel: dict[str, Any]) -> dict[int, int]:
    return {int(item["atom_map"]): int(item.get("count", 1)) for item in panel.get("lone_pairs", []) if "atom_map" in item}


def _lone_pair_key(item: dict[str, Any]) -> tuple[int | None, int] | None:
    if "atom_map" not in item:
        return None
    molecule_index = int(item["molecule_index"]) if item.get("molecule_index") is not None else None
    return molecule_index, int(item["atom_map"])


def _source_lone_pair_keys(panel: dict[str, Any]) -> set[tuple[int | None, int]]:
    keys: set[tuple[int | None, int]] = set()
    for arrow in panel.get("arrows", []):
        if "from_lone_pair_atom_map" not in arrow:
            continue
        molecule_index = int(arrow["from_molecule_index"]) if arrow.get("from_molecule_index") is not None else None
        keys.add((molecule_index, int(arrow["from_lone_pair_atom_map"])))
    return keys


def _matches_lone_pair_source(key: tuple[int | None, int] | None, source_keys: set[tuple[int | None, int]]) -> bool:
    if key is None:
        return False
    molecule_index, atom_map = key
    if key in source_keys or (None, atom_map) in source_keys:
        return True
    return any(source_map == atom_map and (source_molecule is None or molecule_index is None) for source_molecule, source_map in source_keys)


def _visible_lone_pair_count(lp: dict[str, Any], display: dict[str, Any], is_source: bool) -> int:
    for key in ("visible_count", "display_count"):
        if key in lp:
            return max(0, int(lp.get(key) or 0))
    if display.get("publication") and is_source and display.get("lone_pair_display") in {"source", "source_only", "reactive"}:
        return 1
    return max(1, int(lp.get("count", 1)))


def _visible_lone_pair_entries(panel: dict[str, Any], display: dict[str, Any]) -> list[dict[str, Any]]:
    mode = str(display.get("lone_pair_display", "all")).lower()
    source_keys = _source_lone_pair_keys(panel)
    visible: list[dict[str, Any]] = []
    for lp in panel.get("lone_pairs", []):
        key = _lone_pair_key(lp)
        if key is None or lp.get("show") is False or lp.get("visible") is False:
            continue
        forced = lp.get("show") is True or lp.get("always_show") is True
        is_source = _matches_lone_pair_source(key, source_keys)
        if mode in {"none", "hide", "hidden"} and not forced:
            continue
        if mode in {"source", "source_only", "reactive"} and not (forced or is_source):
            continue
        entry = dict(lp)
        count = _visible_lone_pair_count(entry, display, is_source=is_source)
        if count <= 0:
            continue
        entry["_visible_count"] = count
        visible.append(entry)
    return visible


def _has_visible_lone_pair_source(panel: dict[str, Any], arrow: dict[str, Any], display: dict[str, Any]) -> bool:
    if "from_lone_pair_atom_map" not in arrow:
        return True
    molecule_index = int(arrow["from_molecule_index"]) if arrow.get("from_molecule_index") is not None else None
    wanted = (molecule_index, int(arrow["from_lone_pair_atom_map"]))
    for lp in _visible_lone_pair_entries(panel, display):
        if _matches_lone_pair_source(_lone_pair_key(lp), {wanted}):
            return True
    return False


def _object_ref_for_arrow(arrow: dict[str, Any], direction: str, panel_index: int) -> str:
    if direction == "source":
        if "from_lone_pair_atom_map" in arrow:
            if "from_molecule_index" in arrow:
                return f"p{panel_index}:m{int(arrow['from_molecule_index'])}:lp:{int(arrow['from_lone_pair_atom_map'])}:0"
            return f"p{panel_index}:lp:{int(arrow['from_lone_pair_atom_map'])}:0"
        if "from_bond" in arrow:
            a, b = sorted(int(value) for value in arrow["from_bond"])
            if "from_molecule_index" in arrow:
                return f"p{panel_index}:m{int(arrow['from_molecule_index'])}:bond:{a}-{b}"
            return f"p{panel_index}:bond:{a}-{b}"
        if "from_atom_map" in arrow:
            if "from_molecule_index" in arrow:
                return f"p{panel_index}:m{int(arrow['from_molecule_index'])}:atom:{int(arrow['from_atom_map'])}"
            return f"p{panel_index}:atom:{int(arrow['from_atom_map'])}"
        if "from_xy" in arrow:
            return f"p{panel_index}:xy:source"
        return ""
    if "to_bond" in arrow:
        a, b = sorted(int(value) for value in arrow["to_bond"])
        if "to_molecule_index" in arrow:
            return f"p{panel_index}:m{int(arrow['to_molecule_index'])}:bond:{a}-{b}"
        return f"p{panel_index}:bond:{a}-{b}"
    if "to_atom_map" in arrow:
        if "to_molecule_index" in arrow:
            return f"p{panel_index}:m{int(arrow['to_molecule_index'])}:atom:{int(arrow['to_atom_map'])}"
        return f"p{panel_index}:atom:{int(arrow['to_atom_map'])}"
    if "to_xy" in arrow:
        return f"p{panel_index}:xy:target"
    return ""


def _object_ref_for_arrow_with_hints(
    arrow: dict[str, Any],
    direction: str,
    panel_index: int,
    source_hint: int | None = None,
    target_hint: int | None = None,
) -> str:
    hinted = dict(arrow)
    if source_hint is not None and hinted.get("from_molecule_index") is None:
        hinted["from_molecule_index"] = source_hint
    if target_hint is not None and hinted.get("to_molecule_index") is None:
        hinted["to_molecule_index"] = target_hint
    return _object_ref_for_arrow(hinted, direction, panel_index)


def _shorten_curve_endpoints(
    start: tuple[float, float],
    end: tuple[float, float],
    start_gap: float = 8,
    end_gap: float = 13,
) -> tuple[tuple[float, float], tuple[float, float]]:
    sx, sy = start
    ex, ey = end
    dx = ex - sx
    dy = ey - sy
    length = math.hypot(dx, dy)
    if length <= start_gap + end_gap + 1:
        return start, end
    ux = dx / length
    uy = dy / length
    return (sx + ux * start_gap, sy + uy * start_gap), (ex - ux * end_gap, ey - uy * end_gap)


def _distance_point_to_segment(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    px, py = point
    sx, sy = start
    ex, ey = end
    dx = ex - sx
    dy = ey - sy
    length_sq = dx * dx + dy * dy
    if length_sq <= 1e-9:
        return math.hypot(px - sx, py - sy)
    t = max(0.0, min(1.0, ((px - sx) * dx + (py - sy) * dy) / length_sq))
    cx = sx + t * dx
    cy = sy + t * dy
    return math.hypot(px - cx, py - cy)


def _quadratic_point(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    t: float,
) -> tuple[float, float]:
    sx, sy = start
    mx, my = control
    ex, ey = end
    one_minus = 1.0 - t
    return (
        one_minus * one_minus * sx + 2 * one_minus * t * mx + t * t * ex,
        one_minus * one_minus * sy + 2 * one_minus * t * my + t * t * ey,
    )


def _box_penalty(point: tuple[float, float], box: tuple[float, float, float, float], padding: float) -> float:
    px, py = point
    x1, y1, x2, y2 = box
    if x1 - padding <= px <= x2 + padding and y1 - padding <= py <= y2 + padding:
        dx = min(abs(px - (x1 - padding)), abs(px - (x2 + padding)))
        dy = min(abs(py - (y1 - padding)), abs(py - (y2 + padding)))
        return (padding + min(dx, dy) + 1.0) ** 2
    return 0.0


def _curve_overlap_score(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    obstacles: list[dict[str, Any]],
    panel_box: tuple[float, float, float, float],
) -> float:
    x, y, width, height = panel_box
    score = 0.0
    for step in range(4, 17):
        t = step / 20.0
        point = _quadratic_point(start, control, end, t)
        if not (x + 2 <= point[0] <= x + width - 2 and y + 2 <= point[1] <= y + height - 2):
            score += 80.0
        if math.hypot(point[0] - start[0], point[1] - start[1]) < 9:
            continue
        if math.hypot(point[0] - end[0], point[1] - end[1]) < 9:
            continue
        for obstacle in obstacles:
            kind = obstacle.get("kind")
            padding = float(obstacle.get("padding", 5.0))
            if kind == "segment":
                distance = _distance_point_to_segment(point, obstacle["start"], obstacle["end"])
                if distance < padding:
                    score += (padding - distance + 1.0) ** 2
            elif kind == "circle":
                cx, cy = obstacle["center"]
                radius = float(obstacle.get("radius", 8.0)) + padding
                distance = math.hypot(point[0] - cx, point[1] - cy)
                if distance < radius:
                    score += (radius - distance + 1.0) ** 2
            elif kind == "box":
                score += _box_penalty(point, obstacle["box"], padding)
    return score


def _route_arrow_control_point(
    start: tuple[float, float],
    end: tuple[float, float],
    base_curvature: float,
    obstacles: list[dict[str, Any]],
    panel_box: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    sx, sy = start
    ex, ey = end
    nx, ny = _unit_normal(sx, sy, ex, ey)
    sign = -1.0 if base_curvature < 0 else 1.0
    magnitude = max(abs(base_curvature), 16.0)
    candidate_curvatures = [
        base_curvature,
        sign * magnitude * 1.35,
        sign * magnitude * 1.75,
        sign * magnitude * 2.25,
        -sign * magnitude * 1.15,
        -sign * magnitude * 1.6,
        -sign * magnitude * 2.05,
    ]
    best: tuple[float, float, float, float] | None = None
    for curvature in candidate_curvatures:
        mx = (sx + ex) / 2 + nx * curvature
        my = (sy + ey) / 2 + ny * curvature
        score = _curve_overlap_score(start, (mx, my), end, obstacles, panel_box)
        score += abs(curvature - base_curvature) * 0.02
        if best is None or score < best[0]:
            best = (score, mx, my, curvature)
    assert best is not None
    return best[1], best[2], best[3], best[0]


def _auto_curvature(arrow: dict[str, Any], arrow_index: int) -> float:
    if "curvature" in arrow:
        return float(arrow["curvature"])
    source = _arrow_source_kind(arrow)
    target = _arrow_target_kind(arrow)
    base = 24.0
    if source == "bond" and target == "atom":
        base = 16.0
    elif source == "bond" and target == "bond":
        base = 20.0
    elif source == "lone_pair":
        base = 26.0
    side = str(arrow.get("curve_side", "above")).lower()
    sign = 1.0 if side in {"below", "down", "clockwise"} else -1.0
    return sign * base


def _arrow_source_anchor(arrow: dict[str, Any]) -> dict[str, Any]:
    if "from_lone_pair_atom_map" in arrow:
        anchor = {"lone_pair_atom_map": arrow["from_lone_pair_atom_map"]}
        if "from_angle" in arrow:
            anchor["angle"] = arrow["from_angle"]
        if "from_molecule_index" in arrow:
            anchor["molecule_index"] = arrow["from_molecule_index"]
        return anchor
    if "from_bond" in arrow:
        anchor = {"bond": arrow["from_bond"]}
        if "from_molecule_index" in arrow:
            anchor["molecule_index"] = arrow["from_molecule_index"]
        return anchor
    if "from_atom_map" in arrow:
        anchor = {"atom_map": arrow["from_atom_map"]}
        if "from_molecule_index" in arrow:
            anchor["molecule_index"] = arrow["from_molecule_index"]
        return anchor
    if "from_xy" in arrow:
        return {"xy": arrow["from_xy"]}
    return {}


def _arrow_target_anchor(arrow: dict[str, Any]) -> dict[str, Any]:
    if "to_bond" in arrow:
        anchor = {"bond": arrow["to_bond"]}
        if "to_molecule_index" in arrow:
            anchor["molecule_index"] = arrow["to_molecule_index"]
        return anchor
    if "to_atom_map" in arrow:
        anchor = {"atom_map": arrow["to_atom_map"]}
        if "to_molecule_index" in arrow:
            anchor["molecule_index"] = arrow["to_molecule_index"]
        return anchor
    if "to_xy" in arrow:
        return {"xy": arrow["to_xy"]}
    return {}


def _layout_panel_molecules(
    panel: dict[str, Any],
    panel_index: int,
    box: tuple[float, float, float, float],
    warnings: list[str],
    style: dict[str, Any],
    align_targets: dict[int, tuple[float, float]] | None = None,
) -> list[dict[str, Any]]:
    x, y, width, height = box
    molecules = panel.get("molecules", [])
    mol_entries: list[dict[str, Any]] = []
    mol_count = max(len(molecules), 1)
    split_index = _reaction_split_index(panel)
    mol_y_offset = float(panel.get("molecule_y_offset", 42))
    mol_height_fraction = float(panel.get("molecule_height_fraction", 0.54))

    def molecule_box(index: int) -> tuple[float, float, float, float]:
        if split_index is None:
            mol_box_w = (width - 32) / mol_count
            return (x + 16 + index * mol_box_w, y + mol_y_offset, mol_box_w - 10, height * mol_height_fraction)
        left_count = max(split_index, 1)
        right_count = max(mol_count - split_index, 1)
        gap_w = min(max(width * 0.14, 72), 126)
        side_w = (width - 36 - gap_w) / 2
        left_x = x + 18
        right_x = x + 18 + side_w + gap_w
        if index < split_index:
            mol_box_w = side_w / left_count
            return (left_x + index * mol_box_w, y + mol_y_offset, mol_box_w - 8, height * mol_height_fraction)
        mol_box_w = side_w / right_count
        local_index = index - split_index
        return (right_x + local_index * mol_box_w, y + mol_y_offset, mol_box_w - 8, height * mol_height_fraction)

    for index, mol_spec in enumerate(molecules):
        mol, mol_warnings = _parse_molecule(mol_spec)
        warnings.extend(f"Panel {panel_index}: {warning}" for warning in mol_warnings)
        if mol is None:
            continue
        raw_coords = _coords_for(mol)
        mol_box = molecule_box(index)
        map_to_idx = _map_index(mol)
        positions = _transform(raw_coords, mol_box, max_bond_px=float(style["max_bond_px"]))
        if mol_spec.get("flip_x"):
            center_x = mol_box[0] + mol_box[2] / 2
            positions = {atom_idx: (2 * center_x - px, py) for atom_idx, (px, py) in positions.items()}
        if mol_spec.get("flip_y"):
            center_y = mol_box[1] + mol_box[3] / 2
            positions = {atom_idx: (px, 2 * center_y - py) for atom_idx, (px, py) in positions.items()}
        positions = _apply_alignment(positions, map_to_idx, box, align_targets)
        positions_by_map = {map_no: positions[idx] for map_no, idx in map_to_idx.items() if idx in positions}
        idx_to_map = {idx: atom_map for atom_map, idx in map_to_idx.items()}
        neighbor_positions_by_map: dict[int, list[tuple[float, float]]] = {}
        neighbor_positions_by_idx: dict[int, list[tuple[float, float]]] = {}
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            begin_map = idx_to_map.get(begin_idx)
            end_map = idx_to_map.get(end_idx)
            if end_idx in positions:
                neighbor_positions_by_idx.setdefault(begin_idx, []).append(positions[end_idx])
            if begin_idx in positions:
                neighbor_positions_by_idx.setdefault(end_idx, []).append(positions[begin_idx])
            if begin_map and end_idx in positions:
                neighbor_positions_by_map.setdefault(begin_map, []).append(positions[end_idx])
            if end_map and begin_idx in positions:
                neighbor_positions_by_map.setdefault(end_map, []).append(positions[begin_idx])
        mol_entries.append(
            {
                "mol": mol,
                "positions": positions,
                "positions_by_map": positions_by_map,
                "neighbor_positions_by_map": neighbor_positions_by_map,
                "neighbor_positions_by_idx": neighbor_positions_by_idx,
                "map_to_idx": map_to_idx,
                "spec": mol_spec,
                "mol_box": mol_box,
                "mol_index": index,
            }
        )
    return mol_entries


def _arrow_source_kind(arrow: dict[str, Any]) -> str:
    if "from_lone_pair_atom_map" in arrow:
        return "lone_pair"
    if "from_bond" in arrow:
        return "bond"
    if "from_atom_map" in arrow:
        return "atom_center"
    if "from_xy" in arrow:
        return "free_xy"
    return "missing"


def _arrow_target_kind(arrow: dict[str, Any]) -> str:
    if "to_bond" in arrow:
        return "bond"
    if "to_atom_map" in arrow:
        return "atom"
    if "to_xy" in arrow:
        return "free_xy"
    return "missing"


def _panel_atom_maps(panel: dict[str, Any], panel_index: int, validation_warnings: list[str]) -> set[int]:
    atom_maps: set[int] = set()
    for mol_index, mol_spec in enumerate(panel.get("molecules", []), start=1):
        mol, mol_warnings = _parse_molecule(mol_spec)
        validation_warnings.extend(f"Panel {panel_index}, molecule {mol_index}: {warning}" for warning in mol_warnings)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                atom_maps.add(atom_map)
    return atom_maps


def _referenced_atom_maps(arrow: dict[str, Any]) -> list[int]:
    refs: list[int] = []
    for key in ("from_lone_pair_atom_map", "from_atom_map", "to_atom_map"):
        if key in arrow:
            refs.append(int(arrow[key]))
    for key in ("from_bond", "to_bond"):
        if key in arrow:
            refs.extend(int(value) for value in arrow[key])
    return refs


def _reaction_center_atom_maps(panel: dict[str, Any]) -> set[int]:
    maps: set[int] = set()
    for arrow in panel.get("arrows", []):
        maps.update(_referenced_atom_maps(arrow))
    for partial in panel.get("partial_charges", []):
        if "atom_map" in partial:
            maps.add(int(partial["atom_map"]))
    for charge in panel.get("charges", []):
        if "atom_map" in charge:
            maps.add(int(charge["atom_map"]))
    return maps


def _panel_atom_map_set(panel: dict[str, Any], warnings: list[str] | None = None, panel_label: str = "panel") -> set[int]:
    atom_maps: set[int] = set()
    local_warnings = warnings if warnings is not None else []
    for mol_index, mol_spec in enumerate(panel.get("molecules", []), start=1):
        mol, mol_warnings = _parse_molecule(mol_spec)
        local_warnings.extend(f"{panel_label}, molecule {mol_index}: {warning}" for warning in mol_warnings)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                atom_maps.add(int(atom_map))
    return atom_maps


def _panel_mapped_bonds(panel: dict[str, Any]) -> dict[tuple[int, int], str]:
    bonds: dict[tuple[int, int], str] = {}
    for mol_spec in panel.get("molecules", []):
        mol, _ = _parse_molecule(mol_spec)
        if mol is None:
            continue
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom().GetAtomMapNum()
            b = bond.GetEndAtom().GetAtomMapNum()
            if a and b:
                bonds[tuple(sorted((int(a), int(b))))] = str(bond.GetBondType())
    return bonds


def _panel_charge_by_map(panel: dict[str, Any]) -> dict[int, int]:
    charges: dict[int, int] = {}
    for mol_spec in panel.get("molecules", []):
        mol, _ = _parse_molecule(mol_spec)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                charges[int(atom_map)] = int(atom.GetFormalCharge())
    return charges


def _total_formal_charge(panel: dict[str, Any]) -> int:
    return sum(_panel_charge_by_map(panel).values())


def _atom_valence_warnings(panel: dict[str, Any], panel_label: str) -> list[str]:
    warnings: list[str] = []
    max_valence = {"H": 1, "B": 6, "C": 4, "N": 4, "O": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1}
    for mol_index, mol_spec in enumerate(panel.get("molecules", []), start=1):
        mol, _ = _parse_molecule(mol_spec)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in max_valence:
                continue
            explicit_valence = atom.GetExplicitValence()
            if explicit_valence > max_valence[symbol]:
                atom_map = atom.GetAtomMapNum() or atom.GetIdx()
                warnings.append(
                    f"{panel_label}, molecule {mol_index}, atom {atom_map}: explicit valence {explicit_valence} "
                    f"exceeds conservative octet/valence check for {symbol}."
                )
    return warnings


def _graph_edit_atom_maps(edit: dict[str, Any]) -> set[int]:
    maps: set[int] = set()
    for key in ("atom_map", "from_atom_map", "to_atom_map", "donor_atom_map", "acceptor_atom_map"):
        if key in edit:
            maps.add(int(edit[key]))
    for key in ("atom_maps", "bond", "from_bond", "to_bond"):
        if key in edit:
            maps.update(int(value) for value in edit[key])
    return maps


def _validate_graph_edits(
    panel: dict[str, Any],
    panel_label: str,
    atom_maps: set[int],
) -> tuple[list[str], list[str], list[dict[str, Any]]]:
    blocking: list[str] = []
    warnings: list[str] = []
    checks: list[dict[str, Any]] = []
    allowed = {
        "form_bond",
        "break_bond",
        "change_bond",
        "formal_charge",
        "partial_charge",
        "proton_transfer",
        "counterion",
        "radical",
        "resonance",
    }
    for edit_index, edit in enumerate(panel.get("graph_edits", []), start=1):
        edit_type = edit.get("type") or edit.get("operation")
        if not edit_type:
            blocking.append(f"{panel_label}, graph_edit {edit_index}: missing type.")
            continue
        if edit_type not in allowed:
            warnings.append(f"{panel_label}, graph_edit {edit_index}: unknown edit type {edit_type!r}; keep it review-only.")
        referenced_maps = _graph_edit_atom_maps(edit)
        missing_maps = sorted(referenced_maps - atom_maps)
        if missing_maps:
            blocking.append(
                f"{panel_label}, graph_edit {edit_index}: atom-map references not present in panel molecules: {missing_maps}."
            )
        checks.append(
            {
                "panel": panel_label,
                "graph_edit": edit_index,
                "type": edit_type,
                "atom_maps": sorted(referenced_maps),
                "status": "blocked" if missing_maps or not edit_type else "ok",
            }
        )
    return blocking, warnings, checks


def _state_transition_checks(panels: list[dict[str, Any]]) -> dict[str, Any]:
    warnings: list[str] = []
    checks: list[dict[str, Any]] = []
    for index, (left, right) in enumerate(zip(panels, panels[1:]), start=1):
        left_label = left.get("title", f"Panel {index}")
        right_label = right.get("title", f"Panel {index + 1}")
        left_maps = _panel_atom_map_set(left)
        right_maps = _panel_atom_map_set(right)
        missing_in_right = sorted(left_maps - right_maps)
        new_in_right = sorted(right_maps - left_maps)
        if missing_in_right or new_in_right:
            warnings.append(
                f"{left_label} -> {right_label}: atom-map set changes; document omitted proton/counterion/byproduct handling."
            )
        left_charge = _total_formal_charge(left)
        right_charge = _total_formal_charge(right)
        has_charge_context = any(
            item
            for item in (
                left.get("proton_transfers"),
                right.get("proton_transfers"),
                left.get("counterions"),
                right.get("counterions"),
                right.get("charge_explanation"),
            )
        )
        if left_charge != right_charge and not has_charge_context:
            warnings.append(
                f"{left_label} -> {right_label}: total formal charge changes from {left_charge} to {right_charge}; "
                "add proton/counterion/charge_explanation metadata."
            )
        checks.append(
            {
                "from_panel": index,
                "to_panel": index + 1,
                "atom_maps_conserved": not missing_in_right and not new_in_right,
                "formal_charge_before": left_charge,
                "formal_charge_after": right_charge,
                "formal_charge_conserved": left_charge == right_charge,
            }
        )
    return {"warnings": warnings, "checks": checks}


def validate_mechanism_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """Validate whether a mechanism spec has the explicit data a manuscript figure needs."""
    blocking_issues: list[str] = []
    warnings: list[str] = []
    checks: list[dict[str, Any]] = []
    graph_edit_checks: list[dict[str, Any]] = []
    panels = spec.get("panels", [])
    spec_version = str(spec.get("spec_version", "1.0"))
    display = _display_options(spec)
    if not panels:
        blocking_issues.append("Mechanism spec has no panels; provide each elementary step/intermediate explicitly.")
    for panel_index, panel in enumerate(panels, start=1):
        panel_label = panel.get("title", f"Panel {panel_index}")
        molecules = panel.get("molecules", [])
        arrows = panel.get("arrows", [])
        if not molecules:
            blocking_issues.append(f"{panel_label}: missing intermediate/reactant/product molecule.")
        for mol_index, mol_spec in enumerate(molecules, start=1):
            if not (mol_spec.get("smiles") or mol_spec.get("molfile")):
                blocking_issues.append(f"{panel_label}, molecule {mol_index}: provide mapped SMILES or Molfile.")
            mol, mol_warnings = _parse_molecule(mol_spec)
            warnings.extend(f"{panel_label}, molecule {mol_index}: {warning}" for warning in mol_warnings)
            if mol is not None and len(Chem.GetMolFrags(mol)) > 1 and arrows:
                warnings.append(
                    f"{panel_label}, molecule {mol_index}: contains disconnected fragments. "
                    "For mechanism figures, put each reactant/intermediate/byproduct in its own molecule entry so layout and arrows stay chemically readable."
                )
        atom_maps = _panel_atom_maps(panel, panel_index, warnings)
        if arrows and not atom_maps:
            blocking_issues.append(f"{panel_label}: arrows require atom-mapped molecules.")
        if spec_version.startswith("2") and arrows and not panel.get("graph_edits"):
            warnings.append(f"{panel_label}: MechanismSpec v2 should include graph_edits for AI step tracing.")
        edit_blocking, edit_warnings, edit_checks = _validate_graph_edits(panel, panel_label, atom_maps)
        blocking_issues.extend(edit_blocking)
        warnings.extend(edit_warnings)
        graph_edit_checks.extend(edit_checks)
        warnings.extend(_atom_valence_warnings(panel, panel_label))
        for arrow_index, arrow in enumerate(arrows, start=1):
            source = _arrow_source_kind(arrow)
            target = _arrow_target_kind(arrow)
            if source == "missing":
                blocking_issues.append(f"{panel_label}, arrow {arrow_index}: missing electron source anchor.")
            if target == "missing":
                blocking_issues.append(f"{panel_label}, arrow {arrow_index}: missing electron sink anchor.")
            if source == "atom_center":
                warnings.append(
                    f"{panel_label}, arrow {arrow_index}: source is an atom center. "
                    "For publication mechanisms, start arrows at from_lone_pair_atom_map or from_bond."
                )
            if source == "free_xy" or target == "free_xy":
                warnings.append(
                    f"{panel_label}, arrow {arrow_index}: free xy anchors are graphical only; prefer atom-map/bond anchors."
                )
            if "from_lone_pair_atom_map" in arrow and not _has_visible_lone_pair_source(panel, arrow, display):
                blocking_issues.append(
                    f"{panel_label}, arrow {arrow_index}: source lone pair is not visible. "
                    "Add a matching visible lone_pairs entry."
                )
            missing_maps = sorted(set(_referenced_atom_maps(arrow)) - atom_maps)
            if missing_maps:
                blocking_issues.append(
                    f"{panel_label}, arrow {arrow_index}: atom-map anchors not present in panel molecules: {missing_maps}."
                )
            if not arrow.get("kind"):
                warnings.append(f"{panel_label}, arrow {arrow_index}: arrow kind is missing; use electron_pair or single_electron.")
            checks.append(
                {
                    "panel": panel_index,
                    "arrow": arrow_index,
                    "source": source,
                    "target": target,
                    "electron_count": arrow.get("electron_count", 2 if arrow.get("kind") != "single_electron" else 1),
                    "label": arrow.get("label"),
                }
            )
        has_lone_pair_source = any("from_lone_pair_atom_map" in arrow for arrow in arrows)
        if has_lone_pair_source and not _visible_lone_pair_entries(panel, display):
            warnings.append(f"{panel_label}: arrows are present but no lone-pair dots are shown.")
        if arrows and not (panel.get("charges") or panel.get("partial_charges")):
            warnings.append(f"{panel_label}: arrows are present but no formal/partial charge annotations are shown.")
    transition_checks = _state_transition_checks(panels)
    warnings.extend(transition_checks["warnings"])
    return {
        "publication_ready": False,
        "ready_for_human_review": not blocking_issues,
        "blocking_issues": blocking_issues,
        "warnings": warnings,
        "arrow_checks": checks,
        "graph_edit_checks": graph_edit_checks,
        "state_transition_checks": transition_checks["checks"],
    }


def _render_panel(
    panel: dict[str, Any],
    panel_index: int,
    box: tuple[float, float, float, float],
    warnings: list[str],
    style: dict[str, Any],
    display: dict[str, Any] | None = None,
    align_targets: dict[int, tuple[float, float]] | None = None,
) -> tuple[str, list[dict[str, Any]]]:
    display = display or {}
    x, y, width, height = box
    parts: list[str] = [f"<g id='panel-{panel_index}'>"]
    if display.get("show_panel_titles", True):
        parts.append(
            f"<text class='panel-title' x='{x + 12:.1f}' y='{y + 24:.1f}'>{escape(panel.get('title', f'Step {panel_index}'))}</text>"
        )
    molecules = panel.get("molecules", [])
    if not molecules:
        warnings.append(f"Panel {panel_index} has no molecules/intermediates.")
    mol_entries = _layout_panel_molecules(
        panel=panel,
        panel_index=panel_index,
        box=box,
        warnings=warnings,
        style=style,
        align_targets=align_targets,
    )
    mol_count = max(len(molecules), 1)
    split_index = _reaction_split_index(panel)
    repeated_maps = _panel_repeated_atom_maps(panel)
    reaction_center_maps = _reaction_center_atom_maps(panel)
    custom_charge_keys = {
        (int(charge["molecule_index"]) if charge.get("molecule_index") is not None else None, int(charge["atom_map"]))
        for charge in panel.get("charges", [])
        if "atom_map" in charge
    }
    visible_lone_pairs = _visible_lone_pair_entries(panel, display)
    arrow_obstacles: list[dict[str, Any]] = []
    atom_render_contexts: dict[tuple[int, int], dict[str, Any]] = {}
    atom_arrow_avoid_positions: dict[tuple[int, int], list[tuple[float, float]]] = {}
    for arrow in panel.get("arrows", []):
        source_hint = _arrow_source_molecule_hint(arrow, mol_entries)
        target_hint = _arrow_target_molecule_hint(arrow, mol_entries)
        source_pos = _anchor_position(
            _anchor_with_molecule_hint(_arrow_source_anchor(arrow), source_hint),
            mol_entries,
            warnings,
        )
        target_pos = _anchor_position(
            _anchor_with_molecule_hint(_arrow_target_anchor(arrow), target_hint),
            mol_entries,
            warnings,
        )
        for role, key in _arrow_atom_keys_for_annotation(arrow, mol_entries):
            if role == "source" and target_pos:
                atom_arrow_avoid_positions.setdefault(key, []).append(target_pos)
            elif role == "target" and source_pos:
                atom_arrow_avoid_positions.setdefault(key, []).append(source_pos)
    lone_pair_targets: dict[tuple[int | None, int], tuple[float, float]] = {}
    for arrow in panel.get("arrows", []):
        if "from_lone_pair_atom_map" not in arrow:
            continue
        target = _anchor_position(
            _anchor_with_molecule_hint(_arrow_target_anchor(arrow), _arrow_target_molecule_hint(arrow, mol_entries)),
            mol_entries,
            warnings,
        )
        if target:
            molecule_index = int(arrow["from_molecule_index"]) if arrow.get("from_molecule_index") is not None else None
            lone_pair_targets.setdefault((molecule_index, int(arrow["from_lone_pair_atom_map"])), target)
    lone_pair_cache: dict[tuple[int | None, int], list[tuple[float, float, float]]] = {}
    for lp in visible_lone_pairs:
        if "atom_map" not in lp:
            continue
        atom_map = int(lp["atom_map"])
        molecule_index = int(lp["molecule_index"]) if lp.get("molecule_index") is not None else None
        cache_key = (molecule_index, atom_map)
        lone_pair_cache[cache_key] = _lone_pair_positions_for_atom(
            atom_map,
            mol_entries,
            warnings,
            count=int(lp.get("_visible_count", lp.get("count", 1))),
            preferred_target=lone_pair_targets.get(cache_key) or lone_pair_targets.get((None, atom_map)),
            molecule_index=molecule_index,
        )
    for entry in mol_entries:
        index = entry["mol_index"]
        mol_spec = entry["spec"]
        mol = entry["mol"]
        mol_box = entry["mol_box"]
        positions = entry["positions"]
        neighbor_positions_by_idx = entry.get("neighbor_positions_by_idx", {})
        show_carbons = bool(mol_spec.get("show_carbons", False))
        show_atom_maps = bool(mol_spec.get("show_atom_maps", display.get("show_atom_maps", True)))
        atom_labels_by_idx: dict[int, str] = {}
        atom_label_overrides_by_idx: dict[int, str | None] = {}
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            label_override = _atom_label_override(mol_spec, int(atom_map) if atom_map else None)
            atom_label_overrides_by_idx[atom.GetIdx()] = label_override
            atom_position = positions.get(atom.GetIdx())
            atom_neighbors = neighbor_positions_by_idx.get(atom.GetIdx(), [])
            atom_labels_by_idx[atom.GetIdx()] = label_override or _atom_label(
                atom,
                show_carbons=show_carbons,
                force_label=bool(atom_map and atom_map in reaction_center_maps),
                atom_pos=atom_position,
                neighbor_positions=atom_neighbors,
            )
        if mol_count > 1 and index > 0 and panel.get("show_plus", True):
            plus_is_reaction_arrow_slot = split_index is not None and index == split_index
            if plus_is_reaction_arrow_slot:
                pass
            else:
                plus_x = mol_box[0] - 11
                plus_y = y + float(panel.get("operator_y_offset", height * 0.38))
                parts.append(f"<text class='title' x='{plus_x:.1f}' y='{plus_y:.1f}' text-anchor='middle'>+</text>")
        if split_index is not None and index == split_index:
            left_entries = [item for item in mol_entries if item["mol_index"] < split_index]
            right_entries = [item for item in mol_entries if item["mol_index"] >= split_index]
            if left_entries and right_entries:
                arrow_start = max(item["mol_box"][0] + item["mol_box"][2] for item in left_entries) + 16
                arrow_end = min(item["mol_box"][0] for item in right_entries) - 18
                arrow_y = y + float(panel.get("operator_y_offset", height * 0.38))
                if arrow_end > arrow_start + 18:
                    parts.append(
                        f"<path d='M {arrow_start:.1f} {arrow_y:.1f} L {arrow_end:.1f} {arrow_y:.1f}' "
                        "stroke='#111' stroke-width='2.2' marker-end='url(#arrow-step)'/>"
                    )
        for bond in mol.GetBonds():
            begin_map = bond.GetBeginAtom().GetAtomMapNum()
            end_map = bond.GetEndAtom().GetAtomMapNum()
            bond_object_id = ""
            if begin_map and end_map:
                a, b = sorted((int(begin_map), int(end_map)))
                bond_object_id = _scoped_svg_object_id(panel_index, index, "bond", f"{a}-{b}", repeated_maps)
                parts.append(f"<g class='bond-object' id='{escape(bond_object_id)}' data-object-kind='bond' data-atom-maps='{a},{b}'>")
            _draw_bond(
                parts,
                *_shorten_bond_for_labels(
                    positions[bond.GetBeginAtomIdx()],
                    positions[bond.GetEndAtomIdx()],
                    atom_labels_by_idx.get(bond.GetBeginAtomIdx(), ""),
                    atom_labels_by_idx.get(bond.GetEndAtomIdx(), ""),
                ),
                _bond_order(bond),
                style=style,
            )
            arrow_obstacles.append(
                {
                    "kind": "segment",
                    "start": positions[bond.GetBeginAtomIdx()],
                    "end": positions[bond.GetEndAtomIdx()],
                    "padding": 7.0,
                }
            )
            if bond_object_id:
                parts.append("</g>")
        for atom in mol.GetAtoms():
            ax, ay = positions[atom.GetIdx()]
            atom_map = atom.GetAtomMapNum()
            label = atom_labels_by_idx.get(atom.GetIdx(), "")
            label_override = atom_label_overrides_by_idx.get(atom.GetIdx())
            atom_neighbors = neighbor_positions_by_idx.get(atom.GetIdx(), [])
            atom_key = (index + 1, int(atom_map)) if atom_map else None
            annotation_avoid_positions = atom_arrow_avoid_positions.get(atom_key, []) if atom_key else []
            atom_colors = style.get("atom_colors", ATOM_COLORS)
            color = atom_colors.get(atom.GetSymbol(), "#111") if isinstance(atom_colors, dict) else ATOM_COLORS.get(atom.GetSymbol(), "#111")
            atom_group_open = False
            if atom_map:
                atom_object_id = _scoped_svg_object_id(panel_index, index, "atom", str(int(atom_map)), repeated_maps)
                parts.append(
                    f"<g class='atom-object' id='{escape(atom_object_id)}' data-object-kind='atom' "
                    f"data-atom-map='{int(atom_map)}' data-anchor-x='{ax:.1f}' data-anchor-y='{ay:.1f}'>"
                )
                atom_group_open = True
            if label:
                _draw_atom_label(
                    parts,
                    atom,
                    label,
                    (ax, ay),
                    atom_neighbors,
                    color,
                    label_override=label_override,
                    font_size=float(style["atom_font_size"]),
                    avoid_positions=annotation_avoid_positions,
                )
                label_boxes = _atom_label_boxes(
                    atom,
                    label,
                    (ax, ay),
                    atom_neighbors,
                    label_override=label_override,
                    font_size=float(style["atom_font_size"]),
                    avoid_positions=annotation_avoid_positions,
                )
                for label_box in label_boxes:
                    arrow_obstacles.append({"kind": "box", "box": label_box, "padding": 3.0})
                core_label_boxes = _atom_label_core_boxes(
                    atom,
                    label,
                    (ax, ay),
                    neighbor_positions=atom_neighbors,
                    label_override=label_override,
                    font_size=float(style["atom_font_size"]),
                )
            else:
                label_boxes = []
                core_label_boxes = []
                arrow_obstacles.append({"kind": "circle", "center": (ax, ay), "radius": 4.0, "padding": 3.5})
            if atom_key:
                atom_render_contexts[atom_key] = {
                    "position": (ax, ay),
                    "neighbors": atom_neighbors,
                    "atom": atom,
                    "label": label,
                    "label_override": label_override,
                    "label_boxes": label_boxes,
                    "core_label_boxes": core_label_boxes,
                    "color": color,
                    "avoid_positions": annotation_avoid_positions,
                }
            charge = _charge_label(atom.GetFormalCharge()) if display.get("show_formal_charges", True) else ""
            has_custom_charge = (
                (index + 1, int(atom_map)) in custom_charge_keys
                or (None, int(atom_map)) in custom_charge_keys
            )
            if charge and not has_custom_charge:
                _draw_formal_charge_marker(
                    parts,
                    (ax, ay),
                    atom_neighbors,
                    charge,
                    color=color,
                    avoid_positions=annotation_avoid_positions,
                )
            if show_atom_maps and atom_map:
                parts.append(f"<text class='map' x='{ax + 9:.1f}' y='{ay + 17:.1f}'>{atom_map}</text>")
            if atom_group_open:
                parts.append("</g>")
        if display.get("show_molecule_labels", True) and mol_spec.get("label"):
            parts.append(
                f"<text class='subtitle' x='{mol_box[0] + mol_box[2] / 2:.1f}' y='{y + height * 0.64:.1f}' "
                f"text-anchor='middle'>{escape(str(mol_spec['label']))}</text>"
            )
    for lp in visible_lone_pairs:
        atom_map = int(lp.get("atom_map"))
        molecule_index = int(lp["molecule_index"]) if lp.get("molecule_index") is not None else None
        color = lp.get("color", "#111")
        for lp_index, (lp_x, lp_y, lp_angle) in enumerate(lone_pair_cache.get((molecule_index, atom_map), [])):
            if molecule_index is not None:
                object_id = f"p{panel_index}:m{molecule_index}:lp:{atom_map}:{lp_index}"
            else:
                object_id = f"p{panel_index}:lp:{atom_map}:{lp_index}"
            _draw_lone_pair_at(
                parts,
                lp_x,
                lp_y,
                lp_angle,
                color=color,
                object_id=object_id,
            )
    if display.get("show_formal_charges", True):
        for charge in panel.get("charges", []):
            anchor = {"atom_map": int(charge.get("atom_map"))}
            if charge.get("molecule_index") is not None:
                anchor["molecule_index"] = int(charge["molecule_index"])
            molecule_index = int(charge["molecule_index"]) if charge.get("molecule_index") is not None else None
            context_key = (molecule_index, int(charge["atom_map"])) if molecule_index is not None else None
            atom_context = atom_render_contexts.get(context_key) if context_key else None
            if atom_context is None:
                fallback_context = _atom_context_for_anchor(anchor, mol_entries, warnings)
                if fallback_context:
                    pos, neighbors = fallback_context
                    atom_context = {"position": pos, "neighbors": neighbors, "avoid_positions": []}
            if atom_context:
                label = charge.get("label") or charge.get("value") or ""
                _draw_formal_charge_marker(
                    parts,
                    atom_context["position"],
                    atom_context.get("neighbors", []),
                    label,
                    color=charge.get("color", "#111"),
                    avoid_positions=atom_context.get("avoid_positions", []),
                )
    if display.get("show_partial_charges", True):
        for partial in panel.get("partial_charges", []):
            anchor = {"atom_map": int(partial.get("atom_map"))}
            if partial.get("molecule_index") is not None:
                anchor["molecule_index"] = int(partial["molecule_index"])
            pos = _anchor_position(anchor, mol_entries, warnings)
            if pos:
                label = partial.get("label") or ("δ+" if partial.get("sign", "+") == "+" else "δ-")
                parts.append(
                    f"<text class='partial' x='{pos[0] + 14:.1f}' y='{pos[1] + 18:.1f}' fill='{partial.get('color', '#0b63ce')}'>"
                    f"{escape(str(label))}</text>"
                )
    for arrow_index, arrow in enumerate(panel.get("arrows", []), start=1):
        source_hint = _arrow_source_molecule_hint(arrow, mol_entries)
        target_hint = _arrow_target_molecule_hint(arrow, mol_entries)
        start_anchor = _anchor_with_molecule_hint(_arrow_source_anchor(arrow), source_hint)
        end_anchor = _anchor_with_molecule_hint(_arrow_target_anchor(arrow), target_hint)
        end = _anchor_position(end_anchor, mol_entries, warnings)
        start: tuple[float, float] | None
        if "from_lone_pair_atom_map" in arrow and end:
            molecule_index = source_hint
            lp_positions = lone_pair_cache.get((molecule_index, int(arrow["from_lone_pair_atom_map"])), [])
            if not lp_positions and molecule_index is not None:
                lp_positions = lone_pair_cache.get((None, int(arrow["from_lone_pair_atom_map"])), [])
            start = (lp_positions[0][0], lp_positions[0][1]) if lp_positions else None
        else:
            start = _anchor_position(start_anchor, mol_entries, warnings)
        if not start or not end:
            continue
        target_boundary_applied = False
        if "to_atom_map" in arrow:
            context_key = (target_hint, int(arrow["to_atom_map"])) if target_hint is not None else None
            target_context = atom_render_contexts.get(context_key) if context_key else None
            if target_context is None and target_hint is None:
                matches = [
                    context
                    for (mol_index, atom_map), context in atom_render_contexts.items()
                    if atom_map == int(arrow["to_atom_map"])
                ]
                target_context = matches[0] if len(matches) == 1 else None
            if target_context:
                end = _atom_boundary_point(target_context, start)
                target_boundary_applied = True
        if "from_atom_map" in arrow:
            context_key = (source_hint, int(arrow["from_atom_map"])) if source_hint is not None else None
            source_context = atom_render_contexts.get(context_key) if context_key else None
            if source_context is None and source_hint is None:
                matches = [
                    context
                    for (mol_index, atom_map), context in atom_render_contexts.items()
                    if atom_map == int(arrow["from_atom_map"])
                ]
                source_context = matches[0] if len(matches) == 1 else None
            if source_context:
                start = _atom_boundary_point(source_context, end)
        if "start_gap" in arrow:
            start_gap = float(arrow["start_gap"])
        elif "from_lone_pair_atom_map" in arrow:
            start_gap = 2.5
        elif "from_atom_map" in arrow:
            start_gap = 2.0
        elif "from_bond" in arrow:
            start_gap = 1.0
        else:
            start_gap = 8.0
        if "end_gap" in arrow:
            end_gap = float(arrow.get("target_gap", 0.0)) if target_boundary_applied else float(arrow["end_gap"])
        elif "to_atom_map" in arrow and "from_bond" in arrow:
            end_gap = 0.0 if target_boundary_applied else 1.8
        elif "to_atom_map" in arrow:
            end_gap = 0.0 if target_boundary_applied else 2.0
        else:
            end_gap = 9.0
        shortened_start, shortened_end = _shorten_curve_endpoints(start, end, start_gap=start_gap, end_gap=end_gap)
        sx, sy = shortened_start
        ex, ey = shortened_end
        curvature = _auto_curvature(arrow, arrow_index)
        if arrow.get("avoid_molecule_overlap", arrow.get("avoid_molecules", True)):
            mx, my, routed_curvature, overlap_score = _route_arrow_control_point(
                shortened_start,
                shortened_end,
                curvature,
                arrow_obstacles,
                box,
            )
        else:
            nx, ny = _unit_normal(sx, sy, ex, ey)
            routed_curvature = curvature
            mx = (sx + ex) / 2 + nx * curvature
            my = (sy + ey) / 2 + ny * curvature
            overlap_score = _curve_overlap_score(shortened_start, (mx, my), shortened_end, arrow_obstacles, box)
        marker = "arrow-radical" if arrow.get("electron_count") == 1 or arrow.get("kind") == "single_electron" else "arrow-pair"
        color = arrow.get("color", style["arrow_color"])
        parts.append(
            f"<path class='mech-arrow' d='M {sx:.1f} {sy:.1f} Q {mx:.1f} {my:.1f} {ex:.1f} {ey:.1f}' "
            f"stroke='{color}' marker-end='url(#{marker})' data-object-kind='electron_pusher' "
            f"data-source-object='{escape(_object_ref_for_arrow_with_hints(arrow, 'source', panel_index, source_hint, target_hint))}' "
            f"data-target-object='{escape(_object_ref_for_arrow_with_hints(arrow, 'target', panel_index, source_hint, target_hint))}' "
            f"data-electron-count='{arrow.get('electron_count', 2 if arrow.get('kind') != 'single_electron' else 1)}' "
            f"data-routed-curvature='{routed_curvature:.1f}' data-overlap-score='{overlap_score:.2f}'/>"
        )
        if display.get("show_arrow_labels", True) and arrow.get("label"):
            parts.append(f"<text class='callout' x='{mx:.1f}' y='{my - 8:.1f}' text-anchor='middle'>{escape(str(arrow['label']))}</text>")
    if display.get("show_notes", True):
        notes = panel.get("notes", [])
        note_y = y + height - 58
        for note in notes[:3]:
            for line in _wrap_text(str(note), width=58)[:3]:
                parts.append(f"<text class='note' x='{x + width / 2:.1f}' y='{note_y:.1f}' text-anchor='middle'>{escape(line)}</text>")
                note_y += 16
    parts.append("</g>")
    return "\n".join(parts), mol_entries


def mechanism_spec_example() -> dict[str, Any]:
    return {
        "spec_version": "2.0",
        "title": "SN2 mechanism",
        "journal_style": "acs",
        "presentation_mode": "publication",
        "layout": {
            "columns": 1,
            "panel_width": 760,
            "panel_height": 188,
            "top_margin": 10,
            "bottom_margin": 10,
            "align_atom_maps": True,
            "show_atom_maps": False,
            "show_arrow_labels": False,
            "show_molecule_labels": False,
            "show_notes": False,
            "show_footer": False,
            "show_partial_charges": False,
            "lone_pair_display": "source_only",
            "arrow_color": "#111",
            "arrow_width": 1.8,
        },
        "evidence": ["Rule-template SN2 example for renderer validation; not literature proof."],
        "uncertainty": "Educational validation fixture; human review is required for manuscript use.",
        "panels": [
            {
                "title": "Backside attack",
                "starting_state": {"role": "reactants", "molecule_indices": [1, 2], "reaction_smiles": "[OH-:1].[CH3:2][Br:3]"},
                "expected_state": {"role": "products", "molecule_indices": [3, 4], "reaction_smiles": "[OH:1][CH3:2].[Br-:3]"},
                "reaction_arrow_after": 2,
                "molecule_y_offset": 34,
                "molecule_height_fraction": 0.58,
                "operator_y_offset": 94,
                "molecules": [
                    {"smiles": "[OH-:1]", "label": "nucleophile", "side": "reactant", "atom_labels": {"1": "HO"}},
                    {"smiles": "[CH3:2][Br:3]", "label": "alkyl bromide", "side": "reactant", "show_carbons": True},
                    {
                        "smiles": "[OH:1][CH3:2]",
                        "label": "methanol product",
                        "side": "product",
                        "show_carbons": True,
                        "flip_x": True,
                        "atom_labels": {"1": "HO"},
                    },
                    {"smiles": "[Br-:3]", "label": "leaving group", "side": "product"},
                ],
                "lone_pairs": [
                    {"molecule_index": 1, "atom_map": 1, "count": 3, "visible_count": 1},
                ],
                "partial_charges": [
                    {"molecule_index": 2, "atom_map": 2, "label": "δ+"},
                    {"molecule_index": 2, "atom_map": 3, "label": "δ-"},
                ],
                "charges": [
                    {"molecule_index": 1, "atom_map": 1, "label": "−", "color": "#111"},
                    {"molecule_index": 4, "atom_map": 3, "label": "−", "color": "#111"},
                ],
                "arrows": [
                    {
                        "from_lone_pair_atom_map": 1,
                        "from_molecule_index": 1,
                        "to_atom_map": 2,
                        "to_molecule_index": 2,
                        "kind": "electron_pair",
                        "electron_count": 2,
                        "label": "O lp → C",
                        "curvature": -26,
                        "start_gap": 2,
                        "end_gap": 18,
                    },
                    {
                        "from_bond": [2, 3],
                        "from_molecule_index": 2,
                        "to_atom_map": 3,
                        "to_molecule_index": 2,
                        "kind": "electron_pair",
                        "electron_count": 2,
                        "label": "C–Br → Br",
                        "curvature": -22,
                        "start_gap": 1,
                        "end_gap": 5.5,
                    },
                ],
                "graph_edits": [
                    {"type": "form_bond", "atom_maps": [1, 2], "order": 1, "source_arrow": 1},
                    {"type": "break_bond", "atom_maps": [2, 3], "source_arrow": 2},
                    {"type": "formal_charge", "atom_map": 3, "from": 0, "to": -1},
                    {"type": "formal_charge", "atom_map": 1, "from": -1, "to": 0},
                ],
                "bond_changes": {
                    "formed": [{"atom_maps": [1, 2], "order": "SINGLE"}],
                    "broken": [{"atom_maps": [2, 3], "order": "SINGLE"}],
                    "changed": [],
                },
                "charge_changes": [
                    {"atom_map": 1, "from": -1, "to": 0},
                    {"atom_map": 3, "from": 0, "to": -1},
                ],
                "notes": ["Electron pair originates at hydroxide lone pair; leaving-group bond electrons end on bromide."],
            },
        ],
    }


def _molecule_state(mol_spec: dict[str, Any], mol_index: int) -> dict[str, Any]:
    mol, mol_warnings = _parse_molecule(mol_spec)
    state: dict[str, Any] = {
        "molecule_index": mol_index,
        "label": mol_spec.get("label"),
        "smiles": mol_spec.get("smiles"),
        "warnings": mol_warnings,
    }
    if mol is None:
        return state
    atoms = []
    bonds = []
    for atom in mol.GetAtoms():
        atoms.append(
            {
                "atom_map": atom.GetAtomMapNum() or None,
                "symbol": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
                "explicit_valence": atom.GetExplicitValence(),
                "total_hydrogens": atom.GetTotalNumHs(),
            }
        )
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtom().GetAtomMapNum() or bond.GetBeginAtomIdx()
        end = bond.GetEndAtom().GetAtomMapNum() or bond.GetEndAtomIdx()
        bonds.append({"atom_maps": [int(begin), int(end)], "order": str(bond.GetBondType())})
    state.update(
        {
            "atom_maps": sorted(atom["atom_map"] for atom in atoms if atom["atom_map"]),
            "formal_charge": sum(atom["formal_charge"] for atom in atoms),
            "atoms": atoms,
            "bonds": bonds,
        }
    )
    return state


def _panel_state(panel: dict[str, Any], panel_index: int) -> dict[str, Any]:
    molecules = [_molecule_state(mol_spec, index) for index, mol_spec in enumerate(panel.get("molecules", []), start=1)]
    return {
        "panel": panel_index,
        "title": panel.get("title", f"Panel {panel_index}"),
        "molecules": molecules,
        "formal_charge": sum(molecule.get("formal_charge", 0) for molecule in molecules),
        "atom_maps": sorted({atom_map for molecule in molecules for atom_map in molecule.get("atom_maps", [])}),
    }


def _anchor_descriptor(arrow: dict[str, Any], direction: str) -> dict[str, Any]:
    if direction == "source":
        if "from_lone_pair_atom_map" in arrow:
            descriptor = {"kind": "lone_pair", "atom_map": int(arrow["from_lone_pair_atom_map"])}
            if "from_molecule_index" in arrow:
                descriptor["molecule_index"] = int(arrow["from_molecule_index"])
            return descriptor
        if "from_bond" in arrow:
            descriptor = {"kind": "bond", "atom_maps": [int(value) for value in arrow["from_bond"]]}
            if "from_molecule_index" in arrow:
                descriptor["molecule_index"] = int(arrow["from_molecule_index"])
            return descriptor
        if "from_atom_map" in arrow:
            descriptor = {"kind": "atom_center", "atom_map": int(arrow["from_atom_map"])}
            if "from_molecule_index" in arrow:
                descriptor["molecule_index"] = int(arrow["from_molecule_index"])
            return descriptor
        if "from_xy" in arrow:
            return {"kind": "free_xy", "xy": arrow["from_xy"]}
        return {"kind": "missing"}
    if "to_bond" in arrow:
        descriptor = {"kind": "bond", "atom_maps": [int(value) for value in arrow["to_bond"]]}
        if "to_molecule_index" in arrow:
            descriptor["molecule_index"] = int(arrow["to_molecule_index"])
        return descriptor
    if "to_atom_map" in arrow:
        descriptor = {"kind": "atom", "atom_map": int(arrow["to_atom_map"])}
        if "to_molecule_index" in arrow:
            descriptor["molecule_index"] = int(arrow["to_molecule_index"])
        return descriptor
    if "to_xy" in arrow:
        return {"kind": "free_xy", "xy": arrow["to_xy"]}
    return {"kind": "missing"}


def _electron_moves(panel: dict[str, Any]) -> list[dict[str, Any]]:
    moves = []
    for arrow_index, arrow in enumerate(panel.get("arrows", []), start=1):
        moves.append(
            {
                "arrow": arrow_index,
                "kind": arrow.get("kind", "electron_pair"),
                "electron_count": arrow.get("electron_count", 2 if arrow.get("kind") != "single_electron" else 1),
                "source": _anchor_descriptor(arrow, "source"),
                "sink": _anchor_descriptor(arrow, "target"),
                "label": arrow.get("label"),
            }
        )
    return moves


def _inferred_bond_changes(panel: dict[str, Any], next_panel: dict[str, Any] | None) -> dict[str, Any]:
    if next_panel is None:
        return {"formed": [], "broken": [], "changed": []}
    before = _panel_mapped_bonds(panel)
    after = _panel_mapped_bonds(next_panel)
    formed = [
        {"atom_maps": list(key), "order": after[key]}
        for key in sorted(after.keys() - before.keys())
    ]
    broken = [
        {"atom_maps": list(key), "order": before[key]}
        for key in sorted(before.keys() - after.keys())
    ]
    changed = [
        {"atom_maps": list(key), "from": before[key], "to": after[key]}
        for key in sorted(before.keys() & after.keys())
        if before[key] != after[key]
    ]
    return {"formed": formed, "broken": broken, "changed": changed}


def _inferred_charge_changes(panel: dict[str, Any], next_panel: dict[str, Any] | None) -> list[dict[str, Any]]:
    if next_panel is None:
        return []
    before = _panel_charge_by_map(panel)
    after = _panel_charge_by_map(next_panel)
    changes = []
    for atom_map in sorted(before.keys() | after.keys()):
        before_charge = before.get(atom_map)
        after_charge = after.get(atom_map)
        if before_charge != after_charge:
            changes.append({"atom_map": atom_map, "from": before_charge, "to": after_charge})
    return changes


def _trace_state(
    panel: dict[str, Any],
    panel_index: int,
    key: str,
    resolved_panel: dict[str, Any] | None = None,
    resolved_index: int | None = None,
) -> dict[str, Any]:
    resolved = _panel_state(resolved_panel or panel, resolved_index or panel_index)
    declared = panel.get(key)
    if declared:
        return {"declared": declared, "resolved_state": resolved}
    return resolved


def build_mechanism_trace(spec: dict[str, Any], validation: dict[str, Any]) -> dict[str, Any]:
    panels = spec.get("panels", [])
    steps = []
    for panel_index, panel in enumerate(panels, start=1):
        next_panel = panels[panel_index] if panel_index < len(panels) else None
        panel_checks = [check for check in validation.get("arrow_checks", []) if check.get("panel") == panel_index]
        steps.append(
            {
                "step_number": panel_index,
                "title": panel.get("title", f"Panel {panel_index}"),
                "starting_state": _trace_state(panel, panel_index, "starting_state"),
                "electron_moves": _electron_moves(panel),
                "graph_edits": panel.get("graph_edits", []),
                "bond_changes": panel.get("bond_changes") or _inferred_bond_changes(panel, next_panel),
                "charge_changes": panel.get("charge_changes") or _inferred_charge_changes(panel, next_panel),
                "expected_state": _trace_state(
                    panel,
                    panel_index,
                    "expected_state",
                    resolved_panel=next_panel or panel,
                    resolved_index=panel_index + 1 if next_panel else panel_index,
                ),
                "proton_counterion_handling": {
                    "proton_transfers": panel.get("proton_transfers", []),
                    "counterions": panel.get("counterions", []),
                    "charge_explanation": panel.get("charge_explanation"),
                },
                "validation_checks": panel_checks,
                "evidence": panel.get("evidence", spec.get("evidence", [])),
                "uncertainty": panel.get("uncertainty", spec.get("uncertainty", "Mechanistic assignment requires human review.")),
                "human_review_needed": True,
            }
        )
    return {
        "schema": "codex.mechanism_trace.v1",
        "source_spec_version": str(spec.get("spec_version", "1.0")),
        "summary": {
            "step_count": len(steps),
            "ready_for_human_review": validation.get("ready_for_human_review", False),
            "publication_ready": validation.get("publication_ready", False),
        },
        "steps": steps,
        "global_validation": {
            "blocking_issues": validation.get("blocking_issues", []),
            "warnings": validation.get("warnings", []),
            "state_transition_checks": validation.get("state_transition_checks", []),
            "graph_edit_checks": validation.get("graph_edit_checks", []),
        },
    }


def _xml_attr(value: Any) -> str:
    return escape(str(value), quote=True)


def _cdxml_anchor_id(
    arrow: dict[str, Any],
    atom_ids: dict[tuple[int, int], str],
    atom_ids_first: dict[int, str],
    bond_ids: dict[tuple[int, int, int], str],
    bond_ids_first: dict[tuple[int, int], str],
    direction: str,
) -> str | None:
    if direction == "source":
        molecule_index = int(arrow["from_molecule_index"]) if arrow.get("from_molecule_index") is not None else None
        if "from_lone_pair_atom_map" in arrow:
            atom_map = int(arrow["from_lone_pair_atom_map"])
            return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
        if "from_atom_map" in arrow:
            atom_map = int(arrow["from_atom_map"])
            return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
        if "from_bond" in arrow:
            a, b = sorted(int(value) for value in arrow["from_bond"])
            return bond_ids.get((molecule_index, a, b)) if molecule_index is not None else bond_ids_first.get((a, b))
        return None
    molecule_index = int(arrow["to_molecule_index"]) if arrow.get("to_molecule_index") is not None else None
    if "to_atom_map" in arrow:
        atom_map = int(arrow["to_atom_map"])
        return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
    if "to_bond" in arrow:
        a, b = sorted(int(value) for value in arrow["to_bond"])
        return bond_ids.get((molecule_index, a, b)) if molecule_index is not None else bond_ids_first.get((a, b))
    return None


def render_cdxml_document(spec: dict[str, Any], warnings: list[str]) -> str:
    """Render a ChemDraw-like CDXML document with real graph/vector objects.

    The CDXML stays intentionally conservative, but it now writes fragments,
    nodes, bonds, text, reaction-step records, and electron-flow arrow objects
    instead of only embedding the JSON spec as plain text.
    """
    style = _style_for_spec(spec)
    display = _display_options(spec)
    panels = spec.get("panels", [])
    metrics = _canvas_metrics(spec)
    columns = metrics["columns"]
    panel_w = metrics["panel_width"]
    panel_h = metrics["panel_height"]
    top_margin = metrics["top_margin"]
    alignment_by_map: dict[int, tuple[float, float]] = {}
    lines = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        "<CDXML CreationProgram='codex-organic-chem' Name='mechanism-canvas'>",
        "  <fonttable><font id='3' name='Arial'/></fonttable>",
        "  <colortable><color r='1' g='1' b='1'/><color r='0' g='0' b='0'/><color r='0.84' g='0.15' b='0.16'/></colortable>",
        "  <page id='1'>",
        f"    <t id='2' p='36 36'><s font='3' size='14'>{escape(str(spec.get('title', 'Mechanism canvas')))}</s></t>",
    ]
    next_id = 3
    for panel_index, panel in enumerate(panels, start=1):
        row = (panel_index - 1) // columns
        col = (panel_index - 1) % columns
        box = (col * panel_w, top_margin + row * panel_h, panel_w, panel_h)
        panel_warnings: list[str] = []
        entries = _layout_panel_molecules(panel, panel_index, box, panel_warnings, style, alignment_by_map)
        for atom_map, relative in _relative_map_positions(entries, box).items():
            alignment_by_map.setdefault(atom_map, relative)
        panel_atom_ids: dict[tuple[int, int], str] = {}
        panel_atom_ids_first: dict[int, str] = {}
        panel_bond_ids: dict[tuple[int, int, int], str] = {}
        panel_bond_ids_first: dict[tuple[int, int], str] = {}
        fragment_ids: list[str] = []
        arrow_ids: list[str] = []
        next_id += 1
        title_id = str(next_id)
        lines.append(f"    <t id='{title_id}' p='{36 + col * panel_w} {76 + row * panel_h}'><s font='3' size='11'>{escape(str(panel.get('title', f'Panel {panel_index}')))}</s></t>")
        for entry in entries:
            next_id += 1
            fragment_id = str(next_id)
            fragment_ids.append(fragment_id)
            lines.append(f"    <fragment id='{fragment_id}'>")
            mol = entry["mol"]
            positions = entry["positions"]
            atom_idx_to_id: dict[int, str] = {}
            mol_index = entry["mol_index"] + 1
            visible_lone_pairs = _visible_lone_pair_entries(panel, display)
            lone_pair_counts = {
                int(item["atom_map"]): int(item.get("_visible_count", item.get("count", 1)))
                for item in visible_lone_pairs
                if "atom_map" in item and int(item.get("molecule_index", mol_index)) == mol_index
            }
            for atom in mol.GetAtoms():
                next_id += 1
                atom_id = str(next_id)
                atom_idx_to_id[atom.GetIdx()] = atom_id
                atom_map = atom.GetAtomMapNum()
                if atom_map:
                    panel_atom_ids[(mol_index, int(atom_map))] = atom_id
                    panel_atom_ids_first.setdefault(int(atom_map), atom_id)
                ax, ay = positions[atom.GetIdx()]
                attrs = [
                    f"id='{atom_id}'",
                    f"p='{ax:.2f} {ay:.2f}'",
                    f"Element='{_xml_attr(atom.GetSymbol())}'",
                ]
                if atom.GetFormalCharge():
                    attrs.append(f"Charge='{atom.GetFormalCharge()}'")
                if atom_map:
                    attrs.append(f"AtomMap='{atom_map}'")
                    if int(atom_map) in lone_pair_counts:
                        attrs.append(f"NumLonePair='{lone_pair_counts[int(atom_map)]}'")
                lines.append(f"      <n {' '.join(attrs)}/>")
            for bond in mol.GetBonds():
                next_id += 1
                bond_id = str(next_id)
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                a = begin_atom.GetAtomMapNum()
                b = end_atom.GetAtomMapNum()
                if a and b:
                    key = tuple(sorted((int(a), int(b))))
                    panel_bond_ids[(mol_index, key[0], key[1])] = bond_id
                    panel_bond_ids_first.setdefault(key, bond_id)
                lines.append(
                    f"      <b id='{bond_id}' B='{atom_idx_to_id[bond.GetBeginAtomIdx()]}' "
                    f"E='{atom_idx_to_id[bond.GetEndAtomIdx()]}' Order='{_bond_order(bond)}'/>"
                )
            lines.append("    </fragment>")
        for arrow_index, arrow in enumerate(panel.get("arrows", []), start=1):
            resolved_arrow = dict(arrow)
            source_hint = _arrow_source_molecule_hint(arrow, entries)
            target_hint = _arrow_target_molecule_hint(arrow, entries)
            if source_hint is not None and resolved_arrow.get("from_molecule_index") is None:
                resolved_arrow["from_molecule_index"] = source_hint
            if target_hint is not None and resolved_arrow.get("to_molecule_index") is None:
                resolved_arrow["to_molecule_index"] = target_hint
            source_id = _cdxml_anchor_id(
                resolved_arrow,
                panel_atom_ids,
                panel_atom_ids_first,
                panel_bond_ids,
                panel_bond_ids_first,
                "source",
            )
            target_id = _cdxml_anchor_id(
                resolved_arrow,
                panel_atom_ids,
                panel_atom_ids_first,
                panel_bond_ids,
                panel_bond_ids_first,
                "target",
            )
            if not source_id or not target_id:
                continue
            next_id += 1
            arrow_id = str(next_id)
            arrow_ids.append(arrow_id)
            arrow_type = "SingleElectron" if arrow.get("electron_count") == 1 or arrow.get("kind") == "single_electron" else "ElectronPair"
            lines.append(
                f"    <arrow id='{arrow_id}' ArrowType='{arrow_type}' "
                f"ElectronCount='{arrow.get('electron_count', 1 if arrow_type == 'SingleElectron' else 2)}' "
                f"SourceObject='{source_id}' TargetObject='{target_id}' "
                f"Label='{_xml_attr(arrow.get('label', f'arrow {arrow_index}'))}'/>"
            )
        next_id += 1
        lines.append(
            f"    <reactionstep id='{next_id}' Reactants='{' '.join(fragment_ids)}' "
            f"Arrows='{' '.join(arrow_ids)}' StepNumber='{panel_index}'/>"
        )
    if warnings:
        next_id += 1
        warning_text = escape("\n".join(warnings))
        lines.append(f"    <t id='{next_id}' p='36 260'><s font='3' size='8'>Warnings: {warning_text}</s></t>")
    lines.extend(["  </page>", "</CDXML>"])
    return "\n".join(lines) + "\n"


def render_chemdoodle_json(spec: dict[str, Any]) -> dict[str, Any]:
    style = _style_for_spec(spec)
    display = _display_options(spec)
    panels = spec.get("panels", [])
    metrics = _canvas_metrics(spec)
    columns = metrics["columns"]
    panel_w = metrics["panel_width"]
    panel_h = metrics["panel_height"]
    top_margin = metrics["top_margin"]
    molecules_json: list[dict[str, Any]] = []
    shapes: list[dict[str, Any]] = []
    atom_label_overrides: dict[str, str] = {}
    alignment_by_map: dict[int, tuple[float, float]] = {}
    for panel_index, panel in enumerate(panels, start=1):
        row = (panel_index - 1) // columns
        col = (panel_index - 1) % columns
        box = (col * panel_w, top_margin + row * panel_h, panel_w, panel_h)
        entries = _layout_panel_molecules(panel, panel_index, box, [], style, alignment_by_map)
        for atom_map, relative in _relative_map_positions(entries, box).items():
            alignment_by_map.setdefault(atom_map, relative)
        atom_ids_by_map: dict[tuple[int, int], str] = {}
        atom_ids_first_by_map: dict[int, str] = {}
        bond_ids_by_map: dict[tuple[int, int, int], str] = {}
        bond_ids_first_by_map: dict[tuple[int, int], str] = {}
        for entry in entries:
            atoms = []
            bonds = []
            atom_idx_to_order: dict[int, int] = {}
            mol = entry["mol"]
            mol_index = entry["mol_index"] + 1
            visible_lone_pairs = _visible_lone_pair_entries(panel, display)
            lone_pair_counts = {
                int(item["atom_map"]): int(item.get("_visible_count", item.get("count", 1)))
                for item in visible_lone_pairs
                if "atom_map" in item and int(item.get("molecule_index", mol_index)) == mol_index
            }
            for order, atom in enumerate(mol.GetAtoms()):
                atom_id = f"p{panel_index}m{mol_index}a{order}"
                atom_idx_to_order[atom.GetIdx()] = order
                atom_map = atom.GetAtomMapNum()
                if atom_map:
                    atom_ids_by_map[(mol_index, int(atom_map))] = atom_id
                    atom_ids_first_by_map.setdefault(int(atom_map), atom_id)
                x, y = entry["positions"][atom.GetIdx()]
                atom_payload: dict[str, Any] = {"x": round(x, 2), "y": round(y, 2), "i": atom_id}
                if atom.GetSymbol() != "C":
                    atom_payload["l"] = atom.GetSymbol()
                label_override = _atom_label_override(entry["spec"], int(atom_map) if atom_map else None)
                if label_override:
                    atom_label_overrides[atom_id] = label_override
                    atom_payload["display_label"] = label_override
                h_count = atom.GetTotalNumHs()
                if h_count and not label_override:
                    atom_payload["h"] = int(h_count)
                if atom.GetFormalCharge():
                    atom_payload["c"] = atom.GetFormalCharge()
                if atom_map:
                    atom_payload["map"] = int(atom_map)
                    if int(atom_map) in lone_pair_counts:
                        atom_payload["p"] = lone_pair_counts[int(atom_map)]
                        atom_payload["numLonePair"] = lone_pair_counts[int(atom_map)]
                atoms.append(atom_payload)
            for bond_index, bond in enumerate(mol.GetBonds()):
                bond_id = f"p{panel_index}m{entry['mol_index'] + 1}b{bond_index}"
                begin_map = bond.GetBeginAtom().GetAtomMapNum()
                end_map = bond.GetEndAtom().GetAtomMapNum()
                if begin_map and end_map:
                    key = tuple(sorted((int(begin_map), int(end_map))))
                    bond_ids_by_map[(mol_index, key[0], key[1])] = bond_id
                    bond_ids_first_by_map.setdefault(key, bond_id)
                bonds.append(
                    {
                        "b": atom_idx_to_order[bond.GetBeginAtomIdx()],
                        "e": atom_idx_to_order[bond.GetEndAtomIdx()],
                        "i": bond_id,
                        "o": _bond_order(bond),
                    }
                )
            molecules_json.append({"a": atoms, "b": bonds, "panel": panel_index, "label": entry["spec"].get("label")})
        for arrow_index, arrow in enumerate(panel.get("arrows", []), start=1):
            resolved_arrow = dict(arrow)
            source_hint = _arrow_source_molecule_hint(arrow, entries)
            target_hint = _arrow_target_molecule_hint(arrow, entries)
            if source_hint is not None and resolved_arrow.get("from_molecule_index") is None:
                resolved_arrow["from_molecule_index"] = source_hint
            if target_hint is not None and resolved_arrow.get("to_molecule_index") is None:
                resolved_arrow["to_molecule_index"] = target_hint
            source_id = _chemdoodle_anchor_id(
                resolved_arrow,
                atom_ids_by_map,
                atom_ids_first_by_map,
                bond_ids_by_map,
                bond_ids_first_by_map,
                "source",
            )
            target_id = _chemdoodle_anchor_id(
                resolved_arrow,
                atom_ids_by_map,
                atom_ids_first_by_map,
                bond_ids_by_map,
                bond_ids_first_by_map,
                "target",
            )
            if not source_id or not target_id:
                continue
            shapes.append(
                {
                    "i": f"p{panel_index}s{arrow_index}",
                    "t": "Pusher",
                    "o1": source_id,
                    "o2": target_id,
                    "e": arrow.get("electron_count", 2 if arrow.get("kind") != "single_electron" else 1),
                    "panel": panel_index,
                    "label": arrow.get("label"),
                    "source_detail": _arrow_source_kind(arrow),
                    "target_detail": _arrow_target_kind(arrow),
                }
            )
    return {
        "m": molecules_json,
        "s": shapes,
        "metadata": {
            "format": "ChemDoodle JSON",
            "source": "codex-organic-chem mechanism canvas",
            "spec_version": str(spec.get("spec_version", "1.0")),
            "pusher_model": "Pusher shapes reference atom or bond object ids with electron count.",
            "chemdoodle_json_contract": "ChemDoodle JSONInterpreter uses atom field 'p' for numLonePair and Pusher fields 'o1'/'o2'/'e'.",
            "atom_label_overrides": atom_label_overrides,
        },
    }


def render_chemdoodle_html(spec: dict[str, Any], chemdoodle_json: dict[str, Any]) -> str:
    try:
        from .figure_tools import figure_tool_statuses

        chemdoodle_status = figure_tool_statuses()["tools"]["chemdoodle_web"]
        js_path = chemdoodle_status.get("js_path")
        css_path = chemdoodle_status.get("css_path")
    except Exception:  # pragma: no cover - defensive only
        js_path = None
        css_path = None
    json_text = json.dumps(chemdoodle_json, ensure_ascii=False)
    style = _style_for_spec(spec)
    metrics = _canvas_metrics(spec)
    display = _display_options(spec)
    width = metrics["width"]
    height = metrics["height"]
    css_href = css_path or "ChemDoodleWeb.css"
    js_src = js_path or "ChemDoodleWeb.js"
    toolbar_html = (
        "<div class=\"toolbar\">ChemDoodle object-bound mechanism viewer: atoms carry lone-pair counts, Pusher shapes carry source/target object ids.</div>"
        if display["show_chemdoodle_toolbar"]
        else ""
    )
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>{escape(str(spec.get("title", "Mechanism canvas")))}</title>
    <link rel="stylesheet" href="{escape(css_href, quote=True)}" />
    <style>
      body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; }}
      .toolbar {{ border-bottom: 1px solid #ddd; padding: 8px 10px; }}
      #viewer {{ margin: 0; border: 0 !important; display: block; }}
    </style>
  </head>
  <body>
    {toolbar_html}
    <canvas id="viewer" width="{width}" height="{height}"></canvas>
    <script src="{escape(js_src, quote=True)}"></script>
    <script>
      const mechanism = {json_text};
      const content = ChemDoodle.readJSON(JSON.stringify(mechanism));
      const atomLabelOverrides = mechanism.metadata && mechanism.metadata.atom_label_overrides || {{}};
      for (const molecule of content.molecules || []) {{
        for (const atom of molecule.atoms || []) {{
          if (atom.tmpid && atomLabelOverrides[atom.tmpid]) {{
            atom.altLabel = atomLabelOverrides[atom.tmpid];
            atom.implicitH = 0;
          }}
        }}
      }}
      const canvas = new ChemDoodle.ViewerCanvas('viewer', {width}, {height});
      canvas.styles.atoms_displayTerminalCarbonLabels_2D = true;
      canvas.styles.atoms_showAttributedCarbons_2D = true;
      canvas.styles.atoms_font_size_2D = 16;
      canvas.styles.atoms_lonePairDistance_2D = 9;
      canvas.styles.atoms_lonePairSpread_2D = 4.8;
      canvas.styles.atoms_lonePairDiameter_2D = 1.7;
      canvas.styles.atoms_useJMOLColors = true;
      canvas.styles.bonds_width_2D = 1.9;
      canvas.styles.shapes_color = '{escape(str(style["arrow_color"]), quote=True)}';
      canvas.styles.shapes_lineWidth = {float(style["arrow_width"]):.1f};
      canvas.loadContent(content.molecules || [], content.shapes || []);
    </script>
  </body>
</html>
"""


def _chemdoodle_anchor_id(
    arrow: dict[str, Any],
    atom_ids: dict[tuple[int, int], str],
    atom_ids_first: dict[int, str],
    bond_ids: dict[tuple[int, int, int], str],
    bond_ids_first: dict[tuple[int, int], str],
    direction: str,
) -> str | None:
    if direction == "source":
        molecule_index = int(arrow["from_molecule_index"]) if arrow.get("from_molecule_index") is not None else None
        if "from_lone_pair_atom_map" in arrow:
            atom_map = int(arrow["from_lone_pair_atom_map"])
            return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
        if "from_atom_map" in arrow:
            atom_map = int(arrow["from_atom_map"])
            return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
        if "from_bond" in arrow:
            a, b = sorted(int(value) for value in arrow["from_bond"])
            return bond_ids.get((molecule_index, a, b)) if molecule_index is not None else bond_ids_first.get((a, b))
        return None
    molecule_index = int(arrow["to_molecule_index"]) if arrow.get("to_molecule_index") is not None else None
    if "to_atom_map" in arrow:
        atom_map = int(arrow["to_atom_map"])
        return atom_ids.get((molecule_index, atom_map)) if molecule_index is not None else atom_ids_first.get(atom_map)
    if "to_bond" in arrow:
        a, b = sorted(int(value) for value in arrow["to_bond"])
        return bond_ids.get((molecule_index, a, b)) if molecule_index is not None else bond_ids_first.get((a, b))
    return None


def build_ketcher_adapter_payload(spec: dict[str, Any]) -> dict[str, Any]:
    try:
        from .figure_tools import figure_tool_statuses

        ketcher_status = figure_tool_statuses()["tools"]["ketcher"]
    except Exception as exc:  # pragma: no cover - defensive only
        ketcher_status = {"status": "unknown", "message": str(exc)}
    panels = []
    for panel_index, panel in enumerate(spec.get("panels", []), start=1):
        panels.append(
            {
                "panel": panel_index,
                "title": panel.get("title", f"Panel {panel_index}"),
                "mapped_smiles": [mol.get("smiles") for mol in panel.get("molecules", []) if mol.get("smiles")],
                "molfile_count": sum(1 for mol in panel.get("molecules", []) if mol.get("molfile")),
                "atom_mapping_required": bool(panel.get("arrows")),
            }
        )
    return {
        "adapter": "ketcher_input_correction",
        "role": "open-source structure/reaction correction before mechanism reasoning",
        "status": ketcher_status,
        "supported_interop": ["SMILES", "Molfile", "Rxnfile", "CDXML", "KET"],
        "panels": panels,
        "boundaries": [
            "Use Ketcher to correct structures, reactions, and atom maps.",
            "Do not treat Ketcher SVG export alone as validated electron-flow mechanism evidence.",
        ],
    }


def build_marvin_adapter_payload(spec: dict[str, Any]) -> dict[str, Any]:
    try:
        from .figure_tools import figure_tool_statuses

        statuses = figure_tool_statuses()["tools"]
        marvin_status = statuses.get("marvin_js", {})
        marvin_desktop_status = statuses.get("marvin", {})
    except Exception as exc:  # pragma: no cover - defensive only
        marvin_status = {"status": "unknown", "message": str(exc)}
        marvin_desktop_status = {"status": "unknown", "message": str(exc)}
    return {
        "adapter": "marvin_review_target",
        "role": "ChemAxon Marvin review target for CDXML/MOL/RXN/CXSMILES workflows",
        "status": {
            "marvin_js": marvin_status,
            "marvin_desktop": marvin_desktop_status,
        },
        "primary_artifact": "mechanism.cdxml",
        "supported_interop": ["CDXML", "MOL/RXN", "SMILES", "CXSMILES/CXON"],
        "panels": [
            {
                "panel": index,
                "title": panel.get("title", f"Panel {index}"),
                "mapped_smiles": [mol.get("smiles") for mol in panel.get("molecules", []) if mol.get("smiles")],
            }
            for index, panel in enumerate(spec.get("panels", []), start=1)
        ],
        "boundaries": [
            "Marvin JS/Desktop is proprietary ChemAxon software and cannot be bundled without the user's licensed package.",
            "Use mechanism.cdxml plus mapped SMILES/MOL/RXN data for Marvin-side human review when Marvin is installed.",
        ],
    }


def render_mechanism_canvas(spec: dict[str, Any], output_dir: str | None = None) -> dict[str, Any]:
    warnings: list[str] = []
    validation = validate_mechanism_spec(spec)
    warnings.extend(validation["warnings"])
    style = _style_for_spec(spec)
    layout = spec.get("layout", {})
    display = _display_options(spec)
    panels = spec.get("panels", [])
    metrics = _canvas_metrics(spec)
    columns = metrics["columns"]
    panel_w = metrics["panel_width"]
    panel_h = metrics["panel_height"]
    width = metrics["width"]
    height = metrics["height"]
    top_margin = metrics["top_margin"]
    parts = _svg_header(width, height, style=style)
    title = spec.get("title", "Detailed stepwise mechanism")
    subtitle = spec.get("subtitle")
    if display["show_title"]:
        parts.append(f"<text class='title' x='{width / 2:.1f}' y='34' text-anchor='middle'>{escape(str(title))}</text>")
    if subtitle and display["show_subtitle"]:
        parts.append(f"<text class='subtitle' x='{width / 2:.1f}' y='58' text-anchor='middle'>{escape(str(subtitle))}</text>")
    if not panels:
        warnings.append("Mechanism spec has no panels/intermediates.")
    alignment_by_map: dict[int, tuple[float, float]] = {}
    for idx, panel in enumerate(panels):
        row = idx // columns
        col = idx % columns
        box = (col * panel_w, top_margin + row * panel_h, panel_w, panel_h)
        align_targets = alignment_by_map if layout.get("align_atom_maps", True) else None
        panel_svg, mol_entries = _render_panel(
            panel,
            idx + 1,
            box,
            warnings,
            style=style,
            display=display,
            align_targets=align_targets,
        )
        parts.append(panel_svg)
        for atom_map, relative in _relative_map_positions(mol_entries, box).items():
            alignment_by_map.setdefault(atom_map, relative)
        if idx < len(panels) - 1:
            arrow_x = col * panel_w + panel_w - 36
            arrow_y = top_margin + row * panel_h + panel_h * 0.42
            if col < columns - 1:
                parts.append(
                    f"<path d='M {arrow_x - 44:.1f} {arrow_y:.1f} L {arrow_x:.1f} {arrow_y:.1f}' "
                    "stroke='#111' stroke-width='3' marker-end='url(#arrow-step)'/>"
                )
    if display["show_footer"]:
        parts.append(
            f"<text class='callout' x='{width / 2:.1f}' y='{height - 20:.1f}' text-anchor='middle'>"
            "Mechanism canvas generated from explicit intermediates, atom-map anchors, charges, partial charges, and lone-pair annotations."
            "</text>"
        )
    parts.append("</svg>")
    svg = "\n".join(parts)
    cdxml = render_cdxml_document(spec, warnings)
    chemdoodle_json = render_chemdoodle_json(spec)
    chemdoodle_html = render_chemdoodle_html(spec, chemdoodle_json)
    mechanism_trace = build_mechanism_trace(spec, validation)
    ketcher_adapter = build_ketcher_adapter_payload(spec)
    marvin_adapter = build_marvin_adapter_payload(spec)
    status = "ok"
    if validation["blocking_issues"]:
        status = "blocked_for_publication"
    elif warnings:
        status = "ok_with_warnings"
    result = {
        "status": status,
        "spec_version": str(spec.get("spec_version", "1.0")),
        "mechanism_accuracy_mode": "explicit_intermediate_atom_mapped_canvas",
        "journal_style": style["name"],
        "svg": svg,
        "cdxml": cdxml,
        "cdxml_skeleton": cdxml,
        "chemdoodle_json": chemdoodle_json,
        "chemdoodle_json_string": json.dumps(chemdoodle_json, ensure_ascii=False, indent=2),
        "chemdoodle_html": chemdoodle_html,
        "ketcher_adapter": ketcher_adapter,
        "marvin_adapter": marvin_adapter,
        "mechanism_trace": mechanism_trace,
        "mechanism_trace_json": json.dumps(mechanism_trace, ensure_ascii=False, indent=2),
        "publication_checks": validation,
        "warnings": warnings,
        "requirements": [
            "Every chemically meaningful arrow should be anchored to atom_map or bond atom-map pairs.",
            "Every intermediate must be supplied as SMILES/Molfile; the renderer does not infer unknown intermediates.",
            "Show only mechanism-relevant lone pairs by default; use visible_count/all-lone-pair display for teaching figures.",
            "Show formal charges and partial charges explicitly when they clarify the electron flow.",
            "Use the mechanism_trace as an auditable step trace for AI reasoning rather than relying on a single overall reaction.",
            "Use ChemDraw/ChemDoodle/Marvin/Illustrator/Inkscape for final manual typography and journal style cleanup.",
        ],
        "common_literature_workflow": {
            "chemical_canvas": ["ChemDraw/CDXML", "ChemDoodle", "ChemAxon Marvin/Marvin JS"],
            "final_polish": ["ChemDraw style sheets", "Adobe Illustrator", "Inkscape"],
            "why": [
                "Mechanism arrows are graphical electron-flow objects with chemical anchors.",
                "Electron-source lone pairs, charges, radicals, atom mapping, and partial charges are explicit annotations.",
                "Journal graphics are normally vector artwork with journal-specific font and bond settings.",
            ],
        },
    }
    if output_dir:
        from .figure_tools import figure_followup_for_files

        path = Path(output_dir).expanduser().resolve()
        path.mkdir(parents=True, exist_ok=True)
        svg_path = path / "mechanism.svg"
        cdxml_path = path / "mechanism.cdxml"
        spec_path = path / "mechanism.spec.json"
        trace_path = path / "mechanism.trace.json"
        chemdoodle_path = path / "mechanism.chemdoodle.json"
        chemdoodle_html_path = path / "mechanism.chemdoodle.html"
        ketcher_path = path / "mechanism.ketcher.json"
        marvin_path = path / "mechanism.marvin.json"
        svg_path.write_text(svg, encoding="utf-8")
        cdxml_path.write_text(cdxml, encoding="utf-8")
        spec_path.write_text(json.dumps(spec, ensure_ascii=False, indent=2), encoding="utf-8")
        trace_path.write_text(result["mechanism_trace_json"], encoding="utf-8")
        chemdoodle_path.write_text(result["chemdoodle_json_string"], encoding="utf-8")
        chemdoodle_html_path.write_text(chemdoodle_html, encoding="utf-8")
        ketcher_path.write_text(json.dumps(ketcher_adapter, ensure_ascii=False, indent=2), encoding="utf-8")
        marvin_path.write_text(json.dumps(marvin_adapter, ensure_ascii=False, indent=2), encoding="utf-8")
        result["output_files"] = [
            str(svg_path),
            str(cdxml_path),
            str(spec_path),
            str(trace_path),
            str(chemdoodle_path),
            str(chemdoodle_html_path),
            str(ketcher_path),
            str(marvin_path),
        ]
        result["editor_followup"] = figure_followup_for_files(svg_path=str(svg_path), cdxml_path=str(cdxml_path))
    return result


def render_cdxml_skeleton(spec: dict[str, Any], warnings: list[str]) -> str:
    """Backward-compatible alias for callers that used the old name."""
    return render_cdxml_document(spec, warnings)
