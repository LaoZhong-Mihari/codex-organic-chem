"""Composed route/scheme figure renderer.

Molecules are drawn through :mod:`figure_style`, so every structure in a route
shares one bond length, line width, and font size regardless of its size — the
property that makes a multi-step scheme read as a single drawing. Cells are
sized to each molecule's ink and caption instead of a fixed 345x300 box, and
rows are packed from those real sizes, so captions and neighbouring structures
cannot land on top of one another by construction.

Three spec options make a multi-step scheme read like one drawn by hand:

* ``align_structures`` (default on): a structure given as SMILES is laid out
  with the part it shares with one of the previous structures (largest common
  substructure, complete rings only) held at that structure's 2D coordinates,
  so a protected core does not flip or turn from step to step. Molfile input
  keeps its own coordinates. ``"align": false`` on an item opts it out.
* ``subscript_formulas`` (default on): digits that follow a letter or a closing
  bracket in arrow text become subscripts (CH2Cl2 -> CH₂Cl₂).
* ``max_width`` (default off): a row wider than this many pixels is split before
  arrows into the fewest rows that fit, with widths as even as possible; each
  continuation row starts with its arrow.

Every render finishes with :func:`figure_audit.audit_figure`; its result is
returned as ``publication_checks`` and any blocking issue downgrades ``status``.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from .figure_audit import FigureElement, audit_figure, svg_to_png
from .figure_style import _normalize_bond_units, prepare_mol, render_molecule, resolve_preset
from .rdkit_tools import RDKIT_AVAILABLE, normalize_structure

MOL_GAP = 26
ROW_GAP = 64
TOP_MARGIN = 40
LEFT_MARGIN = 40
RIGHT_MARGIN = 40
LABEL_BLOCK = 52
FOOTER_LINE_HEIGHT = 24
MIN_MOL_CELL = 96
ARROW_DEFAULT_WIDTH = 125
CHAR_W = 0.56  # approximate Arial advance per px of font size
CHAR_W_BOLD = 0.6

# Captions under a structure: a bold label (usually the compound number) and a
# grey code line (usually a name). Offsets are from the bottom of the row's
# tallest structure; the defaults reproduce the one-line layout (24 + 19 + 9 = 52).
LABEL_PX = 15
CODE_PX = 13
CAPTION_TOP = 24
LABEL_LINE = 18
CODE_GAP = 19
CODE_LINE = 17
CAPTION_DESCENT = 9
CAPTION_MAX_PX = 220  # captions wrap at the wider of this and their structure

# Layout alignment with earlier structures.
ALIGN_LOOKBACK = 3  # how many earlier structures may serve as the template
ALIGN_MIN_ATOMS = 5
ALIGN_MIN_FRACTION = 0.5  # shared atoms as a fraction of the smaller structure
ALIGN_TIMEOUT_S = 2
RDKIT_DEPICT_BOND = 1.5  # Compute2DCoords bond length; prepare_mol rescales to 1.0
ALIGN_MIN_SEPARATION = 0.45  # bond lengths; closer non-bonded atoms mean a cramped layout
ALIGN_MAX_BOND_RATIO = 1.3

_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_FORMULA_DIGITS = re.compile(r"(?<=[A-Za-z)\]])\d+")
# Tokens whose digits are not atom counts: mechanism labels (SN2, SN1', E1cb) and
# numbered ligands, catalysts or substituents (L1, R2, X1, G2, HG2).
_NOT_FORMULA = re.compile(r"^(?:S[NE]\d'?|E\d(?:cb)?|H?G\d|[LRX]\d+[a-z]?)$")


def route_figure_spec_example() -> dict[str, Any]:
    return {
        "title": "Example protected route",
        "subtitle": "Composed route figure from validated molecule SMILES.",
        "style_preset": "acs",
        "rows": [
            {
                "items": [
                    {
                        "id": "starting_material",
                        "label": "Starting material",
                        "code": "S0",
                        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    },
                    {"type": "arrow", "label": "1 equiv reagent", "sublabel": "monitor conversion"},
                    {
                        "id": "product",
                        "label": "Product candidate",
                        "code": "P",
                        "smiles": "O=C(O)c1ccccc1O",
                    },
                ]
            }
        ],
        "footnotes": ["Example route figure; replace labels and SMILES with reviewed structures."],
    }


def subscript_formulas(text: str) -> str:
    """Write atom counts in reagent text as subscripts: CH2Cl2 -> CH₂Cl₂.

    Only digits right after a letter or a closing bracket change, so amounts,
    temperatures and locants (0 °C, 1 atm, 10 mol%, 2,2-dimethoxypropane) stay
    as written, and so do mechanism labels and numbered ligands (SN2, E1cb, L1).
    """

    def fix(match: re.Match[str]) -> str:
        token = match.group(0)
        if _NOT_FORMULA.match(token.strip("()[]{}.:")):
            return token
        return _FORMULA_DIGITS.sub(lambda digits: digits.group(0).translate(_SUBSCRIPT), token)

    return re.sub(r"[^\s,;/]+", fix, text)


def render_route_figure(
    spec: dict[str, Any],
    output_dir: str | None = None,
    basename: str = "route_figure",
    formats: list[str] | None = None,
) -> dict[str, Any]:
    formats = formats or ["svg"]
    unsupported = sorted(set(formats) - {"svg", "png"})
    if unsupported:
        return {"status": "error", "warnings": [f"Unsupported route figure format(s): {', '.join(unsupported)}"]}
    if not RDKIT_AVAILABLE:
        return {"status": "error", "warnings": ["RDKit is unavailable; route figure rendering requires RDKit."]}

    rows = spec.get("rows")
    if not isinstance(rows, list) or not rows:
        return {"status": "error", "warnings": ["Route figure spec requires a non-empty 'rows' list."]}

    preset = resolve_preset(
        spec.get("style_preset") or spec.get("preset"),
        overrides={"bond_length_px": spec.get("bond_length_px")},
    )
    options = {
        "align": bool(spec.get("align_structures", True)),
        "subscripts": bool(spec.get("subscript_formulas", True)),
    }
    warnings: list[str] = []
    normalized: list[dict[str, Any]] = []
    previous: list[dict[str, Any]] = []
    try:
        prepared_rows = [_prepare_row(row, normalized, warnings, preset, previous, options) for row in rows]
        max_width = int(spec.get("max_width") or 0)
        if max_width > 0:
            prepared_rows = [
                {"items": part, "note": row.get("note") if index == 0 else None}
                for row in prepared_rows
                for index, part in enumerate(_split_row(row["items"], max_width))
            ]
        svg, elements, width, height, row_widths = _compose_svg(spec, prepared_rows)
    except ValueError as exc:
        return {"status": "error", "warnings": [str(exc), *warnings]}

    output_files: list[str] = []
    png_path: Path | None = None
    if output_dir:
        out_dir = Path(output_dir).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        svg_path = out_dir / f"{basename}.svg"
        svg_path.write_text(svg, encoding="utf-8")
        output_files.append(str(svg_path))
        if "png" in formats:
            png_path = out_dir / f"{basename}.png"
            raster_warnings = svg_to_png(svg_path, png_path)
            warnings.extend(raster_warnings)
            if png_path.exists():
                output_files.append(str(png_path))
            else:
                png_path = None

    publication_checks = audit_figure(
        renderer="codex_route_figure_composed_svg",
        width=width,
        height=height,
        elements=elements,
        target_bond_px=float(preset["bond_length_px"]),
        png_path=png_path,
        extra_warnings=[w for w in warnings if "bond length" in w],
    )

    status = "ok"
    if publication_checks["blocking_issues"]:
        status = "blocked_for_publication"
    elif warnings or publication_checks["warnings"]:
        status = "ok_with_warnings"
    if "png" in formats and output_dir and png_path is None:
        status = "ok_with_warnings" if status == "ok" else status

    return {
        "status": status,
        "kind": "route_figure",
        "formats": formats,
        "style_preset": preset["name"],
        "svg": svg,
        "output_files": output_files,
        "molecules": normalized,
        "layout": {"rows": len(row_widths), "row_widths": row_widths, "canvas": [width, height]},
        "publication_checks": publication_checks,
        "warnings": warnings,
    }


def _prepare_row(
    row: dict[str, Any],
    normalized: list[dict[str, Any]],
    warnings: list[str],
    preset: dict[str, Any],
    previous: list[dict[str, Any]],
    options: dict[str, bool],
) -> dict[str, Any]:
    items = row.get("items")
    if not isinstance(items, list) or not items:
        raise ValueError("Each route figure row requires a non-empty 'items' list.")
    prepared_items: list[dict[str, Any]] = []
    for raw in items:
        if not isinstance(raw, dict):
            raise ValueError("Route figure row items must be objects.")
        item_type = raw.get("type") or ("molecule" if raw.get("smiles") or raw.get("molfile") else None)
        if item_type == "molecule":
            prepared = _prepare_molecule(raw, preset, previous if options["align"] else None)
            prepared_items.append(prepared)
            previous.append({"mol": prepared["mol"], "name": prepared["label"] or prepared["id"] or prepared["code"]})
            record = prepared["record"]
            entry = {
                "id": prepared.get("id"),
                "label": prepared.get("label"),
                "code": prepared.get("code"),
                "input_smiles": raw.get("smiles"),
                "canonical_smiles": record["canonical_smiles"],
                "isomeric_smiles": record["isomeric_smiles"],
                "inchi_key": record.get("inchi_key"),
                "warnings": record.get("warnings", []),
                "metadata": record.get("metadata", {}),
            }
            if prepared["alignment"]:
                entry.update(prepared["alignment"])
            normalized.append(entry)
            warnings.extend(record.get("warnings", []))
            warnings.extend(prepared["depiction"].warnings)
        elif item_type == "arrow":
            label = str(raw.get("label", ""))
            sublabel = str(raw.get("sublabel") or raw.get("label2") or "")
            if options["subscripts"]:
                label, sublabel = subscript_formulas(label), subscript_formulas(sublabel)
            width = int(raw.get("width", 0)) or _arrow_width_for(label, sublabel)
            prepared_items.append({"type": "arrow", "label": label, "sublabel": sublabel, "width": width})
        elif item_type == "plus":
            prepared_items.append({"type": "plus", "width": int(raw.get("width", 44))})
        elif item_type == "text":
            text = str(raw.get("text", ""))
            width = int(raw.get("width", 0)) or max(60, int(len(text) * 14 * CHAR_W) + 24)
            prepared_items.append({"type": "text", "text": text, "width": width})
        else:
            raise ValueError("Route figure item requires type 'molecule', 'arrow', 'plus', or 'text'.")
    return {"items": prepared_items, "note": row.get("note")}


def _arrow_width_for(label: str, sublabel: str) -> int:
    # The arrow must clear its own condition labels, or long reagent lines
    # spill onto the neighbouring structures.
    longest = max(len(label), len(sublabel))
    return max(ARROW_DEFAULT_WIDTH, int(longest * 14 * CHAR_W) + 36)


def _text_width(text: str, px: float, bold: bool = False) -> float:
    return len(text) * px * (CHAR_W_BOLD if bold else CHAR_W)


def _wrap_caption(text: str, px: float, max_px: float, bold: bool = False) -> list[str]:
    """Greedy word wrap; a word longer than a line is broken after its hyphens."""
    lines: list[str] = []
    for word in text.split():
        parts = re.findall(r"[^-]*-|[^-]+$", word) if _text_width(word, px, bold) > max_px else [word]
        for index, part in enumerate(parts):
            joined = lines[-1] + ("" if index else " ") + part if lines else part
            if lines and _text_width(joined, px, bold) <= max_px:
                lines[-1] = joined
            else:
                lines.append(part)
    return lines


def _caption_lines(label: str, code: str, max_px: float) -> list[dict[str, Any]]:
    """Caption lines with their baseline offset below the structure row."""
    lines: list[dict[str, Any]] = []
    offset = CAPTION_TOP
    for index, text in enumerate(_wrap_caption(label, LABEL_PX, max_px, bold=True) if label else []):
        offset = CAPTION_TOP + index * LABEL_LINE
        lines.append({"text": text, "dy": offset, "cls": "molLabel", "width": _text_width(text, LABEL_PX, True), "px": LABEL_PX})
    first_code = offset + CODE_GAP if lines else CAPTION_TOP
    for index, text in enumerate(_wrap_caption(code, CODE_PX, max_px) if code else []):
        lines.append(
            {"text": text, "dy": first_code + index * CODE_LINE, "cls": "code", "width": _text_width(text, CODE_PX), "px": CODE_PX}
        )
    return lines


def _prepare_molecule(raw: dict[str, Any], preset: dict[str, Any], previous: list[dict[str, Any]] | None) -> dict[str, Any]:
    smiles = raw.get("smiles")
    molfile = raw.get("molfile")
    name = raw.get("id") or raw.get("label") or "<unnamed>"
    if not smiles and not molfile:
        raise ValueError("Molecule route figure item requires 'smiles' or 'molfile'.")
    record = normalize_structure(smiles=smiles, molfile=molfile, source="route_figure").to_dict()
    if not record.get("canonical_smiles") and not record.get("molfile"):
        raise ValueError(f"Invalid molecule for route figure item {name}.")
    alignment: dict[str, Any] | None = None
    try:
        mol = None
        if molfile:
            mol = prepare_mol(molfile=molfile)  # keep the coordinates the user drew
        elif previous and raw.get("align", True):
            from rdkit import Chem

            base = Chem.MolFromSmiles(record.get("isomeric_smiles") or smiles)
            if base is not None:
                mol, alignment = _align_to_previous(base, previous)
        if mol is None:
            mol = prepare_mol(smiles=record.get("isomeric_smiles") or smiles)
        depiction = render_molecule(mol=mol, preset=preset)
    except ValueError as exc:
        raise ValueError(f"RDKit could not draw route figure item {name}: {exc}") from exc
    label = str(raw.get("label") or raw.get("title") or raw.get("id") or "")
    code = str(raw.get("code") or raw.get("id") or "")
    if code == label:
        code = ""  # an id used for both lines would be printed twice
    captions = _caption_lines(label, code, max(float(depiction.width), CAPTION_MAX_PX))
    caption_width = max((line["width"] for line in captions), default=0.0)
    return {
        "type": "molecule",
        "id": str(raw.get("id") or raw.get("code") or raw.get("label") or ""),
        "label": label,
        "code": code,
        "width": max(depiction.width, MIN_MOL_CELL, math.ceil(caption_width) + 8),
        "height": depiction.height,
        "captions": captions,
        "caption_height": max((line["dy"] for line in captions), default=0) + CAPTION_DESCENT,
        "depiction": depiction,
        "mol": mol,
        "alignment": alignment,
        "record": record,
    }


def _align_to_previous(mol: Any, previous: list[dict[str, Any]]) -> tuple[Any | None, dict[str, Any] | None]:
    """Lay out ``mol`` with the part it shares with an earlier structure held at that structure's coordinates.

    The template is whichever of the last few structures shares the most atoms
    (complete rings only). Returns ``(None, None)`` when nothing substantial is
    shared or the constrained layout comes out cramped or uneven, in which case
    the caller lays the structure out from scratch.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDepictor, rdFMCS

    best: tuple[int, dict[str, Any], str] | None = None
    for entry in reversed(previous[-ALIGN_LOOKBACK:]):
        reference = entry["mol"]
        try:
            result = rdFMCS.FindMCS(
                [reference, mol], timeout=ALIGN_TIMEOUT_S, ringMatchesRingOnly=True, completeRingsOnly=True
            )
        except Exception:  # pragma: no cover - defensive
            continue
        smaller = min(reference.GetNumHeavyAtoms(), mol.GetNumHeavyAtoms())
        if result.numAtoms < ALIGN_MIN_ATOMS or result.numAtoms < ALIGN_MIN_FRACTION * smaller:
            continue
        if best is None or result.numAtoms > best[0]:
            best = (result.numAtoms, entry, result.smartsString)
    if best is None:
        return None, None
    _, entry, smarts = best
    # The template is stored at 1.0 units per bond; the constrained layout is
    # built at RDKit's own 1.5, so scale it up first or the core comes out with
    # shorter bonds than the rest of the molecule.
    template = Chem.Mol(entry["mol"])
    conf = template.GetConformer()
    for idx in range(template.GetNumAtoms()):
        point = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, (point.x * RDKIT_DEPICT_BOND, point.y * RDKIT_DEPICT_BOND, 0.0))
    aligned = Chem.Mol(mol)
    aligned.RemoveAllConformers()
    try:
        pairs = rdDepictor.GenerateDepictionMatching2DStructure(aligned, template, -1, Chem.MolFromSmarts(smarts))
    except Exception:
        return None, None
    if not pairs or aligned.GetNumConformers() == 0:
        return None, None
    _normalize_bond_units(aligned)
    if _layout_is_poor(aligned):
        return None, None
    return aligned, {"aligned_to": entry["name"], "aligned_atoms": len(pairs)}


def _layout_is_poor(mol: Any) -> bool:
    """True if non-bonded atoms nearly touch or bond lengths are uneven."""
    conf = mol.GetConformer()
    points = [conf.GetAtomPosition(idx) for idx in range(mol.GetNumAtoms())]
    lengths = [(points[b.GetBeginAtomIdx()] - points[b.GetEndAtomIdx()]).Length() for b in mol.GetBonds()]
    if lengths and max(lengths) > ALIGN_MAX_BOND_RATIO * max(min(lengths), 1e-6):
        return True
    bonded = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()}
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if frozenset((i, j)) not in bonded and (points[i] - points[j]).Length() < ALIGN_MIN_SEPARATION:
                return True
    return False


def _row_width(items: list[dict[str, Any]]) -> int:
    cursor = LEFT_MARGIN
    for item in items:
        cursor += item["width"] + MOL_GAP if item["type"] == "molecule" else item.get("width", 70)
    if items and items[-1]["type"] == "molecule":
        cursor -= MOL_GAP
    return cursor + RIGHT_MARGIN


def _split_row(items: list[dict[str, Any]], max_width: int) -> list[list[dict[str, Any]]]:
    """Split a row before arrows: fewest rows within ``max_width``, then the most even widths.

    Among splits with the fewest rows, the one with the smallest sum of squared
    gaps to ``max_width`` wins, so no row is left with a single step while
    others are full. A stretch between two arrows that is wider than
    ``max_width`` on its own stays whole on one row.
    """
    if _row_width(items) <= max_width:
        return [items]
    cuts = [0] + [i for i in range(1, len(items)) if items[i]["type"] == "arrow"] + [len(items)]
    if len(cuts) == 2:
        return [items]
    best: dict[int, tuple[int, int, int]] = {0: (0, 0, -1)}  # end cut -> (rows, squared gaps, previous cut)
    for j in range(1, len(cuts)):
        for i in range(j):
            if i not in best:
                continue
            width = _row_width(items[cuts[i] : cuts[j]])
            if width > max_width and j - i > 1:
                continue
            rows, cost, _ = best[i]
            candidate = (rows + 1, cost + (max_width - width) ** 2, i)
            if j not in best or candidate[:2] < best[j][:2]:
                best[j] = candidate
    parts: list[list[dict[str, Any]]] = []
    j = len(cuts) - 1
    while j > 0:
        i = best[j][2]
        parts.append(items[cuts[i] : cuts[j]])
        j = i
    return parts[::-1]


def _compose_svg(
    spec: dict[str, Any], rows: list[dict[str, Any]]
) -> tuple[str, list[FigureElement], int, int, list[int]]:
    row_layouts: list[dict[str, Any]] = []
    for row in rows:
        cursor = LEFT_MARGIN
        max_mol_height = MIN_MOL_CELL
        caption_block = LABEL_BLOCK
        laid_out: list[dict[str, Any]] = []
        for item in row["items"]:
            entry = dict(item)
            entry["x"] = cursor
            laid_out.append(entry)
            if item["type"] == "molecule":
                cursor += item["width"] + MOL_GAP
                max_mol_height = max(max_mol_height, item["height"])
                caption_block = max(caption_block, item["caption_height"])
            else:
                cursor += item.get("width", 70)
        row_layouts.append(
            {
                "items": laid_out,
                "max_mol_height": max_mol_height,
                "caption_block": caption_block,
                "width": _row_width(row["items"]),
                "note": row.get("note"),
            }
        )
    content_width = max(row["width"] for row in row_layouts)

    title = str(spec.get("title", "")).strip()
    subtitle = str(spec.get("subtitle", "")).strip()
    title_block = 92 if title or subtitle else TOP_MARGIN
    footnotes = [str(note) for note in spec.get("footnotes", []) if str(note).strip()]
    footer_height = (len(footnotes) * FOOTER_LINE_HEIGHT + 30) if footnotes else 20
    # Titles and footnotes are part of the figure too: widen the canvas rather than clip them.
    text_widths = [_text_width(title, 24, bold=True), _text_width(subtitle, 15)] + [_text_width(n, 13) for n in footnotes]
    text_width = LEFT_MARGIN + max(text_widths) + RIGHT_MARGIN
    width = int(max(int(spec.get("width", 0) or 0), content_width, math.ceil(text_width)))
    total_rows_height = sum(row["max_mol_height"] + row["caption_block"] for row in row_layouts)
    height = int(spec.get("height", 0) or 0) or (
        title_block + total_rows_height + ROW_GAP * (len(row_layouts) - 1) + footer_height
    )

    elements: list[FigureElement] = []
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n',
        "  <style>\n",
        "    .title { font: 700 24px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "    .subtitle { font: 15px Arial, Helvetica, sans-serif; fill: #444; }\n",
        "    .molLabel { font: 700 15px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "    .code { font: 13px Arial, Helvetica, sans-serif; fill: #555; }\n",
        "    .arrowLabel { font: 14px Arial, Helvetica, sans-serif; fill: #222; }\n",
        "    .risk { font: 13px Arial, Helvetica, sans-serif; fill: #444; }\n",
        "    .rxnArrow { stroke: #111; stroke-width: 1.8; fill: none; }\n",
        "    .arrowHead { stroke: #111; stroke-width: 1.8; fill: none; stroke-linecap: round; stroke-linejoin: round; }\n",
        "    .plus { font: 24px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "  </style>\n",
        f'  <rect x="0" y="0" width="{width}" height="{height}" fill="#fff"/>\n',
    ]
    if title:
        parts.append(f'  <text x="{LEFT_MARGIN}" y="42" class="title">{_xml_escape(title)}</text>\n')
        elements.append(
            FigureElement(
                kind="title",
                label="title",
                box=(LEFT_MARGIN, 24, LEFT_MARGIN + _text_width(title, 24, bold=True), 46),
                collide_strictly=False,
            )
        )
    if subtitle:
        parts.append(f'  <text x="{LEFT_MARGIN}" y="68" class="subtitle">{_xml_escape(subtitle)}</text>\n')

    y = title_block
    for row in row_layouts:
        row_mid_y = y + row["max_mol_height"] / 2
        for item in row["items"]:
            if item["type"] == "molecule":
                depiction = item["depiction"]
                cell_x = item["x"] + (item["width"] - depiction.width) / 2
                item_y = y + (row["max_mol_height"] - depiction.height) / 2
                parts.append(f'  <g transform="translate({cell_x:.1f},{item_y:.1f})">\n')
                parts.append(depiction.svg_body)
                parts.append("\n  </g>\n")
                geometry = depiction.translated_geometry(cell_x, item_y)
                ink = geometry["ink_box"] or (cell_x, item_y, cell_x + depiction.width, item_y + depiction.height)
                elements.append(
                    FigureElement(
                        kind="molecule",
                        label=item.get("label") or item.get("id") or "molecule",
                        box=ink,
                        bond_lengths_px=geometry["bond_lengths_px"],
                        label_boxes=geometry["label_boxes"],
                    )
                )
                center = item["x"] + item["width"] / 2
                base_y = y + row["max_mol_height"]
                for line in item["captions"]:
                    line_y = base_y + line["dy"]
                    parts.append(
                        f'  <text x="{center:.1f}" y="{line_y}" class="{line["cls"]}" '
                        f'text-anchor="middle">{_xml_escape(line["text"])}</text>\n'
                    )
                    half = line["width"] / 2
                    kind = "label" if line["cls"] == "molLabel" else "code"
                    elements.append(
                        FigureElement(
                            kind="caption",
                            label=f"{kind}:{line['text']}",
                            box=(center - half, line_y - line["px"] * 0.85, center + half, line_y + line["px"] * 0.2),
                            collide_strictly=False,
                        )
                    )
            elif item["type"] == "arrow":
                x1 = item["x"] + 5
                x2 = item["x"] + item["width"] - 5
                parts.append(_arrow_svg(x1, row_mid_y, x2, item["label"], item["sublabel"]))
                elements.append(
                    FigureElement(
                        kind="arrow",
                        label=f"arrow:{item['label'] or 'step'}",
                        box=(x1, row_mid_y - 34, x2, row_mid_y + 32),
                        collide_strictly=False,
                    )
                )
            elif item["type"] == "plus":
                parts.append(
                    f'  <text x="{item["x"] + item["width"] / 2:.1f}" y="{row_mid_y + 8:.1f}" '
                    f'class="plus" text-anchor="middle">+</text>\n'
                )
            elif item["type"] == "text":
                parts.append(
                    f'  <text x="{item["x"] + item["width"] / 2:.1f}" y="{row_mid_y + 5:.1f}" class="arrowLabel" '
                    f'text-anchor="middle">{_xml_escape(item["text"])}</text>\n'
                )
        y += row["max_mol_height"] + row["caption_block"] + ROW_GAP

    if footnotes:
        footer_y = height - 24 - (len(footnotes) - 1) * FOOTER_LINE_HEIGHT
        for idx, note in enumerate(footnotes):
            parts.append(
                f'  <text x="{LEFT_MARGIN}" y="{footer_y + idx * FOOTER_LINE_HEIGHT}" class="risk">'
                f"{_xml_escape(note)}</text>\n"
            )
    parts.append("</svg>\n")
    return "".join(parts), elements, width, height, [row["width"] for row in row_layouts]


def _arrow_svg(x1: float, y: float, x2: float, label: str = "", sublabel: str = "") -> str:
    mid = (x1 + x2) / 2
    return (
        f'  <path d="M {x1} {y} L {x2} {y}" class="rxnArrow"/>\n'
        f'  <path d="M {x2 - 9} {y - 5.5} L {x2} {y} L {x2 - 9} {y + 5.5}" class="arrowHead"/>\n'
        f'  <text x="{mid}" y="{y - 14}" class="arrowLabel" text-anchor="middle">{_xml_escape(label)}</text>\n'
        f'  <text x="{mid}" y="{y + 24}" class="arrowLabel" text-anchor="middle">{_xml_escape(sublabel)}</text>\n'
    )


def _xml_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
