"""Composed route/scheme figure renderer.

Molecules are drawn through :mod:`figure_style`, so every structure in a route
shares one bond length, line width, and font size regardless of its size — the
property that makes a multi-step scheme read as a single drawing. Cells are
sized to each molecule's ink instead of a fixed 345x300 box, and rows are packed
from those real sizes, so captions and neighbouring structures cannot land on
top of one another by construction.

Every render finishes with :func:`figure_audit.audit_figure`; its result is
returned as ``publication_checks`` and any blocking issue downgrades ``status``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .figure_audit import FigureElement, audit_figure, svg_to_png
from .figure_style import render_molecule, resolve_preset
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
    warnings: list[str] = []
    normalized: list[dict[str, Any]] = []
    try:
        prepared_rows = [_prepare_row(row, normalized, warnings, preset) for row in rows]
        svg, elements, width, height = _compose_svg(spec, prepared_rows)
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
        "publication_checks": publication_checks,
        "warnings": warnings,
    }


def _prepare_row(
    row: dict[str, Any],
    normalized: list[dict[str, Any]],
    warnings: list[str],
    preset: dict[str, Any],
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
            prepared = _prepare_molecule(raw, preset)
            prepared_items.append(prepared)
            record = prepared["record"]
            normalized.append(
                {
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
            )
            warnings.extend(record.get("warnings", []))
            warnings.extend(prepared["depiction"].warnings)
        elif item_type == "arrow":
            label = str(raw.get("label", ""))
            sublabel = str(raw.get("sublabel") or raw.get("label2") or "")
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


def _prepare_molecule(raw: dict[str, Any], preset: dict[str, Any]) -> dict[str, Any]:
    smiles = raw.get("smiles")
    molfile = raw.get("molfile")
    if not smiles and not molfile:
        raise ValueError("Molecule route figure item requires 'smiles' or 'molfile'.")
    record = normalize_structure(smiles=smiles, molfile=molfile, source="route_figure").to_dict()
    if not record.get("canonical_smiles") and not record.get("molfile"):
        raise ValueError(f"Invalid molecule for route figure item {raw.get('id') or raw.get('label') or '<unnamed>'}.")
    try:
        depiction = render_molecule(
            smiles=record.get("isomeric_smiles") or smiles,
            molfile=molfile,
            preset=preset,
        )
    except ValueError as exc:
        raise ValueError(
            f"RDKit could not draw route figure item {raw.get('id') or raw.get('label') or '<unnamed>'}: {exc}"
        ) from exc
    return {
        "type": "molecule",
        "id": str(raw.get("id") or raw.get("code") or raw.get("label") or ""),
        "label": str(raw.get("label") or raw.get("title") or raw.get("id") or ""),
        "code": str(raw.get("code") or raw.get("id") or ""),
        "width": max(depiction.width, MIN_MOL_CELL),
        "height": depiction.height,
        "depiction": depiction,
        "record": record,
    }


def _compose_svg(
    spec: dict[str, Any], rows: list[dict[str, Any]]
) -> tuple[str, list[FigureElement], int, int]:
    row_layouts: list[dict[str, Any]] = []
    content_width = 0
    for row in rows:
        cursor = LEFT_MARGIN
        max_mol_height = MIN_MOL_CELL
        laid_out: list[dict[str, Any]] = []
        for item in row["items"]:
            entry = dict(item)
            entry["x"] = cursor
            laid_out.append(entry)
            if item["type"] == "molecule":
                cursor += item["width"] + MOL_GAP
                max_mol_height = max(max_mol_height, item["height"])
            else:
                cursor += item.get("width", 70)
        width = cursor - (MOL_GAP if laid_out and laid_out[-1]["type"] == "molecule" else 0) + RIGHT_MARGIN
        content_width = max(content_width, width)
        row_layouts.append({"items": laid_out, "max_mol_height": max_mol_height, "note": row.get("note")})

    title = str(spec.get("title", "")).strip()
    subtitle = str(spec.get("subtitle", "")).strip()
    title_block = 92 if title or subtitle else TOP_MARGIN
    footnotes = [str(note) for note in spec.get("footnotes", []) if str(note).strip()]
    footer_height = (len(footnotes) * FOOTER_LINE_HEIGHT + 30) if footnotes else 20
    width = max(int(spec.get("width", 0) or 0), content_width)
    total_rows_height = sum(row["max_mol_height"] + LABEL_BLOCK for row in row_layouts)
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
                box=(LEFT_MARGIN, 24, LEFT_MARGIN + len(title) * 24 * CHAR_W, 46),
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
                caption_y = y + row["max_mol_height"] + 24
                if item.get("label"):
                    parts.append(
                        f'  <text x="{center:.1f}" y="{caption_y}" class="molLabel" '
                        f'text-anchor="middle">{_xml_escape(item["label"])}</text>\n'
                    )
                    half = len(item["label"]) * 15 * CHAR_W / 2
                    elements.append(
                        FigureElement(
                            kind="caption",
                            label=f"label:{item['label']}",
                            box=(center - half, caption_y - 13, center + half, caption_y + 3),
                            collide_strictly=False,
                        )
                    )
                if item.get("code"):
                    parts.append(
                        f'  <text x="{center:.1f}" y="{caption_y + 19}" class="code" '
                        f'text-anchor="middle">{_xml_escape(item["code"])}</text>\n'
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
        y += row["max_mol_height"] + LABEL_BLOCK + ROW_GAP

    if footnotes:
        footer_y = height - 24 - (len(footnotes) - 1) * FOOTER_LINE_HEIGHT
        for idx, note in enumerate(footnotes):
            parts.append(
                f'  <text x="{LEFT_MARGIN}" y="{footer_y + idx * FOOTER_LINE_HEIGHT}" class="risk">'
                f"{_xml_escape(note)}</text>\n"
            )
    parts.append("</svg>\n")
    return "".join(parts), elements, width, height


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
