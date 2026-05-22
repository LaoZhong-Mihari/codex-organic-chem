from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .rdkit_tools import RDKIT_AVAILABLE, normalize_structure

try:  # pragma: no cover - import availability is environment dependent
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception:  # pragma: no cover
    Chem = None
    Draw = None
    rdDepictor = None
    rdMolDraw2D = None


DEFAULT_MOL_WIDTH = 345
DEFAULT_MOL_HEIGHT = 300
ROW_GAP = 110
TOP_MARGIN = 115
LEFT_MARGIN = 40
RIGHT_MARGIN = 40
LABEL_BLOCK = 58
FOOTER_LINE_HEIGHT = 27


def route_figure_spec_example() -> dict[str, Any]:
    return {
        "title": "Example protected route",
        "subtitle": "Composed route figure from validated molecule SMILES.",
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

    warnings: list[str] = []
    normalized: list[dict[str, Any]] = []
    try:
        prepared_rows = [_prepare_row(row, normalized, warnings) for row in rows]
        svg = _compose_svg(spec, prepared_rows)
    except ValueError as exc:
        return {"status": "error", "warnings": [str(exc), *warnings]}

    output_files: list[str] = []
    if output_dir:
        out_dir = Path(output_dir).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        if "svg" in formats or "png" in formats:
            svg_path = out_dir / f"{basename}.svg"
            svg_path.write_text(svg, encoding="utf-8")
            output_files.append(str(svg_path))
        if "png" in formats:
            png_path = out_dir / f"{basename}.png"
            raster_warnings = _svg_to_png(svg_path, png_path)
            warnings.extend(raster_warnings)
            if png_path.exists():
                output_files.append(str(png_path))

    status = "ok"
    if "png" in formats and output_dir and not any(path.endswith(".png") for path in output_files):
        status = "ok_with_warnings"

    return {
        "status": status,
        "kind": "route_figure",
        "formats": formats,
        "svg": svg,
        "output_files": output_files,
        "molecules": normalized,
        "publication_checks": {
            "renderer": "codex_route_figure_composed_svg",
            "vector_master": True,
            "ready_for_visual_review": True,
            "blocking_issues": [],
        },
        "warnings": warnings,
    }


def _prepare_row(row: dict[str, Any], normalized: list[dict[str, Any]], warnings: list[str]) -> dict[str, Any]:
    items = row.get("items")
    if not isinstance(items, list) or not items:
        raise ValueError("Each route figure row requires a non-empty 'items' list.")
    prepared_items: list[dict[str, Any]] = []
    for raw in items:
        if not isinstance(raw, dict):
            raise ValueError("Route figure row items must be objects.")
        item_type = raw.get("type") or ("molecule" if raw.get("smiles") or raw.get("molfile") else None)
        if item_type == "molecule":
            prepared = _prepare_molecule(raw)
            prepared_items.append(prepared)
            normalized.append(
                {
                    "id": prepared.get("id"),
                    "label": prepared.get("label"),
                    "code": prepared.get("code"),
                    "input_smiles": raw.get("smiles"),
                    "canonical_smiles": prepared["record"]["canonical_smiles"],
                    "isomeric_smiles": prepared["record"]["isomeric_smiles"],
                    "inchi_key": prepared["record"].get("inchi_key"),
                    "warnings": prepared["record"].get("warnings", []),
                    "metadata": prepared["record"].get("metadata", {}),
                }
            )
            warnings.extend(prepared["record"].get("warnings", []))
        elif item_type == "arrow":
            prepared_items.append(
                {
                    "type": "arrow",
                    "label": str(raw.get("label", "")),
                    "sublabel": str(raw.get("sublabel") or raw.get("label2") or ""),
                    "width": int(raw.get("width", 125)),
                }
            )
        elif item_type == "plus":
            prepared_items.append({"type": "plus", "width": int(raw.get("width", 55))})
        elif item_type == "text":
            prepared_items.append({"type": "text", "text": str(raw.get("text", "")), "width": int(raw.get("width", 120))})
        else:
            raise ValueError("Route figure item requires type 'molecule', 'arrow', 'plus', or 'text'.")
    return {"items": prepared_items, "note": row.get("note")}


def _prepare_molecule(raw: dict[str, Any]) -> dict[str, Any]:
    smiles = raw.get("smiles")
    molfile = raw.get("molfile")
    if not smiles and not molfile:
        raise ValueError("Molecule route figure item requires 'smiles' or 'molfile'.")
    record = normalize_structure(smiles=smiles, molfile=molfile, source="route_figure").to_dict()
    if not record.get("canonical_smiles") and not record.get("molfile"):
        raise ValueError(f"Invalid molecule for route figure item {raw.get('id') or raw.get('label') or '<unnamed>'}.")
    mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False) if molfile else Chem.MolFromSmiles(record["isomeric_smiles"])
    if mol is None:
        raise ValueError(f"RDKit could not prepare molecule for route figure item {raw.get('id') or raw.get('label') or '<unnamed>'}.")
    width = int(raw.get("width", DEFAULT_MOL_WIDTH))
    height = int(raw.get("height", DEFAULT_MOL_HEIGHT))
    return {
        "type": "molecule",
        "id": str(raw.get("id") or raw.get("code") or raw.get("label") or ""),
        "label": str(raw.get("label") or raw.get("title") or raw.get("id") or ""),
        "code": str(raw.get("code") or raw.get("id") or ""),
        "width": width,
        "height": height,
        "svg": _molecule_svg_body(mol, width, height),
        "record": record,
    }


def _molecule_svg_body(mol: Any, width: int, height: int) -> str:
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()
    opts.bondLineWidth = 1.65
    opts.minFontSize = 9
    opts.maxFontSize = 16
    opts.padding = 0.07
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    body = re.sub(r"^.*?<!-- END OF HEADER -->", "", svg, flags=re.S)
    body = re.sub(r"<\?xml[^>]*>\s*", "", body).strip()
    body = re.sub(r"<svg[^>]*>", "", body, count=1).strip()
    body = re.sub(r"<rect[^>]*>\s*</rect>", "", body)
    body = re.sub(r"</svg>\s*$", "", body).strip()
    return body


def _compose_svg(spec: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    row_layouts: list[dict[str, Any]] = []
    content_width = 0
    for row in rows:
        cursor = LEFT_MARGIN
        max_mol_height = DEFAULT_MOL_HEIGHT
        laid_out: list[dict[str, Any]] = []
        for item in row["items"]:
            entry = dict(item)
            entry["x"] = cursor
            laid_out.append(entry)
            if item["type"] == "molecule":
                cursor += item["width"] + 25
                max_mol_height = max(max_mol_height, item["height"])
            else:
                cursor += item.get("width", 70)
        width = cursor + RIGHT_MARGIN
        content_width = max(content_width, width)
        row_layouts.append({"items": laid_out, "max_mol_height": max_mol_height, "width": width, "note": row.get("note")})

    title = str(spec.get("title", "")).strip()
    subtitle = str(spec.get("subtitle", "")).strip()
    title_block = 92 if title or subtitle else 35
    footnotes = [str(note) for note in spec.get("footnotes", []) if str(note).strip()]
    footer_height = (len(footnotes) * FOOTER_LINE_HEIGHT + 35) if footnotes else 25
    width = max(int(spec.get("width", 0) or 0), content_width, 900)
    total_rows_height = sum(row["max_mol_height"] + LABEL_BLOCK for row in row_layouts)
    height = int(spec.get("height", 0) or 0) or title_block + total_rows_height + ROW_GAP * (len(row_layouts) - 1) + footer_height

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n',
        "  <style>\n",
        "    .title { font: 700 26px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "    .subtitle { font: 16px Arial, Helvetica, sans-serif; fill: #444; }\n",
        "    .molLabel { font: 700 16px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "    .code { font: 13px Arial, Helvetica, sans-serif; fill: #555; }\n",
        "    .arrowLabel { font: 14px Arial, Helvetica, sans-serif; fill: #222; }\n",
        "    .risk { font: 14px Arial, Helvetica, sans-serif; fill: #444; }\n",
        "    .rxnArrow { stroke: #111; stroke-width: 2; fill: none; }\n",
        "    .arrowHead { stroke: #111; stroke-width: 2; fill: none; stroke-linecap: round; stroke-linejoin: round; }\n",
        "    .plus { font: 26px Arial, Helvetica, sans-serif; fill: #111; }\n",
        "  </style>\n",
        f'  <rect x="0" y="0" width="{width}" height="{height}" fill="#fff"/>\n',
    ]
    if title:
        parts.append(f'  <text x="40" y="44" class="title">{_xml_escape(title)}</text>\n')
    if subtitle:
        parts.append(f'  <text x="40" y="72" class="subtitle">{_xml_escape(subtitle)}</text>\n')

    y = TOP_MARGIN if title or subtitle else 40
    for row in row_layouts:
        row_mid_y = y + row["max_mol_height"] / 2
        for item in row["items"]:
            if item["type"] == "molecule":
                item_y = y + (row["max_mol_height"] - item["height"]) / 2
                parts.append(f'  <g transform="translate({item["x"]},{item_y})">\n')
                parts.append(item["svg"])
                parts.append("\n  </g>\n")
                center = item["x"] + item["width"] / 2
                if item.get("label"):
                    parts.append(
                        f'  <text x="{center}" y="{y + row["max_mol_height"] + 28}" class="molLabel" '
                        f'text-anchor="middle">{_xml_escape(item["label"])}</text>\n'
                    )
                if item.get("code"):
                    parts.append(
                        f'  <text x="{center}" y="{y + row["max_mol_height"] + 48}" class="code" '
                        f'text-anchor="middle">{_xml_escape(item["code"])}</text>\n'
                    )
            elif item["type"] == "arrow":
                x1 = item["x"] + 5
                x2 = item["x"] + item["width"] - 5
                parts.append(_arrow_svg(x1, row_mid_y, x2, item["label"], item["sublabel"]))
            elif item["type"] == "plus":
                parts.append(f'  <text x="{item["x"] + item["width"] / 2}" y="{row_mid_y + 9}" class="plus" text-anchor="middle">+</text>\n')
            elif item["type"] == "text":
                parts.append(
                    f'  <text x="{item["x"] + item["width"] / 2}" y="{row_mid_y + 5}" class="arrowLabel" '
                    f'text-anchor="middle">{_xml_escape(item["text"])}</text>\n'
                )
        y += row["max_mol_height"] + LABEL_BLOCK + ROW_GAP

    if footnotes:
        footer_y = height - 28 - (len(footnotes) - 1) * FOOTER_LINE_HEIGHT
        for idx, note in enumerate(footnotes):
            parts.append(f'  <text x="40" y="{footer_y + idx * FOOTER_LINE_HEIGHT}" class="risk">{_xml_escape(note)}</text>\n')
    parts.append("</svg>\n")
    return "".join(parts)


def _arrow_svg(x1: float, y: float, x2: float, label: str = "", sublabel: str = "") -> str:
    mid = (x1 + x2) / 2
    return (
        f'  <path d="M {x1} {y} L {x2} {y}" class="rxnArrow"/>\n'
        f'  <path d="M {x2 - 10} {y - 6} L {x2} {y} L {x2 - 10} {y + 6}" class="arrowHead"/>\n'
        f'  <text x="{mid}" y="{y - 20}" class="arrowLabel" text-anchor="middle">{_xml_escape(label)}</text>\n'
        f'  <text x="{mid}" y="{y + 28}" class="arrowLabel" text-anchor="middle">{_xml_escape(sublabel)}</text>\n'
    )


def _svg_to_png(svg_path: Path, png_path: Path) -> list[str]:
    commands = [
        ["rsvg-convert", "-o", str(png_path), str(svg_path)],
        ["magick", str(svg_path), str(png_path)],
        ["convert", str(svg_path), str(png_path)],
        ["inkscape", str(svg_path), "--export-type=png", f"--export-filename={png_path}"],
    ]
    chrome = _chrome_path()
    if chrome:
        width, height = _svg_dimensions(svg_path)
        commands.append(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-first-run",
                "--no-default-browser-check",
                f"--window-size={width},{height}",
                f"--screenshot={png_path}",
                svg_path.as_uri(),
            ]
        )
    for cmd in commands:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=30)
        except FileNotFoundError:
            continue
        except Exception:
            continue
        if png_path.exists():
            return []
    return [
        "PNG rasterization was requested, but no SVG rasterizer was available. "
        "Install rsvg-convert, ImageMagick, Inkscape, or Chrome, or use the SVG master."
    ]


def _chrome_path() -> str | None:
    candidates = [
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        "google-chrome",
        "chromium",
        "chrome",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return str(path)
        found = shutil.which(candidate)
        if found:
            return found
    return None


def _svg_dimensions(svg_path: Path) -> tuple[int, int]:
    text = svg_path.read_text(encoding="utf-8", errors="ignore")
    width_match = re.search(r"\bwidth=\"(\d+)\"", text)
    height_match = re.search(r"\bheight=\"(\d+)\"", text)
    width = int(width_match.group(1)) if width_match else 1200
    height = int(height_match.group(1)) if height_match else 800
    return width, height


def _xml_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
