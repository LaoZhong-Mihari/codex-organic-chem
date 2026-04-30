from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any


FIGURE_TOOL_GUIDE: dict[str, dict[str, Any]] = {
    "chemdraw": {
        "category": "chemically_aware_editor",
        "purpose": "Final manuscript chemical drawing, CDX/CDXML editing, journal style sheets, reaction arrows, and manual mechanism review.",
        "install": ["Commercial Revvity Signals ChemDraw license or institutional access."],
        "inputs": ["CDXML/CDX", "Molfile/SDF", "SMILES by paste/import", "SVG/PDF/EPS export"],
        "notes": [
            "Best final target when available, but no public local automation API is assumed.",
            "codex-chem writes SVG plus ChemDraw-like CDXML objects; open the CDXML/SVG in ChemDraw for manual review and style cleanup.",
        ],
    },
    "chemdoodle": {
        "category": "chemically_aware_editor",
        "purpose": "Chemical drawing with electron-pushing arrows, lone pairs/electrons, atom mapping, CDXML/SMILES/RXN interoperability.",
        "install": ["Commercial ChemDoodle Desktop license, or ChemDoodle Web Components for custom web integrations."],
        "inputs": ["CDXML/CDX", "MOL/RXN", "SMILES", "ChemDoodle JSON"],
        "notes": [
            "Good model for future mechanism-editor integration because electron pushers are chemical objects.",
            "codex-chem exports ChemDoodle JSON pusher shapes that reference atom/bond object ids and electron counts.",
        ],
    },
    "marvin": {
        "category": "chemically_aware_editor",
        "purpose": "ChemAxon reaction/mechanism drawing, electron-flow arrows, atom mapping, and CDX/CDXML import/export.",
        "install": ["Commercial ChemAxon Marvin/Marvin JS access or institutional license."],
        "inputs": ["CDX/CDXML", "MOL/RXN", "SMILES", "CXSMILES/CXON"],
        "notes": [
            "Useful reference implementation for explicit electron-flow arrows and web editor semantics.",
        ],
    },
    "ketcher": {
        "category": "open_source_structure_editor",
        "purpose": "Open-source web editor for molecule/reaction correction before Codex reasoning.",
        "install": [
            "Use the EPAM Ketcher open-source project or a local Ketcher deployment.",
            "Set CODEX_CHEM_KETCHER_URL=http://localhost:<port> when a local Ketcher instance is running.",
            "Optionally set CODEX_CHEM_KETCHER_DIST=/path/to/ketcher/static/build for a local static build.",
        ],
        "inputs": ["SMILES", "Molfile/Rxnfile", "reaction sketches depending on deployment"],
        "notes": [
            "Best open-source choice for the input-review/correction gate.",
            "Used as an open-source correction adapter for structures, reactions, atom maps, and KET/CDXML interop.",
        ],
    },
    "inkscape": {
        "category": "vector_polish",
        "purpose": "Open-source SVG cleanup after chemical content has been reviewed.",
        "install": ["brew install --cask inkscape"],
        "inputs": ["SVG", "PDF", "EPS"],
        "notes": ["Does not validate chemistry; use after chem_mechanism_render publication checks pass."],
    },
    "illustrator": {
        "category": "vector_polish",
        "purpose": "Professional final vector layout/polish after chemical content has been reviewed.",
        "install": ["Commercial Adobe Illustrator license."],
        "inputs": ["SVG", "PDF", "EPS"],
        "notes": ["Does not validate chemistry; use only after mechanism arrows/intermediates are chemically reviewed."],
    },
    "openbabel_cdxml": {
        "category": "format_interop",
        "purpose": "Structure-only CDXML conversion through Open Babel when obabel is installed.",
        "install": ["brew install open-babel"],
        "inputs": ["CDXML for simple chemical structures"],
        "notes": [
            "Open Babel CDXML support is useful for simple structures, but mechanism arrows, lone pairs, and rich graphical annotations may not round-trip.",
        ],
    },
}


ACS_PUBLICATION_STANDARD: dict[str, Any] = {
    "source": "ACS Publications author guidelines and graphics preparation guidance",
    "chemdraw_style": {
        "document_settings": "ACS or ACS-1996 when available",
        "chain_angle_degrees": 120,
        "bond_spacing_percent_width": 18,
        "fixed_bond_length_pt": 14.4,
        "fixed_bond_length_cm": 0.508,
        "fixed_bond_length_in": 0.2,
        "bold_width_pt": 2.0,
        "line_width_pt": 0.6,
        "margin_width_pt": 1.6,
        "hash_spacing_pt": 2.5,
        "font_family": "Arial or Helvetica",
        "font_size_pt": 10,
        "units": "points",
        "tolerance_pixels": 5,
        "page_setup": "US Letter, 100% scale",
    },
    "production_size": {
        "prefer_single_column": True,
        "single_column_width_pt": 240,
        "single_column_width_in": 3.33,
        "single_column_width_cm": 8.47,
        "double_column_max_width_pt": 504,
        "double_column_max_width_in": 7.0,
        "structure_block_one_column_cm": 11.3,
        "structure_block_two_column_cm": 23.6,
        "max_depth_pt": 660,
        "max_depth_in": 9.167,
        "minimum_final_lettering_pt": 4.5,
        "minimum_line_width_pt": 0.5,
    },
    "raster_resolution": {
        "black_white_line_art_dpi": 1200,
        "grayscale_dpi": 600,
        "color_rgb_dpi": 300,
    },
    "deliverables": [
        "Keep a vector master: CDXML/CDX, SVG, PDF, or EPS.",
        "Embed or outline fonts in EPS/PDF.",
        "If raster export is required, export at final production size and required dpi.",
    ],
    "scheme_design": [
        "Use compact chemical schemes, not slide/card layouts.",
        "Avoid rounded panels, shadows, large titles, explanatory paragraphs, and decorative colors.",
        "Use black structures/arrows by default; color only when it clarifies chemistry and survives grayscale.",
        "Put prose in the manuscript text or caption, not inside the artwork.",
        "Show only mechanism-relevant lone pairs and charges.",
        "Curved arrows must be anchored to visible electron sources and sinks.",
    ],
}


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _local_integration_path(name: str, *parts: str) -> Path:
    return _project_root() / "integrations" / name / Path(*parts)


APP_DETECTORS = {
    "chemdraw": {
        "env": "CODEX_CHEM_CHEMDRAW_APP",
        "patterns": ["ChemDraw*.app", "Signals ChemDraw*.app"],
        "open_name": "ChemDraw",
    },
    "chemdoodle": {
        "env": "CODEX_CHEM_CHEMDOODLE_APP",
        "patterns": ["ChemDoodle*.app"],
        "open_name": "ChemDoodle",
    },
    "marvin": {
        "env": "CODEX_CHEM_MARVIN_APP",
        "patterns": ["Marvin*.app", "MarvinSketch*.app"],
        "open_name": "MarvinSketch",
    },
    "inkscape": {
        "env": "CODEX_CHEM_INKSCAPE_APP",
        "patterns": ["Inkscape.app"],
        "open_name": "Inkscape",
        "cli": "inkscape",
    },
    "illustrator": {
        "env": "CODEX_CHEM_ILLUSTRATOR_APP",
        "patterns": ["Adobe Illustrator*.app", "Adobe Illustrator.app"],
        "open_name": "Adobe Illustrator",
    },
}


def _candidate_roots() -> list[Path]:
    return [
        Path("/Applications"),
        Path.home() / "Applications",
        Path("/System/Applications"),
    ]


def _find_apps(patterns: list[str]) -> list[Path]:
    found: list[Path] = []
    for root in _candidate_roots():
        if not root.exists():
            continue
        for pattern in patterns:
            try:
                found.extend(path for path in root.glob(pattern) if path.exists())
                found.extend(path for path in root.glob(f"*/{pattern}") if path.exists())
                found.extend(path for path in root.glob(f"*/*/{pattern}") if path.exists())
            except OSError:
                continue
    unique = {str(path.resolve()): path.resolve() for path in found}
    return [unique[key] for key in sorted(unique)]


def _app_status(tool: str, detector: dict[str, Any]) -> dict[str, Any]:
    guide = FIGURE_TOOL_GUIDE[tool]
    env = detector["env"]
    configured = os.environ.get(env)
    candidates: list[Path] = []
    source = "not_found"
    if configured:
        path = Path(configured).expanduser()
        if path.exists():
            candidates = [path.resolve()]
            source = env
    if not candidates:
        candidates = _find_apps(detector["patterns"])
        if candidates:
            source = "applications_folder"
    cli_name = detector.get("cli")
    cli_path = shutil.which(cli_name) if cli_name else None
    if not candidates and cli_path:
        source = "path"
    status = "available" if candidates or cli_path else "unavailable"
    app_path = str(candidates[0]) if candidates else None
    open_target = app_path or detector.get("open_name")
    return {
        "name": tool,
        "status": status,
        "source": source,
        "app_path": app_path,
        "cli_path": cli_path,
        "env_var": env,
        "category": guide["category"],
        "purpose": guide["purpose"],
        "inputs": guide["inputs"],
        "install": guide["install"],
        "notes": guide["notes"],
        "open_cdxml_command_template": f'open -a "{open_target}" "{{cdxml_path}}"' if open_target else None,
        "open_svg_command_template": f'open -a "{open_target}" "{{svg_path}}"' if open_target else None,
    }


def _ketcher_status() -> dict[str, Any]:
    guide = FIGURE_TOOL_GUIDE["ketcher"]
    url = os.environ.get("CODEX_CHEM_KETCHER_URL")
    dist = os.environ.get("CODEX_CHEM_KETCHER_DIST")
    configured_dist_path = Path(dist).expanduser().resolve() if dist else None
    bundled_dist_path = _local_integration_path("ketcher", "dist").resolve()
    package_path = _local_integration_path("ketcher", "package.json").resolve()
    dist_path = configured_dist_path if configured_dist_path and configured_dist_path.exists() else bundled_dist_path
    available = bool(url) or bool(dist_path and dist_path.exists())
    source = "url" if url else "not_found"
    if dist_path and dist_path.exists():
        source = "CODEX_CHEM_KETCHER_DIST" if configured_dist_path and configured_dist_path.exists() else "bundled_integration_dist"
    elif package_path.exists():
        source = "bundled_integration_source"
    return {
        "name": "ketcher",
        "status": "available" if available else "unavailable",
        "source": source,
        "category": guide["category"],
        "purpose": guide["purpose"],
        "inputs": guide["inputs"],
        "install": guide["install"],
        "notes": guide["notes"],
        "url": url,
        "dist_path": str(dist_path) if dist_path and dist_path.exists() else None,
        "package_path": str(package_path) if package_path.exists() else None,
        "build_command": f"cd {package_path.parent} && npm install && npm run build" if package_path.exists() else None,
        "env_vars": ["CODEX_CHEM_KETCHER_URL", "CODEX_CHEM_KETCHER_DIST"],
    }


def _chemdoodle_web_status() -> dict[str, Any]:
    configured = os.environ.get("CODEX_CHEM_CHEMDOODLE_WEB_DIR")
    configured_path = Path(configured).expanduser().resolve() if configured else None
    bundled_path = _local_integration_path("chemdoodle").resolve()
    path = configured_path if configured_path and configured_path.exists() else bundled_path
    js_path = path / "ChemDoodleWeb.js"
    css_path = path / "ChemDoodleWeb.css"
    available = js_path.exists() and css_path.exists()
    return {
        "name": "chemdoodle_web",
        "status": "available" if available else "unavailable",
        "source": "CODEX_CHEM_CHEMDOODLE_WEB_DIR" if configured_path and configured_path.exists() else ("bundled_integration" if path.exists() else "not_found"),
        "category": "chemically_aware_web_renderer",
        "purpose": "Render ChemDoodle JSON with atom-bound lone pairs and Pusher electron-flow arrows.",
        "inputs": ["ChemDoodle JSON", "SMILES", "MOL/RXN depending on page integration"],
        "install": [
            "Download ChemDoodle Web Components from web.chemdoodle.com and place ChemDoodleWeb.js/ChemDoodleWeb.css in integrations/chemdoodle.",
            "Alternatively set CODEX_CHEM_CHEMDOODLE_WEB_DIR=/path/to/ChemDoodleWeb.",
        ],
        "notes": [
            "ChemDoodle JSON is the preferred object-bound output for lone pairs and electron pushers.",
            "ChemDoodle Web Components are GPLv3 unless a proprietary iChemLabs license is used.",
        ],
        "dist_path": str(path) if path.exists() else None,
        "js_path": str(js_path) if js_path.exists() else None,
        "css_path": str(css_path) if css_path.exists() else None,
        "env_vars": ["CODEX_CHEM_CHEMDOODLE_WEB_DIR"],
    }


def _marvin_js_status() -> dict[str, Any]:
    configured = os.environ.get("CODEX_CHEM_MARVIN_JS_DIR")
    configured_path = Path(configured).expanduser().resolve() if configured else None
    bundled_path = _local_integration_path("marvin").resolve()
    path = configured_path if configured_path and configured_path.exists() else bundled_path
    likely_files = [path / "editor.html", path / "index.html", path / "marvinjslauncher.js"]
    available = any(candidate.exists() for candidate in likely_files)
    return {
        "name": "marvin_js",
        "status": "available" if available else "unavailable",
        "source": "CODEX_CHEM_MARVIN_JS_DIR" if configured_path and configured_path.exists() else ("bundled_integration" if path.exists() else "not_found"),
        "category": "chemically_aware_web_editor",
        "purpose": "Optional ChemAxon Marvin JS review target for CDXML/MOL/RXN/CXSMILES workflows.",
        "inputs": ["CDX/CDXML", "MOL/RXN", "SMILES", "CXSMILES/CXON"],
        "install": [
            "Marvin JS requires a ChemAxon license and package download.",
            "Place the licensed Marvin JS package in integrations/marvin or set CODEX_CHEM_MARVIN_JS_DIR.",
        ],
        "notes": [
            "The local assistant can generate Marvin-compatible CDXML/MOL/RXN-facing artifacts, but cannot bundle Marvin JS without a license.",
        ],
        "dist_path": str(path) if path.exists() else None,
        "env_vars": ["CODEX_CHEM_MARVIN_JS_DIR"],
    }


def _openbabel_cdxml_status() -> dict[str, Any]:
    guide = FIGURE_TOOL_GUIDE["openbabel_cdxml"]
    path = shutil.which("obabel")
    return {
        "name": "openbabel_cdxml",
        "status": "available" if path else "unavailable",
        "category": guide["category"],
        "purpose": guide["purpose"],
        "inputs": guide["inputs"],
        "install": guide["install"],
        "notes": guide["notes"],
        "cli_path": path,
        "example_command": 'obabel input.cdxml -O output.sdf',
    }


def figure_tool_statuses() -> dict[str, Any]:
    tools = {name: _app_status(name, detector) for name, detector in APP_DETECTORS.items()}
    tools["ketcher"] = _ketcher_status()
    tools["chemdoodle_web"] = _chemdoodle_web_status()
    tools["marvin_js"] = _marvin_js_status()
    tools["openbabel_cdxml"] = _openbabel_cdxml_status()
    chemically_aware = [
        name
        for name in ["chemdraw", "chemdoodle", "chemdoodle_web", "marvin", "marvin_js", "ketcher"]
        if tools[name]["status"] == "available"
    ]
    vector_polish = [
        name
        for name in ["inkscape", "illustrator"]
        if tools[name]["status"] == "available"
    ]
    recommendations = [
        "Use chem_mechanism_render to create the checked SVG/spec/CDXML package first.",
        "Prefer ChemDraw, ChemDoodle, or Marvin for final chemical editing when installed.",
        "Use Ketcher mainly for open-source molecule/reaction correction before reasoning.",
        "Use Inkscape or Illustrator only for vector polish after chemical validation.",
    ]
    return {
        "summary": {
            "chemically_aware_editors_available": chemically_aware,
            "vector_polish_tools_available": vector_polish,
            "chemdraw_available": tools["chemdraw"]["status"] == "available",
            "open_source_correction_available": tools["ketcher"]["status"] == "available",
            "object_bound_renderer_available": tools["chemdoodle_web"]["status"] == "available",
            "marvin_js_available": tools["marvin_js"]["status"] == "available",
            "format_interop_available": tools["openbabel_cdxml"]["status"] == "available",
        },
        "tools": tools,
        "guide": FIGURE_TOOL_GUIDE,
        "acs_publication_standard": ACS_PUBLICATION_STANDARD,
        "recommendations": recommendations,
        "boundaries": [
            "ChemDraw/ChemDoodle/Marvin support is file/workflow integration, not unattended GUI automation.",
            "CDXML generated by codex-chem now contains chemical/vector objects, but SVG/spec/trace remain authoritative until human review.",
            "Do not trust a publication mechanism until intermediates, atom maps, arrows, lone pairs, and charges pass publication_checks and receive human chemical review.",
        ],
    }


def figure_followup_for_files(svg_path: str | None = None, cdxml_path: str | None = None) -> dict[str, Any]:
    statuses = figure_tool_statuses()
    tools = statuses["tools"]
    commands: list[dict[str, str]] = []
    if cdxml_path:
        for name in ["chemdraw", "chemdoodle", "marvin"]:
            tool = tools[name]
            template = tool.get("open_cdxml_command_template")
            if tool["status"] == "available" and template:
                commands.append(
                    {
                        "tool": name,
                        "format": "cdxml",
                        "command": template.format(cdxml_path=cdxml_path),
                    }
                )
    if svg_path:
        for name in ["chemdraw", "inkscape", "illustrator"]:
            tool = tools[name]
            template = tool.get("open_svg_command_template")
            if tool["status"] == "available" and template:
                commands.append(
                    {
                        "tool": name,
                        "format": "svg",
                        "command": template.format(svg_path=svg_path),
                    }
                )
    return {
        "available_commands": commands,
        "recommended_sequence": [
            "Open mechanism.cdxml or mechanism.svg in a chemically aware editor such as ChemDraw when available.",
            "Use mechanism.chemdoodle.json for ChemDoodle Web Components integrations and mechanism.trace.json for AI step review.",
            "Verify every intermediate, atom map, lone pair, formal charge, partial charge, and curved-arrow endpoint.",
            "Only after chemical review, polish typography/spacing in ChemDraw, Inkscape, or Illustrator.",
        ],
        "tool_summary": statuses["summary"],
        "boundaries": statuses["boundaries"],
    }
