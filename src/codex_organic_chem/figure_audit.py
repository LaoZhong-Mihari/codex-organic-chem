"""Automatic quality audit for rendered chemistry figures.

Every figure this package returns is measured before it reaches the user. The
audit is geometric, not stylistic: it reads the placed elements' real ink boxes
and the drawn bond lengths, so it catches the failures a renderer cannot see in
its own output -- structures overlapping each other, ink running past the canvas
edge, captions sitting on top of a molecule, and bond lengths drifting between
panels of one figure.

Callers do not opt in. ``route_figure`` and ``mechanism_canvas`` both run
:func:`audit_figure` and surface the result as ``publication_checks``, and a
blocking issue downgrades the reported status so a broken figure is never
described as publication-ready.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any

Box = tuple[float, float, float, float]


@dataclass
class FigureElement:
    """One placed thing in a composed figure, in final canvas coordinates."""

    kind: str
    label: str
    box: Box
    bond_lengths_px: list[float] = field(default_factory=list)
    label_boxes: list[Box] = field(default_factory=list)
    #: Text and arrows may legitimately sit close to structures; only
    #: structure-vs-structure contact is treated as a hard failure.
    collide_strictly: bool = True


def _intersection_area(a: Box, b: Box) -> float:
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if dx > 0 and dy > 0 else 0.0


def _area(box: Box) -> float:
    return max(box[2] - box[0], 0.0) * max(box[3] - box[1], 0.0)


def _overlap_checks(elements: list[FigureElement]) -> tuple[list[str], list[str], list[dict[str, Any]]]:
    blocking: list[str] = []
    warnings: list[str] = []
    details: list[dict[str, Any]] = []
    for i, first in enumerate(elements):
        for second in elements[i + 1 :]:
            overlap = _intersection_area(first.box, second.box)
            if overlap <= 0.0:
                continue
            smaller = min(_area(first.box), _area(second.box)) or 1.0
            fraction = overlap / smaller
            entry = {
                "a": first.label or first.kind,
                "b": second.label or second.kind,
                "kinds": [first.kind, second.kind],
                "overlap_fraction": round(fraction, 4),
            }
            details.append(entry)
            strict = first.collide_strictly and second.collide_strictly
            message = (
                f"'{entry['a']}' and '{entry['b']}' overlap by "
                f"{fraction * 100:.0f}% of the smaller element."
            )
            if strict and fraction > 0.02:
                blocking.append(message)
            elif fraction > 0.25:
                warnings.append(message)
    return blocking, warnings, details


def _clipping_checks(elements: list[FigureElement], width: float, height: float) -> list[str]:
    issues: list[str] = []
    for element in elements:
        x0, y0, x1, y1 = element.box
        if x0 < -0.5 or y0 < -0.5 or x1 > width + 0.5 or y1 > height + 0.5:
            issues.append(
                f"'{element.label or element.kind}' extends outside the "
                f"{width:.0f}x{height:.0f} canvas and will be clipped."
            )
    return issues


def _bond_length_checks(
    elements: list[FigureElement], target_px: float | None
) -> tuple[list[str], dict[str, Any]]:
    lengths: list[float] = []
    per_element: list[dict[str, Any]] = []
    for element in elements:
        if not element.bond_lengths_px:
            continue
        element_median = median(element.bond_lengths_px)
        lengths.append(element_median)
        per_element.append({"label": element.label or element.kind, "median_bond_px": round(element_median, 2)})
    summary: dict[str, Any] = {"per_element": per_element}
    if not lengths:
        return [], summary
    low, high = min(lengths), max(lengths)
    summary["min_median_bond_px"] = round(low, 2)
    summary["max_median_bond_px"] = round(high, 2)
    if target_px:
        summary["target_bond_px"] = round(float(target_px), 2)
    warnings: list[str] = []
    # A figure reads as one drawing only if every structure shares a bond length.
    if low > 0 and high / low > 1.12:
        warnings.append(
            f"Bond lengths are inconsistent across the figure "
            f"({low:.1f}px to {high:.1f}px); structures will not look like one drawing."
        )
    if target_px and low > 0 and low < float(target_px) * 0.9:
        warnings.append(
            f"At least one structure was drawn below the target bond length "
            f"({low:.1f}px vs {float(target_px):.1f}px); its cell is too small."
        )
    return warnings, summary


def _label_collision_checks(elements: list[FigureElement]) -> list[str]:
    """Atom labels colliding inside one structure make it chemically unreadable."""
    warnings: list[str] = []
    for element in elements:
        boxes = element.label_boxes
        hits = 0
        for i, first in enumerate(boxes):
            for second in boxes[i + 1 :]:
                if _intersection_area(first, second) > 0.35 * (min(_area(first), _area(second)) or 1.0):
                    hits += 1
        if hits:
            warnings.append(
                f"'{element.label or element.kind}' has {hits} overlapping atom "
                "label(s); the depiction may be unreadable."
            )
    return warnings


def svg_to_png(svg_path: Path, png_path: Path, *, scale: float = 1.0) -> list[str]:
    """Rasterize an SVG with whichever converter is installed."""
    commands = [
        ["rsvg-convert", "-z", str(scale), "-o", str(png_path), str(svg_path)],
        ["magick", "-density", str(96 * scale), str(svg_path), str(png_path)],
        ["convert", "-density", str(96 * scale), str(svg_path), str(png_path)],
        ["inkscape", str(svg_path), "--export-type=png", f"--export-dpi={96 * scale}", f"--export-filename={png_path}"],
    ]
    chrome = _chrome_path()
    if chrome:
        width, height = svg_dimensions(svg_path)
        commands.append(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-first-run",
                "--no-default-browser-check",
                "--hide-scrollbars",
                f"--force-device-scale-factor={scale}",
                "--default-background-color=FFFFFFFF",
                f"--window-size={width},{height}",
                f"--screenshot={png_path}",
                svg_path.as_uri(),
            ]
        )
    for cmd in commands:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=60)
        except FileNotFoundError:
            continue
        except Exception:
            continue
        if png_path.exists() and png_path.stat().st_size > 0:
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
        if Path(candidate).exists():
            return candidate
        found = shutil.which(candidate)
        if found:
            return found
    return None


def svg_dimensions(svg_path: Path) -> tuple[int, int]:
    text = svg_path.read_text(encoding="utf-8", errors="ignore")
    width = re.search(r"\bwidth=['\"](\d+)", text)
    height = re.search(r"\bheight=['\"](\d+)", text)
    return (int(width.group(1)) if width else 1200, int(height.group(1)) if height else 800)


def raster_ink_check(png_path: Path) -> tuple[list[str], dict[str, Any]]:
    """Confirm the rasterized figure actually contains drawn content.

    This is the one check that sees what the user will see. It catches a blank or
    near-blank export, which a purely geometric audit cannot detect.
    """
    try:
        from PIL import Image
    except Exception as exc:  # pragma: no cover - optional dependency
        return [], {"status": "skipped", "reason": f"Pillow unavailable: {exc}"}
    try:
        with Image.open(png_path) as image:
            gray = image.convert("L")
            histogram = gray.histogram()
            total = sum(histogram) or 1
            ink = sum(histogram[:200])
            size = gray.size
    except Exception as exc:
        return [], {"status": "error", "reason": str(exc)}
    fraction = ink / total
    summary = {
        "status": "ok",
        "pixels": total,
        "size": list(size),
        "ink_fraction": round(fraction, 5),
    }
    warnings: list[str] = []
    if fraction < 0.0005:
        warnings.append(
            f"The rasterized figure is effectively blank (ink covers {fraction * 100:.3f}% "
            "of the image); the export did not capture the drawing."
        )
    elif fraction > 0.55:
        warnings.append(
            f"The rasterized figure is unusually dense (ink covers {fraction * 100:.0f}% "
            "of the image); check for overlapping or oversized content."
        )
    return warnings, summary


def audit_figure(
    *,
    renderer: str,
    width: float,
    height: float,
    elements: list[FigureElement],
    target_bond_px: float | None = None,
    png_path: Path | None = None,
    extra_blocking: list[str] | None = None,
    extra_warnings: list[str] | None = None,
) -> dict[str, Any]:
    """Measure a rendered figure and report blocking issues plus warnings."""
    blocking: list[str] = list(extra_blocking or [])
    warnings: list[str] = list(extra_warnings or [])

    overlap_blocking, overlap_warnings, overlap_details = _overlap_checks(elements)
    blocking.extend(overlap_blocking)
    warnings.extend(overlap_warnings)
    blocking.extend(_clipping_checks(elements, width, height))

    bond_warnings, bond_summary = _bond_length_checks(elements, target_bond_px)
    warnings.extend(bond_warnings)
    warnings.extend(_label_collision_checks(elements))

    if not elements:
        blocking.append("The figure contains no measurable drawn elements.")

    raster_summary: dict[str, Any] = {"status": "not_requested"}
    if png_path is not None and png_path.exists():
        raster_warnings, raster_summary = raster_ink_check(png_path)
        warnings.extend(raster_warnings)

    return {
        "renderer": renderer,
        "audit": "automatic_geometry_and_raster",
        "vector_master": True,
        "canvas": [round(float(width), 1), round(float(height), 1)],
        "element_count": len(elements),
        "bond_geometry": bond_summary,
        "overlaps": overlap_details,
        "raster_check": raster_summary,
        "blocking_issues": blocking,
        "warnings": warnings,
        "ready_for_visual_review": not blocking,
        "publication_ready": False,
    }
