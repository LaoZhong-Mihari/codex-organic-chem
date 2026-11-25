"""Image cropping, magnification, and coordinate-grid overlays for dense chemical pages.

Vision-capable assistants frequently misread stereo wedges/hashes, ring sizes,
or substituent attachment points when ten structures share one exam page or
journal Scheme. :func:`crop_image` lets an assistant without a shell/Python tool
either:

1. Overlay a labeled coordinate grid (``grid=True``) onto a page image to read
   exact bounding boxes; and
2. Crop and upscale (``scale=2.0``) one or more ``[x1, y1, x2, y2]`` regions
   into standalone PNG files that it can inspect at full resolution or pass to
   OCSR / human review.
"""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw


def _safe_stem(raw: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(raw).strip()).strip("._")
    return cleaned or fallback


def _parse_boxes(boxes: list[dict[str, Any] | list[float]] | str | None) -> list[dict[str, Any]]:
    if boxes is None:
        return []
    if isinstance(boxes, str):
        stripped = boxes.strip()
        if not stripped:
            return []
        parsed = json.loads(stripped)
        if isinstance(parsed, dict):
            parsed = parsed.get("boxes") or parsed.get("crops") or [parsed]
        boxes = parsed
    if not isinstance(boxes, list):
        raise ValueError("boxes must be a list of {id, box: [x1, y1, x2, y2]} objects or coordinate lists.")
    normalized_list: list[dict[str, Any]] = []
    for idx, entry in enumerate(boxes, start=1):
        if isinstance(entry, (list, tuple)) and len(entry) == 4:
            normalized_list.append({"id": f"crop_{idx}", "box": [float(v) for v in entry]})
        elif isinstance(entry, dict):
            coords = entry.get("box") or entry.get("bbox") or entry.get("rect")
            if coords is None and all(k in entry for k in ("x1", "y1", "x2", "y2")):
                coords = [entry["x1"], entry["y1"], entry["x2"], entry["y2"]]
            if not isinstance(coords, (list, tuple)) or len(coords) != 4:
                raise ValueError(f"Crop {idx} must provide 'box': [x1, y1, x2, y2].")
            cid = entry.get("id") or entry.get("label") or entry.get("name") or f"crop_{idx}"
            normalized_list.append({"id": str(cid), "box": [float(v) for v in coords]})
        else:
            raise ValueError(f"Unsupported crop specification at index {idx}: {entry!r}")
    return normalized_list


def _to_pixel_box(
    coords: list[float],
    width: int,
    height: int,
    normalized: bool | None,
    padding_px: int,
) -> tuple[int, int, int, int]:
    x1, y1, x2, y2 = coords
    is_norm = normalized if normalized is not None else (
        0.0 <= min(x1, y1, x2, y2) and max(x1, y1, x2, y2) <= 1.0
    )
    if is_norm:
        px1, py1, px2, py2 = x1 * width, y1 * height, x2 * width, y2 * height
    else:
        px1, py1, px2, py2 = x1, y1, x2, y2
    left = int(round(min(px1, px2))) - max(0, int(padding_px))
    top = int(round(min(py1, py2))) - max(0, int(padding_px))
    right = int(round(max(px1, px2))) + max(0, int(padding_px))
    bottom = int(round(max(py1, py2))) + max(0, int(padding_px))
    left = max(0, min(left, width - 1))
    top = max(0, min(top, height - 1))
    right = max(left + 1, min(right, width))
    bottom = max(top + 1, min(bottom, height))
    return left, top, right, bottom


def _draw_grid_overlay(img: Image.Image, steps: int = 10) -> Image.Image:
    base = img.convert("RGBA")
    overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = base.size

    for i in range(1, steps):
        frac = i / steps
        x = int(round(frac * w))
        y = int(round(frac * h))
        major = (i % 5 == 0)
        line_color = (220, 38, 38, 120) if major else (37, 99, 235, 80)
        line_w = 2 if major else 1
        draw.line([(x, 0), (x, h)], fill=line_color, width=line_w)
        draw.line([(0, y), (w, y)], fill=line_color, width=line_w)

        x_label = f"{frac:.1f} ({x})"
        y_label = f"{frac:.1f} ({y})"
        draw.rectangle((x + 2, 2, x + 6 + len(x_label) * 6, 16), fill=(255, 255, 255, 215))
        draw.text((x + 4, 3), x_label, fill=(30, 41, 59, 255))
        draw.rectangle((2, y + 2, 6 + len(y_label) * 6, y + 16), fill=(255, 255, 255, 215))
        draw.text((4, y + 3), y_label, fill=(30, 41, 59, 255))

    return Image.alpha_composite(base, overlay).convert("RGB")


def crop_image(
    image_path: str | Path,
    boxes: list[dict[str, Any] | list[float]] | str | None = None,
    output_dir: str | Path | None = None,
    *,
    grid: bool = False,
    scale: float = 2.0,
    padding_px: int = 12,
    normalized: bool | None = None,
) -> dict[str, Any]:
    """Crop and magnify regions of a page/scheme image, or render a coordinate grid.

    Parameters
    ----------
    image_path:
        Source PNG/JPEG/WebP page or scheme image.
    boxes:
        List of ``{"id": "A", "box": [x1, y1, x2, y2]}`` (or JSON string).
        Coordinates in ``[0.0, 1.0]`` are treated as normalized fractions of
        ``(width, height)`` unless ``normalized=False``; larger values are
        treated as pixels.
    output_dir:
        Directory where cropped PNGs / grid image are written. Defaults to a
        persistent temp directory when omitted.
    grid:
        When ``True``, also writes ``grid_overlay.png`` with 0.1-step normalized
        and pixel coordinate rulers along both axes.
    scale:
        Upscaling factor applied to each crop (default ``2.0``, clamped to
        ``[1.0, 4.0]``) so small stereo wedges and atom labels become legible.
    padding_px:
        Extra pixel margin retained around each requested box before upscaling.
    """
    src = Path(str(image_path)).expanduser()
    if not src.is_file():
        return {
            "status": "error",
            "path": str(src),
            "warnings": [f"Image file not found: {src}"],
            "crops": [],
        }

    try:
        parsed_boxes = _parse_boxes(boxes)
    except (ValueError, json.JSONDecodeError) as exc:
        return {
            "status": "error",
            "path": str(src),
            "warnings": [f"Invalid boxes specification: {exc}"],
            "crops": [],
        }

    if not parsed_boxes and not grid:
        grid = True  # Calling with only image_path produces a coordinate ruler + dimensions.

    out_dir = Path(output_dir).expanduser().resolve() if output_dir else Path(
        tempfile.mkdtemp(prefix="codex-chem-crops-")
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    with Image.open(src) as raw_img:
        img = raw_img.convert("RGB")
        width, height = img.size
        grid_file: str | None = None
        if grid:
            grid_img = _draw_grid_overlay(img)
            grid_dest = out_dir / "grid_overlay.png"
            grid_img.save(grid_dest, format="PNG")
            grid_file = str(grid_dest)

        bounded_scale = max(1.0, min(float(scale), 4.0))
        crop_records: list[dict[str, Any]] = []
        used_names: set[str] = set()
        for idx, spec in enumerate(parsed_boxes, start=1):
            px_box = _to_pixel_box(spec["box"], width, height, normalized, padding_px)
            sub = img.crop(px_box)
            if bounded_scale > 1.01:
                new_size = (
                    max(1, int(round(sub.width * bounded_scale))),
                    max(1, int(round(sub.height * bounded_scale))),
                )
                sub = sub.resize(new_size, Image.Resampling.LANCZOS)
            stem = _safe_stem(spec["id"], f"crop_{idx}")
            while stem in used_names:
                stem = f"{stem}_{idx}"
            used_names.add(stem)
            dest = out_dir / f"{stem}.png"
            sub.save(dest, format="PNG")
            crop_records.append(
                {
                    "id": spec["id"],
                    "path": str(dest),
                    "pixel_box": list(px_box),
                    "normalized_box": [
                        round(px_box[0] / width, 4),
                        round(px_box[1] / height, 4),
                        round(px_box[2] / width, 4),
                        round(px_box[3] / height, 4),
                    ],
                    "output_size": [sub.width, sub.height],
                }
            )

    return {
        "status": "ok",
        "source_image": str(src.resolve()),
        "source_size": {"width": width, "height": height},
        "output_dir": str(out_dir),
        "grid_image": grid_file,
        "crops": crop_records,
        "warnings": [],
    }
