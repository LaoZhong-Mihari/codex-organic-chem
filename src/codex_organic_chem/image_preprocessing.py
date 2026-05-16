from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import os
import shutil
import tempfile
from typing import Any


@dataclass(frozen=True)
class OcsrImageVariant:
    name: str
    path: Path
    temporary: bool = False


@dataclass
class OcsrImageVariantSet:
    variants: list[OcsrImageVariant]
    warnings: list[str] = field(default_factory=list)
    cleanup_dir: Path | None = None

    def cleanup(self) -> None:
        if self.cleanup_dir and os.environ.get("CODEX_CHEM_KEEP_PREPROCESSED_IMAGES") != "1":
            shutil.rmtree(self.cleanup_dir, ignore_errors=True)


def build_ocsr_image_variants(image_path: Path) -> OcsrImageVariantSet:
    variants = [OcsrImageVariant("original", image_path, temporary=False)]
    warnings: list[str] = []
    if os.environ.get("CODEX_CHEM_DISABLE_IMAGE_PREPROCESSING") == "1":
        return OcsrImageVariantSet(variants=variants, warnings=warnings)
    try:
        from PIL import Image, ImageFilter, ImageOps
    except Exception as exc:  # pragma: no cover - depends on optional Pillow install
        return OcsrImageVariantSet(
            variants=variants,
            warnings=[f"Pillow unavailable; OCSR image preprocessing skipped: {exc}"],
        )

    try:
        source = ImageOps.exif_transpose(Image.open(image_path))
    except Exception as exc:
        return OcsrImageVariantSet(variants=variants, warnings=[f"OCSR image preprocessing skipped: {exc}"])

    tmp_dir = Path(tempfile.mkdtemp(prefix="codex_chem_ocsr_"))
    seen: set[str] = set()

    def image_key(image: Any) -> str:
        digest = hashlib.sha256()
        digest.update(str(image.mode).encode())
        digest.update(str(image.size).encode())
        digest.update(image.tobytes())
        return digest.hexdigest()

    def save_variant(name: str, image: Any) -> None:
        key = image_key(image)
        if key in seen:
            return
        seen.add(key)
        path = tmp_dir / f"{name}.png"
        image.save(path)
        variants.append(OcsrImageVariant(name, path, temporary=True))

    try:
        rgb = _flatten_to_white(source)
        cropped = _crop_to_content(rgb)
        scaled = _scale_for_ocsr(cropped)
        padded = ImageOps.expand(scaled, border=_padding_for(scaled), fill="white")
        save_variant("white_padded", padded)

        gray = ImageOps.grayscale(padded)
        high_contrast = ImageOps.autocontrast(gray)
        binary = high_contrast.point(lambda pixel: 0 if pixel < 235 else 255, "L").convert("RGB")
        save_variant("high_contrast_binary", binary)

        thickened = binary.convert("L").filter(ImageFilter.MinFilter(3)).convert("RGB")
        save_variant("high_contrast_thick", thickened)

        if _mean_luminance(gray) < 128:
            save_variant("inverted_dark_background", ImageOps.invert(padded.convert("RGB")))
    except Exception as exc:
        warnings.append(f"OCSR image preprocessing failed; using original image only: {exc}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return OcsrImageVariantSet(variants=variants[:1], warnings=warnings)

    if len(variants) == 1:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return OcsrImageVariantSet(variants=variants, warnings=warnings)
    warnings.append(
        "OCSR image preprocessing generated line-art variants: "
        + ", ".join(variant.name for variant in variants[1:])
    )
    return OcsrImageVariantSet(variants=variants, warnings=warnings, cleanup_dir=tmp_dir)


def _flatten_to_white(image: Any) -> Any:
    from PIL import Image

    rgba = image.convert("RGBA")
    background = Image.new("RGBA", rgba.size, "white")
    background.alpha_composite(rgba)
    return background.convert("RGB")


def _crop_to_content(image: Any) -> Any:
    from PIL import ImageOps

    gray = ImageOps.grayscale(image)
    mask = gray.point(lambda pixel: 255 if pixel < 248 else 0, "L")
    bbox = mask.getbbox()
    if not bbox:
        return image
    left, top, right, bottom = bbox
    width, height = image.size
    margin = max(16, int(max(right - left, bottom - top) * 0.18))
    return image.crop(
        (
            max(0, left - margin),
            max(0, top - margin),
            min(width, right + margin),
            min(height, bottom + margin),
        )
    )


def _scale_for_ocsr(image: Any) -> Any:
    longest = max(image.size)
    if longest <= 0:
        return image
    target = 768
    max_target = 1600
    if longest < target:
        scale = min(4.0, target / longest)
    elif longest > max_target:
        scale = max_target / longest
    else:
        return image
    size = (max(1, int(image.size[0] * scale)), max(1, int(image.size[1] * scale)))
    return image.resize(size)


def _padding_for(image: Any) -> int:
    return max(24, int(max(image.size) * 0.08))


def _mean_luminance(image: Any) -> float:
    histogram = image.histogram()
    total = sum(histogram) or 1
    return sum(value * count for value, count in enumerate(histogram)) / total
