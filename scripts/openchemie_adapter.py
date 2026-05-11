#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
from pathlib import Path


def _molecule_payload(input_path: str) -> dict:
    with open(os.devnull, "w", encoding="utf-8") as sink, contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
        import numpy as np
        from PIL import Image
        from openchemie import OpenChemIE

        image = np.array(Image.open(input_path).convert("RGB"))
        model = OpenChemIE(device="cpu")
        result = model.molscribe.predict_images([image], return_confidence=True, batch_size=1)[0]
    return {
        "tool": "openchemie_molscribe",
        "smiles": result.get("smiles"),
        "molfile": result.get("molfile"),
        "confidence": result.get("confidence"),
        "adapter_warnings": [],
    }


def _reaction_payload(input_path: str) -> dict:
    with open(os.devnull, "w", encoding="utf-8") as sink, contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
        from PIL import Image
        from openchemie import OpenChemIE

        image = Image.open(input_path).convert("RGB")
        model = OpenChemIE(device="cpu")
        result = model.extract_reactions_from_figures([image], batch_size=1, molscribe=True, ocr=True)[0]
    candidates = []
    for reaction in result.get("reactions", []):
        candidates.append(
            {
                "tool": "openchemie_rxnscribe",
                "reaction": reaction,
                "adapter_warnings": [
                    "OpenChemIE reaction output contains structured roles/bboxes; convert to reaction SMILES after manual review."
                ],
            }
        )
    return {
        "tool": "openchemie_rxnscribe",
        "candidates": candidates,
        "adapter_warnings": [] if candidates else ["OpenChemIE produced no reaction candidates."],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Adapter from OpenChemIE to codex-chem JSON lines.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--mode", choices=("molecule", "reaction"), default="molecule")
    args = parser.parse_args()

    input_path = str(Path(args.input).expanduser().resolve())
    try:
        payload = _reaction_payload(input_path) if args.mode == "reaction" else _molecule_payload(input_path)
    except Exception as exc:
        payload = {
            "tool": "openchemie",
            "adapter_warnings": [f"OpenChemIE adapter failed: {exc}"],
        }
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
