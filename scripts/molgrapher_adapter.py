#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description="Adapter from MolGrapher to codex-chem JSON lines.")
    parser.add_argument("--input", required=True)
    args = parser.parse_args()

    os.environ.setdefault("DYLD_FALLBACK_LIBRARY_PATH", "/opt/homebrew/lib")

    with contextlib.redirect_stdout(sys.stderr):
        from molgrapher.models.molgrapher_model import MolgrapherModel

        model = MolgrapherModel(
            {
                "force_cpu": True,
                "force_no_multiprocessing": True,
                "num_threads_pytorch": 4,
                "num_processes_mp": 1,
                "chunk_size": 1,
                "assign_stereo": False,
                "align_rdkit_output": False,
                "remove_captions": False,
                "save_mol_folder": "",
                "preprocess": True,
                "clean": False,
                "visualize": False,
                "visualize_rdkit": False,
            }
        )
        annotations = model.predict_batch([args.input])
    candidates = []
    warnings: list[str] = []
    for annotation in annotations:
        smiles = annotation.get("smi")
        if not smiles:
            continue
        candidates.append(
            {
                "tool": "molgrapher",
                "smiles": smiles,
                "confidence": annotation.get("conf"),
                "labels": {
                    "abbreviations": annotation.get("abbreviations", []),
                    "abbreviations_ocr": annotation.get("abbreviations_ocr", []),
                    "file_info": annotation.get("file-info", {}),
                },
            }
        )
    if not candidates:
        warnings.append("MolGrapher produced no SMILES annotations.")
    print(json.dumps({"tool": "molgrapher", "candidates": candidates, "adapter_warnings": warnings}))


if __name__ == "__main__":
    main()
