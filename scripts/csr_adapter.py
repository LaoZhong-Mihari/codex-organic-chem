#!/usr/bin/env python3

import argparse
import contextlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


DEFAULT_CSR_PYTHON = (
    Path.home()
    / ".local"
    / "share"
    / "codex-organic-chem"
    / "ocsr-tools"
    / "csr-osx64"
    / "bin"
    / "python"
)


def _bootstrap_csr_python() -> None:
    if os.environ.get("CODEX_CHEM_CSR_ADAPTER_BOOTSTRAPPED") == "1":
        return
    if not DEFAULT_CSR_PYTHON.exists():
        return
    if Path(sys.executable).resolve() == DEFAULT_CSR_PYTHON.resolve():
        return
    os.environ["CODEX_CHEM_CSR_ADAPTER_BOOTSTRAPPED"] = "1"
    os.execv(str(DEFAULT_CSR_PYTHON), [str(DEFAULT_CSR_PYTHON), __file__, *sys.argv[1:]])


def _labels(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if item is not None]
    return [str(value)]


def _candidate(entry: Any) -> Optional[Dict[str, Any]]:
    labels = []  # type: List[str]
    smiles = None  # type: Optional[str]
    if isinstance(entry, tuple) and len(entry) >= 2:
        first, second = entry[0], entry[1]
        if isinstance(first, (list, tuple, set)):
            labels = _labels(first)
            smiles = str(second) if second else None
        else:
            smiles = str(first) if first else None
            labels = _labels(second)
    elif isinstance(entry, dict):
        smiles = entry.get("smiles") or entry.get("canonical_smiles")
        labels = _labels(entry.get("labels") or entry.get("label"))
    elif hasattr(entry, "smiles"):
        smiles = str(getattr(entry, "smiles") or "")
        labels = _labels(getattr(entry, "labels", None) or getattr(entry, "label", None))
    if not smiles:
        return None
    return {
        "tool": "chemschematicresolver",
        "smiles": smiles,
        "confidence": None,
        "labels": labels,
        "adapter_warnings": [
            "ChemSchematicResolver output is legacy OCSR/label resolution; verify labels and R-group attachment manually."
        ],
    }


def main() -> None:
    _bootstrap_csr_python()
    parser = argparse.ArgumentParser(description="Adapter from ChemSchematicResolver to codex-chem JSON lines.")
    parser.add_argument("--input", required=True)
    args = parser.parse_args()

    input_path = str(Path(args.input).expanduser().resolve())
    warnings = []  # type: List[str]
    candidates = []  # type: List[Dict[str, Any]]
    try:
        with contextlib.redirect_stdout(sys.stderr):
            import chemschematicresolver as csr

            results = csr.extract_image(input_path)
        for entry in results or []:
            candidate = _candidate(entry)
            if candidate:
                candidates.append(candidate)
        if not candidates:
            warnings.append("ChemSchematicResolver produced no structure candidates.")
    except Exception as exc:
        warnings.append(f"ChemSchematicResolver adapter failed: {exc}")
    print(json.dumps({"tool": "chemschematicresolver", "candidates": candidates, "adapter_warnings": warnings}))


if __name__ == "__main__":
    main()
