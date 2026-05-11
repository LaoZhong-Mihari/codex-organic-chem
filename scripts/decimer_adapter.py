#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import sys
from statistics import mean


def _confidence_value(tokens: object) -> float | None:
    if not isinstance(tokens, list):
        return None
    values: list[float] = []
    for item in tokens:
        if not isinstance(item, tuple) or len(item) != 2:
            continue
        try:
            values.append(float(item[1]))
        except (TypeError, ValueError):
            continue
    if not values:
        return None
    return float(mean(values))


def main() -> None:
    parser = argparse.ArgumentParser(description="Adapter from DECIMER Image Transformer to codex-chem JSON lines.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--hand-drawn", action="store_true")
    args = parser.parse_args()

    with contextlib.redirect_stdout(sys.stderr):
        from DECIMER import predict_SMILES

        smiles, token_confidence = predict_SMILES(args.input, confidence=True, hand_drawn=args.hand_drawn)
    payload = {
        "tool": "decimer_handdrawn" if args.hand_drawn else "decimer",
        "smiles": smiles,
        "confidence": _confidence_value(token_confidence),
        "adapter_warnings": [],
    }
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
