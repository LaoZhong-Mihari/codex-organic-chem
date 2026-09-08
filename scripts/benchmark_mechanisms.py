"""Reproducible graph/render/round-trip checks; not a statistical model benchmark.

Run: uv run python scripts/benchmark_mechanisms.py --output-dir outputs/mechanisms
"""

from __future__ import annotations

import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

from rdkit import RDLogger

from codex_organic_chem.service import (
    chem_mechanism_render,
    chem_mechanism_validate,
    chem_normalize_structure,
    chem_reaction_analyze,
)

FIXTURES = Path(__file__).resolve().parents[1] / "tests/fixtures/mechanisms"


def negative_specs(sn2: dict) -> dict[str, dict]:
    cases = {}

    def case(name):
        s = deepcopy(sn2)
        cases[name] = s
        return s["panels"][0]

    case("missing_leaving_group_arrow")["arrows"].pop()
    case("attack_wrong_atom")["arrows"][0]["to_atom_map"] = 3
    case("hydrogen_charge_imbalance")["molecules"][2]["smiles"] = "[OH2+:1][CH3:2]"
    case("nonexistent_source_bond")["arrows"][1]["from_bond"] = [1, 3]
    case("arrowhead_electron_mismatch")["arrows"][0]["electron_count"] = 1
    case("wrong_graph_edit")["graph_edits"][0]["order"] = 2
    case("duplicate_atom_map")["molecules"][0]["smiles"] = "[OH-:2]"
    case("wrong_state_side")["arrows"][0]["from_molecule_index"] = 3
    case("missing_product")["expected_state"]["molecule_indices"] = [3]
    case("fake_charge_label")["charges"][0]["label"] = "+"
    case("wrong_element_identity")["molecules"][3]["smiles"] = "[Cl-:3]"
    case("wrong_isotope_identity")["molecules"][3]["smiles"] = "[81Br-:3]"
    case("electron_source_overdraw")["arrows"].extend(
        [deepcopy(sn2["panels"][0]["arrows"][1])] * 2
    )
    return cases


def benchmark(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "method": "GPT-6-authored curated fixtures plus deterministic negative controls; no independent model sampling or old-model A/B",
        "positive": {},
        "negative": {},
    }
    for path in sorted(FIXTURES.glob("*.json")):
        spec = json.loads(path.read_text())
        normals = [
            chem_normalize_structure(smiles=m["smiles"])
            for p in spec["panels"]
            for m in p["molecules"]
        ]
        result = chem_mechanism_render(spec, str(output_dir / path.stem))
        bundle = json.loads(
            ET.fromstring(result["svg"])
            .find("{http://www.w3.org/2000/svg}metadata")
            .text
        )
        digest = hashlib.sha256(
            json.dumps(
                bundle["spec"],
                sort_keys=True,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        balances = []
        for panel in spec["panels"]:
            molecules = panel["molecules"]
            before = ".".join(
                molecules[i - 1]["smiles"]
                for i in panel["starting_state"]["molecule_indices"]
            )
            after = ".".join(
                molecules[i - 1]["smiles"]
                for i in panel["expected_state"]["molecule_indices"]
            )
            balances.append(
                chem_reaction_analyze(f"{before}>>{after}", mode="sanity_check")[
                    "analysis"
                ]["mass_charge_balance"]
            )
        roundtrip = (
            bundle["spec"] == spec
            and bundle["trace"] == result["mechanism_trace"]
            and digest == bundle["trace"]["source_spec_sha256"]
            and chem_mechanism_validate(bundle["spec"])["semantic_checks"]["status"]
            == "valid"
        )
        report["positive"][path.stem] = {
            "render_status": result["status"],
            "semantic_status": result["semantic_checks"]["status"],
            "blocking_issues": result["publication_checks"]["blocking_issues"],
            "warnings": result["warnings"],
            "roundtrip_exact": roundtrip,
            "step_count": len(result["semantic_checks"]["steps"]),
            "balances": balances,
            "normalization": normals,
            "figure_audit": result["publication_checks"]["figure_audit"],
            "source_spec_sha256": digest,
        }
    sn2 = json.loads((FIXTURES / "sn2.json").read_text())
    for name, spec in negative_specs(sn2).items():
        result = chem_mechanism_validate(spec)
        report["negative"][name] = {
            "blocked": bool(result["blocking_issues"]),
            "issues": result["blocking_issues"],
        }
    report["summary"] = {
        "positive_total": len(report["positive"]),
        "semantic_valid": sum(
            r["semantic_status"] == "valid" for r in report["positive"].values()
        ),
        "figure_unblocked": sum(
            not r["blocking_issues"] for r in report["positive"].values()
        ),
        "roundtrip_exact": sum(
            r["roundtrip_exact"] for r in report["positive"].values()
        ),
        "negative_total": len(report["negative"]),
        "negative_blocked": sum(r["blocked"] for r in report["negative"].values()),
    }
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.warning")
    result = benchmark(args.output_dir)
    print(json.dumps(result["summary"], indent=2))
