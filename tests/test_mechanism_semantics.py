import hashlib
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import pytest

from codex_organic_chem.mechanism_semantics import parse_molecule
from codex_organic_chem.service import (
    chem_mechanism_draft,
    chem_mechanism_render,
    chem_mechanism_validate,
)

FIXTURES = Path(__file__).parent / "fixtures/mechanisms"


def fixture(name="sn2"):
    return json.loads((FIXTURES / f"{name}.json").read_text())


@pytest.mark.parametrize(
    "name", ["sn2", "e2", "proton_transfer", "carbonyl_addition", "homolysis"]
)
def test_curated_elementary_steps_replay_and_round_trip(name):
    spec = fixture(name)
    result = chem_mechanism_render(spec)
    assert result["semantic_checks"]["status"] == "valid"
    assert not result["publication_checks"]["blocking_issues"]
    root = ET.fromstring(result["svg"])
    bundle = json.loads(root.find("{http://www.w3.org/2000/svg}metadata").text)
    assert bundle["spec"] == spec
    assert bundle["trace"] == result["mechanism_trace"]
    encoded = json.dumps(
        spec, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode()
    assert hashlib.sha256(encoded).hexdigest() == bundle["trace"]["source_spec_sha256"]
    assert (
        chem_mechanism_validate(bundle["spec"])["semantic_checks"]["status"] == "valid"
    )
    for step in bundle["trace"]["steps"]:
        start, end = (
            step["starting_state"]["resolved_state"],
            step["expected_state"]["resolved_state"],
        )
        assert start["atom_maps"] == end["atom_maps"]
        assert start["formal_charge"] == end["formal_charge"]
        # A map is unique WITHIN a state, even when repeated on the other side.
        assert len(start["atom_maps"]) == sum(
            len(m["atom_maps"]) for m in start["molecules"]
        )
        check = step["semantic_validation"]
        assert check["arrow_predicted_changes"] == check["observed_changes"]


@pytest.mark.parametrize(
    "mutate, issue",
    [
        (lambda p: p["arrows"].pop(), "do not generate"),
        (lambda p: p["arrows"][0].update(to_atom_map=3), "do not generate"),
        (lambda p: p["arrows"][1].update(from_bond=[1, 3]), "source bond"),
        (lambda p: p["arrows"][0].update(electron_count=1), "electron_count"),
        (lambda p: p["graph_edits"][0].update(order=2), "actual bond changes"),
        (lambda p: p["graph_edits"][0].update(source_arrow=2), "source_arrow"),
        (lambda p: p["molecules"][2].update(smiles="[OH2+:1][CH3:2]"), "not conserved"),
        (lambda p: p["molecules"][3].update(smiles="[Cl-:3]"), "identity changed"),
        (lambda p: p["molecules"][0].update(smiles="[OH-:2]"), "duplicate atom map"),
        (
            lambda p: p["arrows"][0].update(from_molecule_index=3),
            "starting-state anchor",
        ),
        (
            lambda p: p["expected_state"].update(molecule_indices=[3]),
            "reaction_arrow_after split",
        ),
        (
            lambda p: p["starting_state"].update(molecule_indices=[1, 2, 3, 4]),
            "disjoint",
        ),
        (lambda p: p["charges"][0].update(label="+"), "displayed charge"),
        (
            lambda p: p["arrows"].extend([deepcopy(p["arrows"][1])] * 2),
            "more electrons",
        ),
        (
            lambda p: p["expected_state"].update(
                reaction_smiles="[OH-:1].[CH3:2][Br:3]"
            ),
            "declared reaction_smiles",
        ),
    ],
)
def test_chemically_inconsistent_figures_are_blocked(mutate, issue):
    spec = fixture()
    mutate(spec["panels"][0])
    result = chem_mechanism_validate(spec)
    assert any(issue in error for error in result["blocking_issues"]), result
    assert not result["ready_for_human_review"]


def test_atom_mapping_preserves_transferred_hydrogen_in_smiles_and_molfile():
    mol = parse_molecule({"smiles": "[O:1]([H:2])[H:3]"})
    assert [a.GetAtomMapNum() for a in mol.GetAtoms()] == [1, 2, 3]
    from rdkit import Chem

    roundtrip = parse_molecule({"molfile": Chem.MolToMolBlock(mol)})
    assert [a.GetAtomMapNum() for a in roundtrip.GetAtoms()] == [1, 2, 3]


def test_trace_uses_actual_diff_instead_of_manual_summary():
    spec = fixture()
    spec["panels"][0]["bond_changes"] = {"formed": [{"atom_maps": [2, 99]}]}
    spec["panels"][0]["charge_changes"] = [{"atom_map": 3, "from": 0, "to": 1}]
    trace = chem_mechanism_render(spec)["mechanism_trace"]["steps"][0]
    assert trace["bond_changes"]["formed"] == [{"atom_maps": [1, 2], "order": "SINGLE"}]
    assert trace["charge_changes"][-1] == {"atom_map": 3, "from": 0, "to": -1}


def test_last_panel_does_not_masquerade_as_its_own_product():
    spec = fixture()
    spec["spec_version"] = "2.0"
    p = spec["panels"][0]
    p.pop("starting_state")
    p.pop("expected_state")
    p.pop("reaction_arrow_after")
    p["molecules"] = p["molecules"][:2]
    for m in p["molecules"]:
        m.pop("side", None)
    p["charges"] = p["charges"][:1]
    r = chem_mechanism_render(spec)
    assert r["semantic_checks"]["status"] == "incomplete"
    assert r["mechanism_trace"]["steps"][0]["expected_state"] == {
        "status": "not_supplied"
    }


def test_prose_does_not_excuse_missing_charge_or_species():
    spec = fixture()
    p = spec["panels"][0]
    p["molecules"][3]["smiles"] = "[Br:3]"
    p["charge_explanation"] = "The solvent handles charge."
    assert chem_mechanism_validate(spec)["semantic_checks"]["status"] == "invalid"


def test_sequential_steps_cannot_teleport_to_different_intermediate():
    spec = fixture("carbonyl_addition")
    spec["panels"][1]["molecules"][0]["smiles"] = "[C:1](#[N:2])[CH:3]=[O:4]"
    result = chem_mechanism_validate(spec)
    assert any("previous expected state" in e for e in result["blocking_issues"])


def test_half_homolysis_is_not_a_complete_step():
    spec = fixture("homolysis")
    spec["panels"][0]["arrows"].pop()
    r = chem_mechanism_validate(spec)
    assert any("paired fishhooks" in e for e in r["blocking_issues"])


def test_radical_dots_and_export_populations_remain_visible():
    r = chem_mechanism_render(fixture("homolysis"))
    root = ET.fromstring(r["svg"])
    dots = [
        c
        for c in root.iter("{http://www.w3.org/2000/svg}path")
        if "stroke-width:0.0px" in c.get("style", "")
    ]
    assert len(dots) == 2
    for dot in dots:
        assert dot.get("d") and "Z" in dot.get("d")
    atoms = r["mechanism_trace"]["steps"][0]["expected_state"]["resolved_state"][
        "molecules"
    ]
    assert sum(a["radical_electrons"] for m in atoms for a in m["atoms"]) == 2
    assert sum(a.get("r", 0) for m in r["chemdoodle_json"]["m"] for a in m["a"]) == 2


@pytest.mark.parametrize(
    "smiles, issue",
    [
        ("[CH3:1][C@H:2]([OH:3])[CH2:4][CH3:5]", "stereo"),
        ("[CH3:1]/[CH:2]=[CH:3]/[CH3:4]", "stereo"),
        ("[13CH4:1]", "isotope"),
    ],
)
def test_canvas_does_not_silently_erase_unsupported_chemical_labels(smiles, issue):
    s = {"panels": [{"molecules": [{"smiles": smiles}]}]}
    assert any(issue in e for e in chem_mechanism_validate(s)["blocking_issues"])


def test_aromatic_ring_drawing_keeps_alternating_bonds():
    s = {
        "presentation_mode": "publication",
        "panels": [
            {"molecules": [{"smiles": "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1"}]}
        ],
    }
    r = chem_mechanism_render(s)
    orders = [b["o"] for b in r["chemdoodle_json"]["m"][0]["b"]]
    assert sorted(orders) == [1, 1, 1, 2, 2, 2]


def test_rule_draft_is_only_a_reaction_preview():
    r = chem_mechanism_draft("CBr.[OH-]>>CO.[Br-]")
    assert "M 130" not in r["steps"][0]["rendered_svg"]


def test_display_metadata_cannot_change_an_atom_or_create_electrons():
    spec = fixture()
    spec["panels"][0]["molecules"][0]["atom_labels"] = {"1": "Cl"}
    assert any(
        "atom label override" in e
        for e in chem_mechanism_validate(spec)["blocking_issues"]
    )
    spec = fixture()
    spec["panels"][0]["lone_pairs"][0]["count"] = 4
    assert any(
        "lone-pair annotation" in e
        for e in chem_mechanism_validate(spec)["blocking_issues"]
    )


def test_strict_spec_does_not_skip_an_unchecked_panel():
    spec = fixture()
    spec["panels"].append({"molecules": [{"smiles": "[CH4:99]"}]})
    assert chem_mechanism_validate(spec)["semantic_checks"]["status"] == "invalid"


def test_validate_cli_reports_semantics():
    process = subprocess.run(
        [
            sys.executable,
            "-m",
            "codex_organic_chem.cli",
            "mechanism-validate",
            "--spec",
            str(FIXTURES / "e2.json"),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(process.stdout)["semantic_checks"]["status"] == "valid"
