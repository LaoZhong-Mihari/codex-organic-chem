import copy
import json
import re
from pathlib import Path

import pytest
from rdkit import Chem

from codex_organic_chem import kit
from codex_organic_chem import route_check as rc
from codex_organic_chem.cli import main as cli_main

FIXTURES = Path(__file__).parent / "fixtures" / "routes"
BOC = "C(=O)OC(C)(C)C"
TBS = "[Si](C)(C)C(C)(C)C"


def load(name: str) -> dict:
    return json.loads((FIXTURES / f"{name}.json").read_text(encoding="utf-8"))


def codes(report: dict, severity: str | None = None) -> set[tuple]:
    return {
        (issue.get("step"), issue["code"])
        for issue in report["issues"]
        if severity is None or issue["severity"] == severity
    }


def with_compound(route: dict, key: str, smiles: str) -> dict:
    mutated = copy.deepcopy(route)
    mutated["compounds"][key]["smiles"] = smiles
    return mutated


def test_tables_are_internally_consistent():
    for table in (rc.FUNCTIONAL_GROUPS, rc.PROTECTING_GROUPS):
        for name, smarts in table.items():
            assert Chem.MolFromSmarts(smarts) is not None, name
    known = set(rc.FUNCTIONAL_GROUPS) | set(rc.PROTECTING_GROUPS) | {
        "stereocenter_between_two_acceptors",
        "epimerizable_alpha_stereocenter",
        "allylic_stereocenter",
        "beta_hydroxy_carbonyl",
        "beta_leaving_group_carbonyl",
    }
    for name, pattern, risks in rc.REAGENT_RULES:
        re.compile(pattern)
        assert set(risks) <= known, name
    assert rc._SELECTIVITY_GROUPS <= known
    for _, _, words in rc.NAMED_DELTAS:
        re.compile(words)


@pytest.mark.parametrize(
    "name", ["roche_polyketide_fragment", "hydroxyproline_double_inversion", "biaryl_amine_amide"]
)
def test_correct_routes_have_no_blocking_issue_and_reach_the_target(name):
    report = kit.route_check(load(name))
    assert report["n_blocking"] == 0
    assert report["final_vs_target"]["relation"] == "identical"
    assert report["n_steps"] == 8 and report["longest_linear_sequence"] == 8
    for step in report["steps"]:
        assert step["stereo"]["epimerized_without_bond_change"] == []
        assert step["stereo"]["ambiguous"] == []


def test_clean_polyketide_route_is_fully_consistent():
    report = kit.route_check(load("roche_polyketide_fragment"))
    assert report["status"] == "consistent" and report["issues"] == []
    assert report["stereo_summary"]["centres_created_in_steps"] == 2
    assert report["protecting_groups_in_final_product"] == {"acetonide": 1}
    deltas = [step["atom_balance"]["delta_product_minus_reactants"] for step in report["steps"]]
    assert deltas[0] == "+6C +14H +1Si" and deltas[-1] == "-6C -14H -1Si"
    assert report["steps"][2]["atom_balance"]["status"] == "closed"


def test_inversions_are_followed_by_spatial_arrangement():
    report = kit.route_check(load("hydroxyproline_double_inversion"))
    inverted_steps = {step["step"] for step in report["steps"] if step["stereo"]["inverted_with_substituent_replaced"]}
    assert inverted_steps == {2, 5}  # Mitsunobu and azide displacement
    assert report["stereo_summary"]["centres_inverted_in_steps"] == 2


def test_wrong_wedges_are_caught_at_the_step_where_they_appear():
    route = load("hydroxyproline_double_inversion")
    mitsunobu_retained = with_compound(
        route, "4", f"COC(=O)[C@@H]1C[C@@H](OC(=O)c2ccc(cc2)[N+](=O)[O-])CN1{BOC}"
    )
    assert (2, "retention_under_inverting_conditions") in codes(kit.route_check(mitsunobu_retained), "review")
    mesylate_inverted = with_compound(route, "6", f"COC(=O)[C@@H]1C[C@@H](OS(C)(=O)=O)CN1{BOC}")
    assert (4, "inversion_where_only_the_heteroatom_was_derivatised") in codes(
        kit.route_check(mesylate_inverted), "review"
    )
    azide_retained = with_compound(route, "7", f"COC(=O)[C@@H]1C[C@H](N=[N+]=[N-])CN1{BOC}")
    assert (5, "retention_under_inverting_conditions") in codes(kit.route_check(azide_retained), "review")
    # an untouched centre flipped in a step that does nothing to it
    alpha_flipped = with_compound(route, "9", f"COC(=O)[C@H]1C[C@@H](NC(C)=O)CN1{BOC}")
    assert (7, "configuration_changed_without_bond_change") in codes(kit.route_check(alpha_flipped), "review")


def test_wrong_final_product_is_blocking_and_located():
    route = load("hydroxyproline_double_inversion")
    wrong = with_compound(route, "10", f"OC(=O)[C@@H]1C[C@H](NC(C)=O)CN1{BOC}")
    report = kit.route_check(wrong)
    assert report["status"] == "blocked"
    assert (None, "final_product_is_not_the_target") in codes(report, "blocking")
    assert report["final_vs_target"]["relation"] == "diastereomers"
    assert report["final_vs_target"]["n_differing"] == 1
    forgot_deprotection = copy.deepcopy(load("roche_polyketide_fragment"))
    forgot_deprotection["steps"].pop()
    report = kit.route_check(forgot_deprotection)
    assert report["final_vs_target"]["relation"] == "different_formula"
    assert report["final_vs_target"]["formula_delta_b_minus_a"] == "-6C -14H -1Si"


def test_miscounted_atoms_are_flagged_against_the_conditions():
    route = load("roche_polyketide_fragment")
    # one methyl missing from the TBS group: no reagent explains +C5H12Si
    short_silyl = with_compound(route, "2", "COC(=O)[C@@H](C)CO[Si](C)(C)C(C)C")
    assert (1, "unexplained_formula_change") in codes(kit.route_check(short_silyl), "review")
    # An extra CH2 in the alcohol: the formula change then equals that of a
    # methyl-ester reduction, so only the atom mapping can see the error.
    homologated = with_compound(route, "5", f"OCC/C=C/[C@@H](C)CO{TBS}")
    report = kit.route_check(homologated)
    assert report["steps"][3]["atom_balance"]["status"] == "consistent_with_conditions"
    assert (4, "atoms_relocated") in codes(report, "review")


def test_declared_balance_must_close_exactly():
    route = load("roche_polyketide_fragment")
    broken = with_compound(route, "4", f"CCOC(=O)/C=C/[C@@H](C)CCO{TBS}")
    report = kit.route_check(broken)
    assert (3, "atom_balance_does_not_close") in codes(report, "blocking")
    assert report["steps"][2]["atom_balance"]["residual"] == "+1C +2H"


def test_regiochemistry_drift_and_undefined_centres_are_reported():
    route = load("roche_polyketide_fragment")
    shifted = with_compound(route, "7", f"OC[C@H](O)C[C@@H](C)CO{TBS}")  # 1,2- instead of 1,3-diol
    assert (7, "atoms_relocated") in codes(kit.route_check(shifted), "review")
    undefined = with_compound(route, "7", f"OCCC(O)[C@@H](C)CO{TBS}")
    assert (6, "undefined_stereo_in_product") in codes(kit.route_check(undefined), "review")
    racemic_grignard = kit.route_check(load("biaryl_amine_amide"))
    assert (2, "undefined_stereo_in_product") in codes(racemic_grignard, "review")


def test_bystander_groups_are_checked_against_the_reagent():
    route = {
        "target": "OC[C@@H](C)CO[Si](C)(C)C(C)(C)C",
        "compounds": {
            "1": "COC(=O)[C@@H](C)CO[Si](C)(C)C(C)(C)C",
            "2": "OC[C@@H](C)CO[Si](C)(C)C(C)(C)C",
        },
        "steps": [
            {"step": 1, "reactants": ["1"], "product": "2", "conditions": "LiAlH4, THF, 0 °C"},
        ],
    }
    assert codes(kit.route_check(route), "review") == set()
    # fluoride would take the silyl ether off; the product is drawn with it intact
    route["steps"][0]["conditions"] = "TBAF, THF, rt"
    report = kit.route_check(route)
    step_codes = {code for _, code in codes(report, "review")}
    assert "bystander_group_vs_reagent" in step_codes
    assert {"group": "silyl_ether", "count": 1} == {
        k: v for k, v in report["steps"][0]["compatibility_checks"][0].items() if k in {"group", "count"}
    }
    two_ketones = {
        "compounds": {"1": "CC(=O)CCC(=O)C1CC1", "2": "CC(O)CCC(=O)C1CC1"},
        "steps": [{"step": 1, "reactants": ["1"], "product": "2", "conditions": "NaBH4, MeOH, 0 °C"}],
    }
    report = kit.route_check(two_ketones)
    found = {code for _, code in codes(report, "review")}
    assert {"selectivity_among_equal_groups", "undefined_stereo_in_product"} <= found
    assert "bystander_group_vs_reagent" not in found  # one question per group, not two
    assert report["steps"][0]["functional_groups"]["bystanders_that_must_survive"] == {"ketone": 1}


def test_delta_interpretation_needs_support_from_the_conditions():
    assert rc.interpret_delta({"C": 6, "H": 14, "Si": 1}, "TBSCl, imidazole")["status"] == "consistent_with_conditions"
    assert rc.interpret_delta({"C": 5, "H": 12, "Si": 1}, "TBSCl, imidazole")["status"] == "unexplained"
    assert (
        rc.interpret_delta({"C": 1, "H": 2}, "NaBH4, MeOH")["status"]
        == "named_but_conditions_do_not_mention_a_matching_reagent"
    )
    combined = rc.interpret_delta({"C": 7, "H": 4}, "BnBr, NaH then DMP")
    assert combined["status"] == "explained_only_as_a_combination"
    assert rc.interpret_delta({"C": -5, "H": -8, "O": -2}, "TFA, CH2Cl2")["status"] == "consistent_with_conditions"


def test_functional_group_and_liability_inventory():
    groups = kit.functional_groups(f"CC(=O)OC[C@H](O{TBS})CN=[N+]=[N-]")
    assert groups["functional_groups"] == {"azide": 1, "ester": 1}
    assert groups["protecting_groups"] == {"TBS_ether": 1, "acetate_ester": 1, "silyl_ether": 1}
    assert kit.functional_groups("C[C@H](C(C)=O)C(=O)OC")["liabilities"] == {"stereocenter_between_two_acceptors": 1}
    aldol = kit.functional_groups("C=CC[C@H](C)C(=O)CC(O)C")["liabilities"]
    assert aldol == {"beta_hydroxy_carbonyl": 1, "epimerizable_alpha_stereocenter": 1}


def test_malformed_routes_come_back_as_data():
    assert kit.route_check({"steps": []})["status"] == "error"
    report = kit.route_check(
        {
            "compounds": {"1": "CCO", "2": "C1CC"},
            "steps": [{"step": 1, "reactants": ["1"], "product": "2", "conditions": "DMP"}],
        }
    )
    assert report["status"] == "blocked"
    assert {code for _, code in codes(report, "blocking")} >= {"unparseable_structure", "unknown_product"}


def test_route_becomes_a_drawable_scheme(tmp_path):
    route = load("roche_polyketide_fragment")
    spec = kit.route_to_figure_spec(route)
    arrows = [item for row in spec["rows"] for item in row["items"] if item.get("type") == "arrow"]
    assert len(arrows) == 8
    # Reagents above the arrow, solvent and temperature below; the comma in a locant survives.
    assert (arrows[6]["label"], arrows[6]["sublabel"]) == ("2,2-dimethoxypropane, PPTS", "CH2Cl2, rt")
    first = next(item for row in spec["rows"] for item in row["items"] if item.get("smiles"))
    assert (first["label"], first["code"]) == ("1", "Roche ester")  # number in bold, name underneath
    assert spec["max_width"] and len(spec["rows"]) == 1  # the renderer wraps the row at arrows
    from codex_organic_chem.service import chem_route_check

    report = chem_route_check(route, figure_output_dir=str(tmp_path))
    assert report["figure"]["status"] in {"ok", "ok_with_warnings"}
    assert (tmp_path / "route.svg").exists()


def test_long_compound_names_no_longer_block_the_route_figure(tmp_path):
    # Compound 1's name is wider than its structure; it used to run off the canvas.
    from codex_organic_chem.service import chem_route_check

    report = chem_route_check(load("hydroxyproline_double_inversion"), figure_output_dir=str(tmp_path))
    assert report["figure"]["status"] in {"ok", "ok_with_warnings"}
    layout_warnings = [w for w in report["figure"]["warnings"] if "PNG rasterization" not in w]
    assert layout_warnings == []


def test_conditions_split_into_reagents_and_settings():
    assert rc.split_conditions("Boc2O, Et3N, CH2Cl2, rt") == ("Boc2O, Et3N", "CH2Cl2, rt")
    assert rc.split_conditions("Pd(PPh3)4, K2CO3, dioxane/H2O, 90 °C") == ("Pd(PPh3)4, K2CO3", "dioxane/H2O, 90 °C")
    assert rc.split_conditions("H2 (1 atm), Pd/C, MeOH, rt") == ("H2 (1 atm), Pd/C", "MeOH, rt")
    assert rc.split_conditions("LiOH, THF/H2O 3:1, 0 °C, 2 h") == ("LiOH", "THF/H2O 3:1, 0 °C, 2 h")
    assert rc.split_conditions("toluene, reflux") == ("", "toluene, reflux")
    assert rc.split_conditions("") == ("", "")


def test_cli_exposes_the_audit(tmp_path, capsys):
    cli_main(["route-check", "--example"])
    example = json.loads(capsys.readouterr().out)
    path = tmp_path / "route.json"
    path.write_text(json.dumps(example), encoding="utf-8")
    cli_main(["route-check", "--route", str(path), "--issues-only"])
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "consistent" and "steps" not in report
    cli_main(["compare", "--a", "C[C@H](O)CC", "--b", "C[C@@H](O)CC"])
    assert json.loads(capsys.readouterr().out)["relation"] == "enantiomers"
    cli_main(["enumerate", "--smiles", "CC1CCCCC1", "--task", "substitution"])
    assert json.loads(capsys.readouterr().out)["n_including_stereoisomers"] == 12
