"""Regression tests for Phase 3 chemistry engine accuracy and performance fixes."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from unittest.mock import patch

from codex_organic_chem import kit
from codex_organic_chem.reaction import ReactionParts, mass_charge_balance
from codex_organic_chem.semiempirical import _embed_geometry, xtb_thermo_record
from codex_organic_chem.service import chem_mechanism_render, chem_mechanism_validate
from codex_organic_chem.utils import ToolStatus

ROUTE_FIXTURES = Path(__file__).parent / "fixtures" / "routes"


def _load_route(name: str) -> dict:
    return json.loads((ROUTE_FIXTURES / f"{name}.json").read_text(encoding="utf-8"))


def test_face_accessibility_assigns_re_si_on_achiral_and_chiral_substrates():
    # Achiral prochiral ketone (butan-2-one): C2 is carbonyl carbon, must get Re/Si labels.
    butanone = kit.face("CCC(=O)C", atom=2)
    assert butanone["status"] == "ok"
    assert set(butanone["shielding"]) == {"Re", "Si"}
    assert "Re/Si undefined" not in butanone["face_definition"]

    # Chiral ester where carbonyl substituents tie at first shell and differ via chiral carbon.
    ester = kit.face("CC(=O)OC[C@H](C)Cl", atom=1)
    assert ester["status"] == "ok"
    assert set(ester["shielding"]) == {"Re", "Si"}


def test_symmetry_environments_counts_diastereotopic_13c_signals():
    # L-valine has two diastereotopic methyls -> 5 distinct carbon NMR signals (not 4).
    valine = kit.symmetry_environments("CC(C)[C@H](N)C(=O)O")
    assert valine["n_carbon_environments"] == 5

    # (-)-Menthol has diastereotopic isopropyl methyls -> 10 distinct carbon signals.
    menthol = kit.symmetry_environments("CC1CC[C@H](C(C)C)[C@@H](O)C1")
    assert menthol["n_carbon_environments"] == 10


def test_structure_compare_handles_allenes_isotopes_and_chiral_tautomers():
    # Enantiomeric 2,3-pentadienes
    allenes = kit.compare_structures("CC=[C@]=CC", "CC=[C@@]=CC")
    assert allenes["relation"] == "enantiomers"
    audit = kit.stereo_audit("CC=[C@]=CC")
    assert audit["n_stereocenters"] == 1
    assert audit["chiral"] is True

    # Isotopologues and isotopomers
    assert kit.compare_structures("C[13CH2]O", "CCO")["relation"] == "isotopologues"
    assert kit.compare_structures("[2H]CCO", "CCO[2H]")["relation"] == "isotopomers"

    # Tautomers with inverted spectator stereocentre must not silently pass as plain "tautomers"
    taut = kit.compare_structures("C[C@H](Cl)CC(=O)CC(C)=O", "C[C@@H](Cl)CC(O)=CC(C)=O")
    assert taut["relation"] == "tautomers_different_stereo"


def test_route_check_catches_both_former_benchmark_misses_and_single_letter_id():
    # Miss 1: deleting a Boc methyl in compound 12 of biaryl_amine_amide (N-Boc -> N-CO2iPr)
    biaryl = _load_route("biaryl_amine_amide")
    bad_biaryl = copy.deepcopy(biaryl)
    bad_biaryl["compounds"]["12"]["smiles"] = "C[C@@H](NC(=O)[C@@H]1CCCN1C(=O)OC(C)C)c1ccc(cc1)-c1ccccc1"
    rep1 = kit.route_check(bad_biaryl)
    review_codes_1 = {(i.get("step"), i["code"]) for i in rep1["issues"] if i["severity"] == "review"}
    assert (7, "unexplained_formula_change") in review_codes_1

    # Miss 2: inverting C2 of epoxide 6 in roche_polyketide_fragment ((E)-alkene -> cis-epoxide)
    roche = _load_route("roche_polyketide_fragment")
    bad_roche = copy.deepcopy(roche)
    bad_roche["compounds"]["6"]["smiles"] = "OC[C@@H]1O[C@@H]1[C@@H](C)CO[Si](C)(C)C(C)(C)C"
    rep2 = kit.route_check(bad_roche)
    review_codes_2 = {(i.get("step"), i["code"]) for i in rep2["issues"] if i["severity"] == "review"}
    assert (5, "stereospecific_alkene_addition_mismatch") in review_codes_2

    # Missing single-letter ID "C" must not become inline methane
    typo_route = {
        "compounds": {"A": "CCO", "B": "CC=O"},
        "steps": [{"step": 1, "reactants": ["A"], "product": "C", "conditions": "DMP, CH2Cl2"}],
    }
    rep3 = kit.route_check(typo_route)
    assert any(i["code"] == "unknown_product" for i in rep3["issues"])


def test_route_check_meso_bis_tbdps_and_enone_and_polyol_scaling():
    tbdps = "[Si](c1ccccc1)(c1ccccc1)C(C)(C)C"
    meso_route = {
        "target": f"CC(=O)O[C@@H](CO{tbdps})[C@H](O)CO{tbdps}",
        "compounds": {
            "1": f"O[C@H](CO{tbdps})[C@@H](O)CO{tbdps}",
            "2": f"CC(=O)O[C@@H](CO{tbdps})[C@H](O)CO{tbdps}",
        },
        "steps": [{"step": 1, "reactants": ["1"], "product": "2", "conditions": "Ac2O, lipase, THF"}],
    }
    rep = kit.route_check(meso_route)
    assert rep["steps"][0]["stereo"]["epimerized_without_bond_change"] == []
    assert rep["steps"][0]["stereo"]["inverted_where_only_heteroatom_derivatised"] == []

    # Alpha,beta-unsaturated ester is an ester + alkene, not an enone; MVK is an enone.
    assert "enone" not in kit.functional_groups("CCOC(=O)/C=C/C")["functional_groups"]
    assert kit.functional_groups("C=CC(C)=O")["functional_groups"]["enone"] == 1

    # Tri-acetylation of glycerol (+6C +6H +3O) is recognized under Ac2O conditions.
    triol_route = {
        "target": "CC(=O)OCC(OC(C)=O)COC(C)=O",
        "compounds": {"1": "OCC(O)CO", "2": "CC(=O)OCC(OC(C)=O)COC(C)=O"},
        "steps": [{"step": 1, "reactants": ["1"], "product": "2", "conditions": "Ac2O, pyridine, DMAP"}],
    }
    assert kit.route_check(triol_route)["issues"] == []


def test_mechanism_semantics_and_scene_support_spectator_aromatics_and_pi_attack():
    # 1) Benzyl bromide + NaOH SN2 with spectator phenyl ring and spectator [Na+] counterion
    benzylic_sn2 = {
        "spec_version": "2.1",
        "panels": [
            {
                "title": "Benzylic SN2",
                "reaction_arrow_after": 3,
                "starting_state": {"molecule_indices": [1, 2, 3]},
                "expected_state": {"molecule_indices": [4, 5, 6]},
                "molecules": [
                    {"smiles": "[OH-:1]"},
                    {"smiles": "[cH:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[CH2:8][Br:9]"},
                    {"smiles": "[Na+:10]"},
                    {"smiles": "[cH:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[CH2:8][OH:1]"},
                    {"smiles": "[Br-:9]"},
                    {"smiles": "[Na+:10]"},
                ],
                "lone_pairs": [{"atom_map": 1, "molecule_index": 1, "count": 3, "visible_count": 1}],
                "arrows": [
                    {"from_lone_pair_atom_map": 1, "to_atom_map": 8, "kind": "electron_pair"},
                    {"from_bond": [8, 9], "to_atom_map": 9, "kind": "electron_pair"},
                ],
                "graph_edits": [
                    {"type": "form_bond", "atom_maps": [1, 8], "order": 1, "source_arrow": 1},
                    {"type": "break_bond", "atom_maps": [8, 9], "source_arrow": 2},
                    {"type": "formal_charge", "atom_map": 1, "from": -1, "to": 0},
                    {"type": "formal_charge", "atom_map": 9, "from": 0, "to": -1},
                ],
            }
        ],
    }
    rendered_sn2 = chem_mechanism_render(benzylic_sn2)
    assert rendered_sn2["semantic_checks"]["status"] == "valid"
    assert rendered_sn2["status"] == "ok"

    # 2) Intermolecular alkene pi-attack on H-Br using cross-molecule to_bond: [2, 3]
    alkene_hbr = {
        "spec_version": "2.1",
        "panels": [
            {
                "title": "Alkene protonation",
                "reaction_arrow_after": 2,
                "starting_state": {"molecule_indices": [1, 2]},
                "expected_state": {"molecule_indices": [3, 4]},
                "molecules": [
                    {"smiles": "[CH2:1]=[CH2:2]"},
                    {"smiles": "[H:3][Br:4]"},
                    {"smiles": "[CH2+:1][CH2:2][H:3]"},
                    {"smiles": "[Br-:4]"},
                ],
                "arrows": [
                    {"from_bond": [1, 2], "to_bond": [2, 3], "kind": "electron_pair"},
                    {"from_bond": [3, 4], "to_atom_map": 4, "kind": "electron_pair"},
                ],
                "graph_edits": [
                    {"type": "change_bond", "atom_maps": [1, 2], "from": 2, "to": 1, "source_arrow": 1},
                    {"type": "form_bond", "atom_maps": [2, 3], "order": 1, "source_arrow": 1},
                    {"type": "break_bond", "atom_maps": [3, 4], "source_arrow": 2},
                    {"type": "formal_charge", "atom_map": 1, "from": 0, "to": 1},
                    {"type": "formal_charge", "atom_map": 4, "from": 0, "to": -1},
                ],
            }
        ],
    }
    rendered_pi = chem_mechanism_render(alkene_hbr)
    assert rendered_pi["semantic_checks"]["status"] == "valid"
    assert rendered_pi["status"] == "ok"

    # Also allow to_atom_map: 3 when graph_edits disambiguates the attaching carbon
    alkene_hbr_atom = copy.deepcopy(alkene_hbr)
    alkene_hbr_atom["panels"][0]["arrows"][0] = {
        "from_bond": [1, 2],
        "to_atom_map": 3,
        "kind": "electron_pair",
    }
    assert chem_mechanism_validate(alkene_hbr_atom)["semantic_checks"]["status"] == "valid"

    # 3) Active aromatic SEAr addition (benzene + Br2 -> Wheland sigma-complex + Br-)
    sear = {
        "spec_version": "2.1",
        "panels": [
            {
                "title": "SEAr sigma-complex formation",
                "reaction_arrow_after": 2,
                "starting_state": {"molecule_indices": [1, 2]},
                "expected_state": {"molecule_indices": [3, 4]},
                "molecules": [
                    {"smiles": "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1"},
                    {"smiles": "[Br:7][Br:8]"},
                    {"smiles": "[CH:1]1([Br:7])[CH+:2][CH:3]=[CH:4][CH:5]=[CH:6]1"},
                    {"smiles": "[Br-:8]"},
                ],
                "arrows": [
                    {"from_bond": [1, 2], "to_bond": [1, 7], "kind": "electron_pair"},
                    {"from_bond": [7, 8], "to_atom_map": 8, "kind": "electron_pair"},
                ],
                "graph_edits": [
                    {"type": "change_bond", "atom_maps": [1, 2], "from": 2, "to": 1, "source_arrow": 1},
                    {"type": "form_bond", "atom_maps": [1, 7], "order": 1, "source_arrow": 1},
                    {"type": "break_bond", "atom_maps": [7, 8], "source_arrow": 2},
                    {"type": "formal_charge", "atom_map": 2, "from": 0, "to": 1},
                    {"type": "formal_charge", "atom_map": 8, "from": 0, "to": -1},
                ],
            }
        ],
    }
    sear_res = chem_mechanism_render(sear)
    assert sear_res["semantic_checks"]["status"] == "valid", sear_res["semantic_checks"]
    assert sear_res["status"] == "ok"


def test_semiempirical_boron_uff_fallback_and_failed_opt_blocks_hessian():
    # Phenylboronic acid lacks MMFF94 boron parameters; UFF fallback must succeed without warning
    _mol, xyz, charge, unpaired, warnings = _embed_geometry("OB(O)c1ccccc1")
    assert charge == 0 and unpaired == 0
    assert "B " in xyz
    assert not any("started from the raw embedding" in w for w in warnings)

    # Failed --opt in xtb_thermo_record must return status="error" without running --hess on input.xyz
    calls = []

    def fake_run(cmd, timeout_s=60, cwd=None):
        calls.append(cmd)
        return 1, "", "opt failed"

    with (
        patch(
            "codex_organic_chem.semiempirical.executable_status",
            return_value=ToolStatus(name="xtb", status="available", version="6.7.1"),
        ),
        patch("codex_organic_chem.semiempirical.run_command", side_effect=fake_run),
    ):
        rec = xtb_thermo_record("CCO")
    assert rec.status == "error"
    assert len(calls) == 1  # Only --opt was attempted; --hess was NOT called on unoptimized input.xyz


def test_reaction_balance_does_not_borrow_charge_from_excluded_reagent_atoms():
    bal = mass_charge_balance(ReactionParts(reactants=["[Cu+2]"], reagents=["[Cl-]"], products=["[Cu+]"]))
    assert bal["reagents_counted_in_balance"] is False
    assert bal["status"] == "unbalanced"
    assert bal["charge_delta"] == -1
