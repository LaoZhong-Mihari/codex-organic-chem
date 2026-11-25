"""Regression tests for Phase 4 figure rendering, OCSR, CDXML, and workbench improvements."""

from __future__ import annotations

import importlib.util
import time
import xml.etree.ElementTree as ET
from pathlib import Path

from codex_organic_chem import kit
from codex_organic_chem import structure_review as sr
from codex_organic_chem.route_figure import _text_width, _wrap_caption, render_route_figure, subscript_formulas
from codex_organic_chem.scheme_ocsr import compare_smiles, resolve_placeholders
from codex_organic_chem.service import chem_mechanism_render


def _load_script_module(name: str):
    script_path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, script_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_resolve_placeholders_preserves_stereochemistry_and_mem_connectivity():
    # Replacing [*:1] at slot 1 of [C@H] must NOT invert the stereocentre!
    # C[C@H](O)C(=O)OC is methyl (S)-lactate; C[C@H]([*:1])C(=O)OC with 1=OH must give the SAME configuration.
    res = resolve_placeholders("C[C@H]([*:1])C(=O)OC", {"1": "OH"})
    assert res["status"] == "resolved"
    cmp = kit.compare(res["smiles"], "C[C@H](O)C(=O)OC")
    assert cmp["relation"] == "identical", (res["smiles"], cmp)

    # Replacing [*:1] with H on an explicit quaternary chiral centre preserves handedness
    res_h = resolve_placeholders("CC[C@]([*:1])(C)O", {"1": "H"})
    assert res_h["status"] == "resolved"
    assert kit.compare(res_h["smiles"], "CC[C@H](C)O")["relation"] == "identical"

    # R-O-MEM (RO-[*:1] with 1=MEM) must be an acetal R-O-CH2-O-CH2CH2-OCH3, NOT a peroxide R-O-O-...
    mem_ether = resolve_placeholders("CCO[*:1]", {"1": "MEM"})
    assert mem_ether["status"] == "resolved"
    assert "OO" not in mem_ether["smiles"]
    assert kit.compare(mem_ether["smiles"], "CCOCOCH2CH2OC".replace("CH2", "C"))["relation"] == "identical"

    # Expanded organic abbreviations (OTBS, Boc, OTf, Bpin)
    assert resolve_placeholders("OCC[*:1]", {"1": "OTBS"})["status"] == "resolved"
    assert resolve_placeholders("NCC[*:1]", {"1": "NHBoc"})["status"] == "resolved"
    assert resolve_placeholders("c1ccccc1[*:1]", {"1": "Bpin"})["status"] == "resolved"

    # Stereo-aware OCSR benchmark comparison
    stereo_cmp = compare_smiles("C[C@H](O)CC", "C[C@@H](O)CC")
    assert stereo_cmp["exact_graph_match"] is True
    assert stereo_cmp["exact_stereo_match"] is False


def test_molscribe_adapter_chiral_direction_and_hydride_labels():
    mod = _load_script_module("molscribe_adapter")
    s1, _, w1 = mod._molscribe_graph_to_smiles(
        {
            "atoms": [
                {"atom_symbol": "C"},
                {"atom_symbol": "[C@H]"},
                {"atom_symbol": "OH"},
                {"atom_symbol": "[CF3]"},
            ],
            "bonds": [
                {"bond_type": "single", "endpoint_atoms": [0, 1]},
                {"bond_type": "single", "endpoint_atoms": [1, 2]},
                {"bond_type": "single", "endpoint_atoms": [1, 3]},
            ],
        }
    )
    s2, _, w2 = mod._molscribe_graph_to_smiles(
        {
            "atoms": [
                {"atom_symbol": "C"},
                {"atom_symbol": "[C@@H]"},
                {"atom_symbol": "OH"},
                {"atom_symbol": "[CF3]"},
            ],
            "bonds": [
                {"bond_type": "single", "endpoint_atoms": [0, 1]},
                {"bond_type": "single", "endpoint_atoms": [1, 2]},
                {"bond_type": "single", "endpoint_atoms": [1, 3]},
            ],
        }
    )
    assert w1 == [] and w2 == []
    assert kit.compare(s1, s2)["relation"] == "enantiomers"
    assert kit.compare(s1, "C[C@H](O)C(F)(F)F")["relation"] == "identical"

    # MEMOCH2 has one more carbon than MEMO
    memo, _, _ = mod._molscribe_graph_to_smiles(
        {"atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[MEMO]"}], "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}]}
    )
    memoch2, _, _ = mod._molscribe_graph_to_smiles(
        {"atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[MEMOCH2]"}], "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}]}
    )
    assert memo != memoch2


def test_route_figure_precatalyst_subscripts_cjk_wrapping_notes_and_dimerization():
    assert subscript_formulas("XPhos-Pd-G3, K3PO4") == "XPhos-Pd-G3, K₃PO₄"
    assert subscript_formulas("Hoveyda-HG2, CH2Cl2") == "Hoveyda-HG2, CH₂Cl₂"

    # Wide CJK characters are wider than ASCII and wrap without spaces
    assert _text_width("对甲苯磺酸", 14) > _text_width("abcde", 14) * 1.5
    wrapped = _wrap_caption("对甲苯磺酸吡啶盐催化的缩酮保护反应中间体", 13, 120)
    assert len(wrapped) >= 2

    # row["note"] is rendered into the SVG
    fig = render_route_figure(
        {
            "rows": [
                {
                    "note": "Fragment A preparation:",
                    "items": [{"label": "1", "smiles": "CCO"}, {"type": "arrow", "label": "DMP"}, {"label": "2", "smiles": "CC=O"}],
                }
            ]
        }
    )
    assert fig["status"] == "ok"
    assert "Fragment A preparation:" in fig["svg"]

    # Homodimerization preserves second equivalent after '+'
    dimer_spec = kit.route_to_figure_spec(
        {
            "compounds": {"1": "C=CCc1ccccc1", "2": "c1ccc(CC=CCCc2ccccc2)cc1"},
            "steps": [
                {"step": 1, "reactants": ["1"], "product": "1", "conditions": "purify"},
                {"step": 2, "reactants": ["1", "1"], "product": "2", "conditions": "Grubbs-II, CH2Cl2"},
            ],
        }
    )
    plus_items = [it for row in dimer_spec["rows"] for it in row["items"] if it.get("type") == "plus"]
    assert len(plus_items) == 1


def test_mechanism_teaching_annotations_diels_alder_and_cdxml_compliance():
    # Teaching mode renders partial_charges and extra visible lone_pairs, and CDXML uses integer Element
    teaching = {
        "spec_version": "2.1",
        "presentation_mode": "teaching",
        "panels": [
            {
                "title": "SN2 teaching",
                "reaction_arrow_after": 2,
                "starting_state": {"molecule_indices": [1, 2]},
                "expected_state": {"molecule_indices": [3, 4]},
                "molecules": [
                    {"smiles": "[OH-:1]"},
                    {"smiles": "[CH3:2][Br:3]"},
                    {"smiles": "[OH:1][CH3:2]"},
                    {"smiles": "[Br-:3]"},
                ],
                "lone_pairs": [{"atom_map": 1, "molecule_index": 1, "count": 3, "visible_count": 2}],
                "partial_charges": [{"atom_map": 2, "molecule_index": 2, "label": "δ+"}],
                "arrows": [
                    {"from_lone_pair_atom_map": 1, "to_atom_map": 2, "kind": "electron_pair"},
                    {"from_bond": [2, 3], "to_atom_map": 3, "kind": "electron_pair"},
                ],
                "graph_edits": [
                    {"type": "form_bond", "atom_maps": [1, 2], "order": 1, "source_arrow": 1},
                    {"type": "break_bond", "atom_maps": [2, 3], "source_arrow": 2},
                    {"type": "formal_charge", "atom_map": 1, "from": -1, "to": 0},
                    {"type": "formal_charge", "atom_map": 3, "from": 0, "to": -1},
                ],
            }
        ],
    }
    res = chem_mechanism_render(teaching)
    assert res["status"] == "ok"
    assert res["svg"].count("class='lone-pair-object'") == 2
    assert "δ+" in res["svg"]
    cdxml_root = ET.fromstring(res["cdxml"])
    nodes = list(cdxml_root.iter("n"))
    assert any(n.get("Element") == "8" for n in nodes)
    assert not any(n.get("Element") == "C" or n.get("Element") == "O" for n in nodes)
    arrows = list(cdxml_root.iter("arrow"))
    assert arrows and arrows[0].get("BoundingBox") and arrows[0].get("Head3D")
    rsteps = list(cdxml_root.iter("reactionstep"))
    assert rsteps and rsteps[0].get("ReactionStepProducts")


def test_review_session_memory_pruning_clears_old_terminal_assets():
    sess = sr.create_review_session(items=[{"id": "1", "smiles": "CCO"}])
    sid = sess["session_id"]
    with sr._SESSIONS_LOCK:
        obj = sr._SESSIONS[sid]
        obj.assets["image-1.png"] = ("image/png", b"x" * 1024)
        obj.status = "completed"
        obj.completed_at = time.time() - 1200
        sr._prune_sessions_locked()
        assert obj.assets == {}
