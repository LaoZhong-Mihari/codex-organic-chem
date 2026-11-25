from __future__ import annotations

import json
from pathlib import Path
from PIL import Image

from codex_organic_chem import kit, service
from codex_organic_chem.cli import main as cli_main


def test_mechanism_steps_check_catches_skeleton_change_radicals_and_unfinished_cascade():
    # Clean Michael addition to Wieland-Miescher precursor
    clean = kit.mechanism_steps_check(
        reactants="C=CCC1C(=O)CCCC1=O.CC(=O)C=C",
        products="CC(=O)CCC1(CC=C)C(=O)CCCC1=O",
        conditions="NaOH, H2O, rt",
        steps=[
            {"step": 1, "type": "proton_transfer", "subtype": "acid_base_proton_transfer", "intermediate_smiles": "C=CCC1=C([O-])CCCC1=O"},
            {"step": 2, "type": "addition", "subtype": "nucleophilic_addition", "intermediate_smiles": "C=CCC1(CC=C(C)[O-])C(=O)CCCC1=O"},
            {"step": 3, "type": "proton_transfer", "subtype": "acid_base_proton_transfer", "intermediate_smiles": "CC(=O)CCC1(CC=C)C(=O)CCCC1=O"},
        ],
    )
    assert clean["status"] == "consistent"
    assert clean["n_blocking"] == 0
    assert clean["final_vs_product"]["relation"] == "identical"

    # ART-005 failure mode: stopping at the initial Michael adduct when the target is a rearranged isomer
    cascade_miss = kit.mechanism_steps_check(
        reactants="C=CCC1C(=O)CCCC1=O.CC(=O)C=C",
        products="C(CCC1C(=O)C(CCC1=O)CC=C)(=O)C",
        conditions="NaOH, H2O, rt",
        steps=[
            {"step": 1, "type": "proton_transfer", "subtype": "acid_base_proton_transfer", "intermediate_smiles": "C=CCC1=C([O-])CCCC1=O"},
            {"step": 2, "type": "addition", "subtype": "nucleophilic_addition", "intermediate_smiles": "C=CCC1(CC=C(C)[O-])C(=O)CCCC1=O"},
            {"step": 3, "type": "proton_transfer", "subtype": "acid_base_proton_transfer", "intermediate_smiles": "CC(=O)CCC1(CC=C)C(=O)CCCC1=O"},
        ],
    )
    assert cascade_miss["status"] == "blocked"
    codes = {i["code"] for i in cascade_miss["issues"]}
    assert "incomplete_cascade_constitutional_isomer" in codes

    # Bogert-Cook R4 failure mode: proton transfer step moves -OH to adjacent carbon
    bad_pt = kit.mechanism_steps_check(
        reactants="CC1(C)CCCC1O",
        products="CC1=C(C)CCC1",
        conditions="H3PO4, heat",
        steps=[
            {"step": 1, "type": "proton_transfer", "subtype": "acid_base_proton_transfer", "intermediate_smiles": "CC1(C)CCC([OH2+])C1"},
        ],
    )
    assert any(i["code"] == "proton_transfer_changes_skeleton" for i in bad_pt["issues"])

    # Polar mechanism accidentally writing a radical atom instead of a cation
    bad_rad = kit.mechanism_steps_check(
        reactants="CC(C)(O)C",
        products="CC(=C)C",
        conditions="H2SO4",
        steps=[
            {"step": 1, "type": "elimination", "subtype": "leaving_group_elimination", "intermediate_smiles": "C[C](C)C"},
        ],
    )
    assert any(i["code"] == "unexpected_radical_intermediate" for i in bad_rad["issues"])


def test_crop_image_and_parse_scheme_inline_boxes(tmp_path):
    page = tmp_path / "exam_page.png"
    img = Image.new("RGB", (400, 300), (255, 255, 255))
    img.save(page)

    grid_only = kit.crop_image(page, output_dir=tmp_path / "out_grid")
    assert grid_only["status"] == "ok"
    assert Path(grid_only["grid_image"]).is_file()
    assert grid_only["source_size"] == {"width": 400, "height": 300}

    boxes_json = json.dumps([{"id": "A", "box": [0.1, 0.2, 0.4, 0.6]}, {"id": "B", "box": [200, 50, 350, 200]}])
    cropped = service.chem_crop_image(str(page), boxes=boxes_json, output_dir=str(tmp_path / "crops"), scale=2.0, grid=True)
    assert cropped["status"] == "ok"
    assert len(cropped["crops"]) == 2
    assert all(Path(c["path"]).is_file() for c in cropped["crops"])

    # Passing JSON bounding boxes directly to chem_parse_scheme no longer raises Errno 63 File name too long
    long_boxes = json.dumps([{"id": f"opt_{i}", "box": [10 * i, 10, 40 + 10 * i, 60]} for i in range(1, 12)])
    scheme_res = service.chem_parse_scheme(str(page), crops=long_boxes, include_drawings=False)
    assert scheme_res["status"] == "ok"
    assert len(scheme_res["compounds"]) == 11
    assert [c["compound"] for c in scheme_res["compounds"]] == [str(i) for i in range(1, 12)]


def test_reaction_analyze_mass_charge_balance_alias_and_cli(capsys):
    res = service.chem_reaction_analyze("CBr.[OH-]>>CO.[Br-]", mode="mass_charge_balance", include_drawings=False)
    assert res["analysis"]["mass_charge_balance"]["status"] == "balanced"
    assert "<svg" not in json.dumps(res) and "M  END" not in json.dumps(res)

    cli_main(
        [
            "mechanism-steps-check",
            "--reaction",
            "C=O.[C-]#N.C#N>>OCC#N.[C-]#N",
            "--steps",
            '["[O-]CC#N", "OCC#N"]',
            "--issues-only",
        ]
    )
    out = json.loads(capsys.readouterr().out)
    assert out["status"] == "consistent" and "steps" not in out


def test_review_off_and_ask_modes_across_legacy_and_review_tools(monkeypatch):
    monkeypatch.setenv("CODEX_CHEM_REVIEW", "off")
    off_input = service.chem_input_review(smiles="CCO")
    assert off_input["status"] == "review_disabled"
    assert off_input["confirmation_required"] is False
    assert off_input["blocked_until_confirmed"] == []

    off_synth = service.chem_synthesis_suggest("CC(=O)Oc1ccccc1C(=O)O", literature=False)
    assert off_synth["status"] == "awaiting_route_decision"

    off_mech = service.chem_mechanism_draft("CBr.[OH-]>>CO.[Br-]", quality="publication")
    assert off_mech["status"] == "ok"

    monkeypatch.setenv("CODEX_CHEM_REVIEW", "ask")
    ask_input = service.chem_input_review(smiles="CCO")
    assert ask_input["status"] == "ask_user_first"
    assert "Open it?" in ask_input["user_prompt"]


def test_every_mcp_tool_has_docstring_and_skill_resources_exist():
    import asyncio
    from unittest.mock import patch
    from codex_organic_chem import mcp_server

    captured = {}

    def _capture_run(self, *args, **kwargs):
        captured["server"] = self

    with patch("mcp.server.fastmcp.FastMCP.run", _capture_run):
        mcp_server.main()

    server = captured["server"]
    tools = asyncio.run(server.list_tools())
    assert len(tools) >= 28
    missing_docs = [t.name for t in tools if not (t.description and t.description.strip())]
    assert missing_docs == [], f"MCP tools missing descriptions: {missing_docs}"

    resources = asyncio.run(server.list_resources())
    uris = {str(r.uri) for r in resources}
    assert "chem://skill/SKILL.md" in uris
    assert "chem://skill/references/kit-api.md" in uris
