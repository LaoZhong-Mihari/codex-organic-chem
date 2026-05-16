import importlib.util
import shlex
import sys
from pathlib import Path

from codex_organic_chem.scheme_ocsr import rank_candidates, resolve_placeholders, route_consistency_warnings
from codex_organic_chem.ocsr import _parse_candidate_lines
from codex_organic_chem.service import (
    chem_draw,
    chem_input_review,
    chem_normalize_structure,
    chem_ocsr_benchmark,
    chem_parse_image,
    chem_parse_scheme,
)


def test_normalize_aspirin():
    result = chem_normalize_structure(smiles="CC(=O)Oc1ccccc1C(=O)O")
    assert result["canonical_smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert result["metadata"]["formula"] == "C9H8O4"
    assert result["confidence"] > 0.5


def test_invalid_smiles_returns_warning():
    result = chem_normalize_structure(smiles="not-a-smiles")
    assert result["confidence"] == 0.0
    assert result["warnings"]


def test_draw_svg():
    result = chem_draw(smiles="c1ccccc1C(=O)O", output="svg")
    assert result["status"] == "ok"
    assert "<svg" in result["data"]


def test_missing_image_returns_structured_error(tmp_path):
    result = chem_parse_image(str(tmp_path / "missing.png"))
    assert result["candidates"] == []
    assert result["warnings"]


def test_custom_ocsr_adapter_handles_image_paths_with_spaces(tmp_path, monkeypatch):
    image = tmp_path / "image with spaces.png"
    image.write_bytes(b"not a real png; adapter only checks path plumbing")
    adapter = tmp_path / "adapter.py"
    adapter.write_text(
        "\n".join(
            [
                "import argparse, json, pathlib",
                "parser = argparse.ArgumentParser()",
                "parser.add_argument('--input', required=True)",
                "args = parser.parse_args()",
                "assert pathlib.Path(args.input).exists()",
                "print(json.dumps({'smiles': 'CCO', 'adapter_warnings': ['adapter note']}))",
            ]
        ),
        encoding="utf-8",
    )
    command = f"{shlex.quote(sys.executable)} {shlex.quote(str(adapter))} --input {{input}}"
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", command)

    result = chem_parse_image(str(image), kind="molecule")

    assert result["status"] == "awaiting_user_confirmation"
    assert result["candidates"][0]["canonical_smiles"] == "CCO"
    assert "adapter note" in result["warnings"]


def test_json_adapter_warning_without_candidates_is_not_smiles():
    candidates, warnings = _parse_candidate_lines(
        '{"tool": "chemschematicresolver", "candidates": [], "adapter_warnings": ["no candidates"]}',
        "chemschematicresolver",
    )

    assert candidates == []
    assert warnings == ["no candidates"]


def test_multi_adapter_ensemble_prefers_sanitized_low_dummy_candidate(tmp_path, monkeypatch):
    image = tmp_path / "scheme crop.png"
    image.write_bytes(b"adapter only checks path plumbing")
    bad_adapter = tmp_path / "bad_adapter.py"
    good_adapter = tmp_path / "good_adapter.py"
    bad_adapter.write_text("print('*.*')\n", encoding="utf-8")
    good_adapter.write_text("print('CCO')\n", encoding="utf-8")
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(bad_adapter))}")
    monkeypatch.setenv("CODEX_CHEM_DECIMER_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(good_adapter))}")

    result = chem_parse_image(str(image), kind="molecule")

    assert result["candidates"][0]["canonical_smiles"] == "CCO"
    assert result["candidates"][0]["metadata"]["ocsr_tool"] == "decimer"
    assert result["ranked_candidates"][0]["metadata"]["rank"] == 1


def test_multi_adapter_ensemble_penalizes_collapsed_tiny_candidates(tmp_path, monkeypatch):
    image = tmp_path / "scheme crop.png"
    image.write_bytes(b"adapter only checks path plumbing")
    collapsed_adapter = tmp_path / "collapsed_adapter.py"
    larger_adapter = tmp_path / "larger_adapter.py"
    collapsed_adapter.write_text("print('C')\n", encoding="utf-8")
    larger_adapter.write_text("print('*OCCC1(C)CCC(C(C)C)C1=C.I')\n", encoding="utf-8")
    monkeypatch.setenv("CODEX_CHEM_MOLGRAPHER_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(collapsed_adapter))}")
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(larger_adapter))}")

    result = chem_parse_image(str(image), kind="molecule")

    assert result["candidates"][0]["canonical_smiles"] != "C"
    assert result["candidates"][0]["metadata"]["ocsr_tool"] == "molscribe"


def test_multi_adapter_ensemble_penalizes_low_confidence_label_artifacts(tmp_path, monkeypatch):
    image = tmp_path / "scheme crop.png"
    image.write_bytes(b"adapter only checks path plumbing")
    artifact_adapter = tmp_path / "artifact_adapter.py"
    placeholder_adapter = tmp_path / "placeholder_adapter.py"
    artifact_adapter.write_text(
        "import json\nprint(json.dumps({'tool': 'molgrapher', 'smiles': 'CC(F)=CB', 'confidence': 0.0}))\n",
        encoding="utf-8",
    )
    placeholder_adapter.write_text(
        "import json\nprint(json.dumps({'tool': 'molscribe', 'smiles': '*C=C*', 'confidence': 0.85}))\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("CODEX_CHEM_MOLGRAPHER_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(artifact_adapter))}")
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(placeholder_adapter))}")

    result = chem_parse_image(str(image), kind="molecule")

    assert result["candidates"][0]["metadata"]["ocsr_tool"] == "molscribe"


def test_ocsr_preprocessing_recovers_faint_skeletal_line_art(tmp_path, monkeypatch):
    from PIL import Image, ImageDraw

    image = tmp_path / "faint skeletal crop.png"
    canvas = Image.new("RGB", (120, 70), "white")
    draw = ImageDraw.Draw(canvas)
    draw.line([(12, 45), (38, 25), (64, 45), (90, 25), (110, 36)], fill=(218, 218, 218), width=1)
    draw.line([(64, 45), (64, 62)], fill=(218, 218, 218), width=1)
    canvas.save(image)
    raw_black_pixels = sum(canvas.convert("L").histogram()[:80])
    assert raw_black_pixels == 0
    assert max(canvas.size) < 500

    adapter = tmp_path / "contrast_sensitive_adapter.py"
    adapter.write_text(
        "\n".join(
            [
                "import argparse, json",
                "from PIL import Image",
                "parser = argparse.ArgumentParser()",
                "parser.add_argument('--input', required=True)",
                "args = parser.parse_args()",
                "img = Image.open(args.input).convert('L')",
                "black_pixels = sum(img.histogram()[:80])",
                "if black_pixels > 80 and max(img.size) >= 500:",
                "    print(json.dumps({'tool': 'contrast-sensitive', 'smiles': 'CC(C)CO', 'confidence': 0.86}))",
                "else:",
                "    print(json.dumps({'tool': 'contrast-sensitive', 'smiles': 'C', 'confidence': 0.25}))",
            ]
        ),
        encoding="utf-8",
    )
    for env_var in (
        "CODEX_CHEM_MOLSCRIBE_CMD",
        "CODEX_CHEM_DECIMER_CMD",
        "CODEX_CHEM_MOLGRAPHER_CMD",
        "CODEX_CHEM_OPENCHEMIE_CMD",
        "CODEX_CHEM_CSR_CMD",
        "CODEX_CHEM_RXNSCRIBE_CMD",
    ):
        monkeypatch.delenv(env_var, raising=False)
    monkeypatch.setenv("CODEX_CHEM_DECIMER_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(adapter))} --input {{input}}")

    result = chem_parse_image(str(image), kind="molecule")

    assert result["candidates"][0]["canonical_smiles"] == "CC(C)CO"
    assert result["candidates"][0]["metadata"]["ocsr_tool"] == "contrast-sensitive"
    assert result["candidates"][0]["metadata"]["image_variant"] in {"high_contrast_binary", "high_contrast_thick"}
    assert any("line-art variants" in warning for warning in result["warnings"])


def test_candidate_ranking_prefers_original_crop_for_stereo_tie():
    original = {
        "canonical_smiles": "[H]C(C)(OS(C)(=O)=O)C(=O)c1ccc(C(C)(C)C)cc1",
        "isomeric_smiles": "[H][C@@](C)(OS(C)(=O)=O)C(=O)c1ccc(C(C)(C)C)cc1",
        "confidence": 0.868,
        "metadata": {"image_variant": "original"},
        "warnings": [],
    }
    preprocessed = {
        "canonical_smiles": "[H]C(C)(OS(C)(=O)=O)C(=O)c1ccc(C(C)(C)C)cc1",
        "isomeric_smiles": "[H][C@](C)(OS(C)(=O)=O)C(=O)c1ccc(C(C)(C)C)cc1",
        "confidence": 0.879,
        "metadata": {"image_variant": "white_padded"},
        "warnings": [],
    }

    ranked = rank_candidates([preprocessed, original])

    assert ranked[0]["metadata"]["image_variant"] == "original"
    assert ranked[0]["isomeric_smiles"] == original["isomeric_smiles"]


def test_ocsr_benchmark_reports_exact_match_with_fake_adapter(tmp_path, monkeypatch):
    image_dir = tmp_path / "crops"
    image_dir.mkdir()
    (image_dir / "A_source_7.png").write_bytes(b"fake crop")
    gold = tmp_path / "gold.tsv"
    gold.write_text("label\tcompound\tsmiles\nA\t7\tCCO\n", encoding="utf-8")
    adapter = tmp_path / "adapter.py"
    adapter.write_text("print('CCO')\n", encoding="utf-8")
    monkeypatch.setenv("CODEX_CHEM_DECIMER_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(adapter))}")

    result = chem_ocsr_benchmark(str(gold), str(image_dir))

    assert result["summary"]["gold_count"] == 1
    assert result["summary"]["exact_graph_matches"] == 1
    assert result["rows"][0]["best"]["tool"] == "decimer"


def test_parse_scheme_resolves_gold_map_placeholders(tmp_path, monkeypatch):
    image = tmp_path / "scheme.png"
    image.write_bytes(b"fake scheme")
    crop_dir = tmp_path / "crops"
    crop_dir.mkdir()
    (crop_dir / "C_source_11.png").write_bytes(b"fake crop")
    legend = tmp_path / "legend.json"
    legend.write_text('{"11": {"1": "OTs"}}', encoding="utf-8")
    adapter = tmp_path / "adapter.py"
    adapter.write_text("print('[*:1]CCOCOCCOC')\n", encoding="utf-8")
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", f"{shlex.quote(sys.executable)} {shlex.quote(str(adapter))}")

    result = chem_parse_scheme(str(image), str(crop_dir), gold_map=str(legend))

    resolved = result["compounds"][0]["resolved_candidates"][0]
    assert resolved["status"] == "resolved"
    assert "S(=O)(=O)" in resolved["smiles"]
    assert "*" not in resolved["smiles"]


def test_molscribe_adapter_expands_common_abbreviations_and_chirality():
    adapter_path = Path(__file__).resolve().parents[1] / "scripts" / "molscribe_adapter.py"
    spec = importlib.util.spec_from_file_location("molscribe_adapter", adapter_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    smiles, molfile, warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [
                {"atom_symbol": "[Me]"},
                {"atom_symbol": "[C@H]"},
                {"atom_symbol": "O"},
                {"atom_symbol": "[CO2Et]"},
            ],
            "bonds": [
                {"bond_type": "single", "endpoint_atoms": [0, 1]},
                {"bond_type": "single", "endpoint_atoms": [1, 2]},
                {"bond_type": "solid wedge", "endpoint_atoms": [1, 3]},
            ],
        }
    )

    assert molfile
    assert warnings == []
    assert "@" in smiles
    assert "C(=O)OCC" in smiles or "CCOC(=O)" in smiles


def test_molscribe_adapter_expands_sulfonyl_leaving_groups_and_protecting_groups():
    adapter_path = Path(__file__).resolve().parents[1] / "scripts" / "molscribe_adapter.py"
    spec = importlib.util.spec_from_file_location("molscribe_adapter", adapter_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    mesylate_smiles, _, mesylate_warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[OMs]"}],
            "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}],
        }
    )
    tosylhydrazide_smiles, _, tosylhydrazide_warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [{"atom_symbol": "N"}, {"atom_symbol": "[NHTs]"}],
            "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}],
        }
    )

    assert mesylate_warnings == []
    assert "S(C)(=O)=O" in mesylate_smiles
    assert mesylate_smiles.startswith("CO")
    assert tosylhydrazide_warnings == []
    assert "Cc1ccc(S(=O)(=O)NN)" in tosylhydrazide_smiles


def test_molscribe_adapter_expands_scheme_abbreviations_and_preserves_placeholders():
    adapter_path = Path(__file__).resolve().parents[1] / "scripts" / "molscribe_adapter.py"
    spec = importlib.util.spec_from_file_location("molscribe_adapter", adapter_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    memo_smiles, _, memo_warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[MEMO]"}],
            "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}],
        }
    )
    sulfonyl_smiles, _, sulfonyl_warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[SO2Tol]"}],
            "bonds": [{"bond_type": "single", "endpoint_atoms": [0, 1]}],
        }
    )
    placeholder_smiles, _, placeholder_warnings = module._molscribe_graph_to_smiles(
        {
            "atoms": [{"atom_symbol": "C"}, {"atom_symbol": "[RO]"}, {"atom_symbol": "[X]"}],
            "bonds": [
                {"bond_type": "single", "endpoint_atoms": [0, 1]},
                {"bond_type": "single", "endpoint_atoms": [0, 2]},
            ],
        }
    )

    assert memo_warnings == []
    assert "COCCOCOC" in memo_smiles
    assert sulfonyl_warnings == []
    assert "Cc1ccc(S(C)(=O)=O)" in sulfonyl_smiles
    assert placeholder_warnings == []
    assert placeholder_smiles.count("*") == 2


def test_placeholder_resolver_expands_leaving_group_and_memo():
    result = resolve_placeholders("CC([*:1])CCOCOCCOC", {"1": "OTs"})

    assert result["status"] == "resolved"
    assert "S(=O)(=O)" in result["smiles"]
    assert "c1ccc(C)cc1" in result["smiles"]
    assert result["unresolved_labels"] == []


def test_placeholder_resolver_replaces_h_without_dummy():
    result = resolve_placeholders("CC1([*:1])CCCCC1", {"1": "H"})

    assert result["status"] == "resolved"
    assert "*" not in result["smiles"]


def test_placeholder_resolver_distinguishes_13_14_16_definitions():
    tosyl = resolve_placeholders("CC([*:1])", {"1": "SO2Tol"})["smiles"]
    carbonyl = resolve_placeholders("CC(=[*:1])", {"1": "O"})["smiles"]
    aldehyde = resolve_placeholders("CC([*:1])", {"1": "CHO"})["smiles"]

    assert tosyl != carbonyl != aldehyde
    assert "S(=O)(=O)" in tosyl
    assert carbonyl == "CC=O"
    assert "C=O" in aldehyde or "CC(C)=O" in aldehyde


def test_route_consistency_flags_aromatic_relation_jump():
    warnings = route_consistency_warnings(
        [
            {"label": "A", "smiles": "CCc1ccc(CCC)cc1"},
            {"label": "B", "smiles": "CCc1cccc(CCC)c1"},
        ]
    )

    assert any("aromatic substitution relation changed" in warning for warning in warnings)


def test_input_review_blocks_until_confirmation():
    result = chem_input_review(smiles="CCO")
    assert result["status"] == "awaiting_user_confirmation"
    assert result["confirmation_required"] is True
    assert "svg" in result["preview"]
