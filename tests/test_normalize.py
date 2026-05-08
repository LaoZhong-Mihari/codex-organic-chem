import importlib.util
import shlex
import sys
from pathlib import Path

from codex_organic_chem.service import chem_draw, chem_input_review, chem_normalize_structure, chem_parse_image


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
                "print(json.dumps({'smiles': 'CCO'}))",
            ]
        ),
        encoding="utf-8",
    )
    command = f"{shlex.quote(sys.executable)} {shlex.quote(str(adapter))} --input {{input}}"
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", command)

    result = chem_parse_image(str(image), kind="molecule")

    assert result["status"] == "awaiting_user_confirmation"
    assert result["candidates"][0]["canonical_smiles"] == "CCO"


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


def test_input_review_blocks_until_confirmation():
    result = chem_input_review(smiles="CCO")
    assert result["status"] == "awaiting_user_confirmation"
    assert result["confirmation_required"] is True
    assert "svg" in result["preview"]
