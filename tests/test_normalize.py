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


def test_input_review_blocks_until_confirmation():
    result = chem_input_review(smiles="CCO")
    assert result["status"] == "awaiting_user_confirmation"
    assert result["confirmation_required"] is True
    assert "svg" in result["preview"]
