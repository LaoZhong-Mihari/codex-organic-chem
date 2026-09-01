import shlex
import sys

from codex_organic_chem.scheme_ocsr import _has_low_confidence_label_artifact_atoms, candidate_score
from codex_organic_chem.service import chem_parse_image


def _adapter(tmp_path, name, line):
    script = tmp_path / f"{name}.py"
    script.write_text(f"print({line!r})\n", encoding="utf-8")
    return f"{shlex.quote(sys.executable)} {shlex.quote(str(script))}"


def test_multi_engine_agreement_outranks_lone_engine(tmp_path, monkeypatch):
    image = tmp_path / "crop.png"
    image.write_bytes(b"plumbing only")
    # Two engines agree on ibuprofen (spelled differently); one alone disagrees.
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", _adapter(tmp_path, "a", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"))
    monkeypatch.setenv("CODEX_CHEM_DECIMER_CMD", _adapter(tmp_path, "b", "OC(=O)C(C)c1ccc(CC(C)C)cc1"))
    monkeypatch.setenv("CODEX_CHEM_MOLGRAPHER_CMD", _adapter(tmp_path, "c", "CC(C)Cc1ccc(CC(=O)O)cc1"))

    result = chem_parse_image(str(image), kind="molecule")

    top = result["candidates"][0]
    assert top["canonical_smiles"] == "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
    assert top["metadata"]["tool_agreement"] == 2
    assert top["metadata"]["consensus"] == "multi_engine"


def test_single_uncorroborated_read_warns(tmp_path, monkeypatch):
    image = tmp_path / "crop.png"
    image.write_bytes(b"plumbing only")
    monkeypatch.setenv("CODEX_CHEM_MOLSCRIBE_CMD", _adapter(tmp_path, "a", "CCO"))
    monkeypatch.setenv("CODEX_CHEM_DECIMER_CMD", _adapter(tmp_path, "b", "CCN"))

    result = chem_parse_image(str(image), kind="molecule")

    assert result["candidates"][0]["metadata"]["consensus"] in {"single_read", "multi_variant_single_engine"}
    if result["candidates"][0]["metadata"]["consensus"] == "single_read":
        assert any("uncorroborated" in warning for warning in result["warnings"])


def test_halogens_are_not_label_artifacts():
    assert not _has_low_confidence_label_artifact_atoms("Clc1ccc(Br)cc1F")
    assert not _has_low_confidence_label_artifact_atoms("Ic1ccccc1")
    assert _has_low_confidence_label_artifact_atoms("[Pd]")


def test_tool_agreement_raises_candidate_score():
    base = {"canonical_smiles": "CCO", "confidence": 0.8, "metadata": {"tool_agreement": 1, "variant_agreement": 1}}
    corroborated = {
        "canonical_smiles": "CCO",
        "confidence": 0.8,
        "metadata": {"tool_agreement": 3, "variant_agreement": 2},
    }
    assert candidate_score(corroborated) > candidate_score(base)


def test_missing_image_keeps_full_result_schema():
    result = chem_parse_image("/nonexistent/path/image.png", kind="molecule")
    assert result["status"] == "no_candidates"
    assert result["candidates"] == []
    assert result["ranked_candidates"] == []
    assert result["confirmation_required"] is False
    assert result["warnings"]
