import pytest

from codex_organic_chem import kit
from codex_organic_chem import model_kit as mk
from codex_organic_chem import service
from codex_organic_chem.literature import _abstract_from_inverted_index, _normalize_doi, literature_search
from codex_organic_chem.utils import executable_status

CAMPHOR = "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2"  # carbonyl carbon is atom 8
NORCAMPHOR = "O=C1C[C@@H]2CC[C@H]1C2"  # carbonyl carbon is atom 1


def test_norcamphor_carbonyl_is_open_on_the_exo_face():
    result = kit.face(NORCAMPHOR, 1)
    open_face = result["more_open_face"]
    shielded_face = "Re" if open_face == "Si" else "Si"
    assert result["shielding_difference"] > 0.1
    # endo face is blocked by the two-carbon bridge (atoms 4, 5); exo only sees C7
    assert {a["heavy_atom_index"] for a in result["shielding_atoms"][shielded_face]} == {4, 5}
    assert [a["heavy_atom_index"] for a in result["shielding_atoms"][open_face]] == [7]


def test_camphor_exo_face_is_blocked_by_the_bridge_methyl():
    result = kit.face(CAMPHOR, 8)
    open_face = result["more_open_face"]
    assert open_face is not None
    shielded_face = "Re" if open_face == "Si" else "Si"
    # atom 0 is a methyl on C7, which hangs over the exo face of the carbonyl
    assert result["shielding_atoms"][shielded_face][0]["heavy_atom_index"] == 0


def test_symmetric_ketone_has_no_facial_bias_and_sp3_centres_are_rejected():
    acetone = kit.face("CC(=O)C", 1)
    assert acetone["more_open_face"] is None and "Re/Si is undefined" in acetone["face_definition"]
    with pytest.raises(kit.StructureError):
        kit.face("CCO", 1)
    assert mk.model_kit("face", "CCO", atoms=[1])["status"] == "error"


def test_dihedral_populations():
    butane = kit.dihedral("CCCC", [0, 1, 2, 3])
    assert butane["anti_periplanar"]["population"] > butane["gauche"]["population"] > 0.1
    dibromide = kit.dihedral("Br[C@@H]1CCCC[C@H]1Br", [0, 1, 6, 7])
    assert dibromide["anti_periplanar"]["population"] > 0.5  # diaxial conformer is reachable


def test_distance_reports_closest_approach_and_its_cost():
    result = kit.distance("C=CCCCCCC=C", 0, 8, threshold=4.0)
    assert result["distance_in_lowest_conformer"] > result["closest_approach"]["distance"]
    assert result["closest_approach"]["relative_energy_kcal"] >= 0
    rigid = kit.distance("c1ccccc1C#Cc1ccccc1", 0, 13)
    assert rigid["population_within_threshold"] == 0 and "No sampled conformer" in rigid["interpretation"]


def test_planarity_without_xtb():
    def label(smiles):
        report = kit.planarity(smiles, refine="none")
        return {s["selection"]: s["label"] for s in report["selections"]}

    assert label("c1ccccc1")["all_heavy_atoms"] == "planar"
    assert label("C1CCCCC1")["all_heavy_atoms"] == "nonplanar"
    assert label("C1=CC=CC=CC=C1")["all_heavy_atoms"] == "nonplanar"  # cyclooctatetraene tub
    biphenyl = label("c1ccccc1-c1ccccc1")
    assert biphenyl["ring_system_1"] == "planar" and biphenyl["all_heavy_atoms"] == "nonplanar"


def test_impossible_bridgehead_configuration_does_not_embed():
    assert kit.embed_check("C1C[C@H]2CC[C@@H]1C2")["embeds"] is True
    impossible = kit.embed_check("C1C[C@H]2CC[C@H]1C2")
    assert impossible["embeds"] is False and impossible["embeds_without_stereo"] is True


def test_isomer_energies_rank_known_pairs_and_refuse_non_isomers():
    decalin = kit.isomer_energies(["C1CC[C@H]2CCCC[C@@H]2C1", "C1CC[C@H]2CCCC[C@H]2C1"])
    assert decalin["most_stable"] == "C1CC[C@H]2CCCC[C@@H]2C1"  # trans
    assert decalin["gap_to_next_kcal"] > 1.0
    butene = kit.isomer_energies(["C/C=C/C", "C/C=C\\C"])
    assert butene["most_stable"] == "C/C=C/C" and butene["decisive"] is False
    with pytest.raises(kit.StructureError):
        kit.isomer_energies(["CCO", "CCN"])


@pytest.mark.skipif(
    executable_status("xtb", ("--version",)).status != "available", reason="xTB is not installed"
)
def test_xtb_refinement_keeps_the_known_ordering():
    decalin = kit.isomer_energies(
        ["C1CC[C@H]2CCCC[C@@H]2C1", "C1CC[C@H]2CCCC[C@H]2C1"], method="xtb"
    )
    assert decalin["most_stable"] == "C1CC[C@H]2CCCC[C@@H]2C1"
    assert all(row["method"].startswith("GFN2-xTB") for row in decalin["isomers"])
    assert kit.planarity("C1=CC=CC=CC=C1", refine="xtb")["geometry"].startswith("GFN2-xTB")


def test_undefined_stereo_is_reported_before_3d_conclusions():
    ensemble = kit.build_ensemble("CC(O)C(C)Cl", n_confs=4)
    assert any("undefined" in warning for warning in ensemble.warnings)
    assert kit.build_ensemble("C[C@H](O)[C@@H](C)Cl", n_confs=4).warnings == []


def test_atom_references_by_map_number():
    by_index = kit.distance("OCCCCO", 0, 5)
    by_map = kit.distance("[OH:7]CCCC[OH:9]", 7, 9, by_map=True)
    assert by_index["atoms"] == by_map["atoms"] == [0, 5]
    table = kit.atom_table("[OH:7]CC")
    assert table["atoms"][0]["map"] == 7 and table["atoms"][1]["neighbors"] == [0, 2]


def test_service_wrappers_return_errors_as_data():
    assert service.chem_compare_structures("C1CC")["status"] == "error"
    assert service.chem_enumerate("CCO", task="nonsense")["status"] == "error"
    assert service.chem_model_kit("distance", "CCO", atoms=[0])["status"] == "error"
    assert service.chem_enumerate("CC1CCCCC1", task="substitution")["n_constitutional_isomers"] == 5
    assert service.chem_compare_structures("CC(O)CC")["n_undefined"] == 1


def test_python_tool_is_off_unless_explicitly_enabled(monkeypatch):
    monkeypatch.delenv(service.PYTHON_ENV_FLAG, raising=False)
    assert service.chem_python("print(1)")["status"] == "disabled"
    monkeypatch.setenv(service.PYTHON_ENV_FLAG, "1")
    result = service.chem_python(
        "from codex_organic_chem import kit; print(kit.compare('CCO', 'COC')['relation'])"
    )
    assert result["status"] == "ok" and result["stdout"].strip() == "constitutional_isomers"
    assert service.chem_python("raise SystemExit(3)")["exit_code"] == 3
    assert service.chem_python("import time; time.sleep(5)", timeout_s=1)["status"] == "timeout"


def test_literature_helpers_and_offline_behaviour(monkeypatch):
    assert _normalize_doi("https://doi.org/10.1038/S41467-026-78679-3") == "10.1038/s41467-026-78679-3"
    assert _abstract_from_inverted_index({"selective": [1], "Highly": [0], "hydroboration": [2]}) == (
        "Highly selective hydroboration"
    )
    assert literature_search("x", source="scopus")["status"] == "error"

    import codex_organic_chem.literature as literature

    def offline(url, timeout_s):
        raise OSError("network disabled in test")

    monkeypatch.setattr(literature, "_get_json", offline)
    result = literature_search("alkyne hydroboration", rows=3)
    assert result["status"] == "error" and result["items"] == []
    assert len(result["warnings"]) == 2  # both sources reported, nothing invented

    def fake(url, timeout_s):
        if "openalex" in url:
            return {
                "results": [
                    {
                        "doi": "https://doi.org/10.1000/abc",
                        "title": "A",
                        "publication_year": 2024,
                        "abstract_inverted_index": {"Text": [0]},
                        "cited_by_count": 3,
                    }
                ]
            }
        return {"message": {"items": [{"DOI": "10.1000/ABC", "title": ["A"]}, {"DOI": "10.1000/xyz", "title": ["B"]}]}}

    monkeypatch.setattr(literature, "_get_json", fake)
    merged = literature_search("q", rows=5)
    assert [item["doi"] for item in merged["items"]] == ["10.1000/abc", "10.1000/xyz"]
    assert merged["items"][0]["abstract"] == "Text" and merged["source"] == "OpenAlex+Crossref"
    assert any("no abstract" in warning for warning in merged["warnings"])


def test_cis_trans_of_two_ring_substituents():
    # 1,2-dimethylcyclohexane: the achiral (meso) isomer is cis, the chiral pair is trans
    isomers = kit.enumerate_stereoisomers("CC1CCCCC1C")
    (meso,) = isomers["achiral_isomers"]
    assert kit.cis_trans(meso, 1, 6)["relationship"] == "cis"
    for chiral in isomers["enantiomer_pairs"][0]:
        assert kit.cis_trans(chiral, 1, 6)["relationship"] == "trans"
    # natural trans-4-hydroxy-L-proline is (2S,4R); its C4 epimer is cis
    trans = kit.cis_trans("OC(=O)[C@@H]1C[C@@H](O)CN1", 3, 5)
    assert trans["relationship"] == "trans" and trans["population_cis"] == 0.0
    assert kit.cis_trans("OC(=O)[C@@H]1C[C@H](O)CN1", 3, 5)["relationship"] == "cis"
    assert kit.cis_trans("C[C@H]1CC[C@@H](C)CC1", 1, 4)["relationship"] == "cis"
    with pytest.raises(kit.StructureError):
        kit.cis_trans("CC1CCCCC1", 1, 3)  # C3 is a CH2
    with pytest.raises(kit.StructureError):
        kit.cis_trans("CC(C)CC", 1, 2)  # not a ring
    assert mk.model_kit("cis_trans", "C[C@H]1CC[C@@H](C)CC1", atoms=[1, 4])["relationship"] == "cis"
    assert service.chem_model_kit("cis_trans", "C[C@H]1CC[C@H](C)CC1", atoms=[1, 4])["relationship"] == "trans"
