from codex_organic_chem.service import chem_compute


def _record(result, method):
    return next(item for item in result["records"] if item["method"] == method)


def test_descriptor_record_for_ethanol():
    result = chem_compute("CCO", tasks=["descriptors"])
    rec = _record(result, "rdkit_descriptors")
    assert rec["status"] == "ok"
    assert rec["results"]["formula"] == "C2H6O"
    assert rec["results"]["h_bond_donors"] == 1


def test_charges_record_for_ethanol():
    result = chem_compute("CCO", tasks=["charges"])
    rec = _record(result, "rdkit_gasteiger_charges")
    assert rec["status"] == "ok"
    assert rec["results"]["atom_charges"]


def test_xtb_unavailable_or_ok_is_structured():
    result = chem_compute("CCO", tasks=["xtb_opt"])
    rec = _record(result, "xtb_opt")
    assert rec["status"] in {"unavailable", "ok", "error"}
    assert rec["method"] == "xtb_opt"


def test_crest_record_is_structured():
    result = chem_compute("CCO", tasks=["crest"])
    rec = _record(result, "crest_conformer_search")
    assert rec["status"] in {"unavailable", "available"}
    assert rec["method"] == "crest_conformer_search"
