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


def test_conformer_record_for_butanol():
    result = chem_compute("CCCCO", tasks=["conformers"], num_confs=4)
    rec = _record(result, "rdkit_conformers")
    assert rec["status"] == "ok"


# xTB and CREST records are covered in tests/test_semiempirical.py, which asserts
# on the returned numbers. The status-membership assertions that used to live here
# were satisfied by a CREST stub that never executed anything.
