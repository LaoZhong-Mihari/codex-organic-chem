import random

import pytest
from rdkit import Chem

from codex_organic_chem import kit
from codex_organic_chem import structure_compare as sc


STEREO_MOLECULES = [
    "C[C@H](O)[C@@H](N)C(=O)O",
    "C[C@@H]1CC[C@H](C(C)C)C[C@H]1O",
    "O=C1C[C@@H]2CC[C@]1(C)C2(C)C",
    "C/C=C/[C@H](O)C/C=C\\C",
    "N[C@@H](Cc1ccccc1)C(=O)N[C@@H](C)C(=O)O",
    "C[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](CO)O1",
]


@pytest.mark.parametrize("smiles", STEREO_MOLECULES)
def test_identity_is_independent_of_how_the_smiles_is_written(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for _ in range(10):
        rewritten = Chem.MolToSmiles(mol, doRandom=True)
        assert sc.compare_structures(smiles, rewritten)["relation"] == "identical"


@pytest.mark.parametrize("smiles", STEREO_MOLECULES)
def test_single_inverted_centre_is_located_on_the_right_atom(smiles):
    mol = sc.parse(smiles)
    centres = [e["atom_index"] for e in sc.stereo_elements(mol) if e["kind"] == "tetrahedral"]
    for centre in centres:
        flipped = Chem.Mol(mol)
        flipped.GetAtomWithIdx(centre).InvertChirality()
        order = list(range(flipped.GetNumAtoms()))
        random.Random(centre).shuffle(order)
        scrambled = Chem.MolToSmiles(Chem.RenumberAtoms(flipped, order), canonical=False)
        result = sc.compare_structures(smiles, scrambled)
        if len(centres) == 1:
            # inverting the only centre gives the mirror image (E/Z is achiral)
            assert result["relation"] == "enantiomers"
            continue
        assert result["relation"] == "diastereomers"
        assert result["n_differing"] == 1
        assert result["differing_elements"][0]["atom_index_a"] == centre


def test_named_relationships():
    l_tartaric = "O[C@H](C(=O)O)[C@H](O)C(=O)O"
    d_tartaric = "O[C@@H](C(=O)O)[C@@H](O)C(=O)O"
    meso = "O[C@H](C(=O)O)[C@@H](O)C(=O)O"
    assert kit.compare(l_tartaric, d_tartaric)["relation"] == "enantiomers"
    assert kit.compare(l_tartaric, meso)["relation"] == "diastereomers"
    assert kit.compare(meso, "O[C@@H](C(=O)O)[C@H](O)C(=O)O")["relation"] == "identical"
    assert kit.compare("C/C=C/C(=O)O", "C/C=C\\C(=O)O")["relation"] == "diastereomers"
    assert kit.compare("CC(=O)CC(C)=O", "CC(O)=CC(C)=O")["relation"] == "tautomers"
    assert kit.compare("CCO", "COC")["relation"] == "constitutional_isomers"
    assert kit.compare("C[C@H](O)CC", "CC(O)CC")["relation"] == "same_constitution_stereo_undefined"
    different = kit.compare("CCO", "CCN")
    assert different["relation"] == "different_formula"
    assert different["formula_delta_b_minus_a"] == "+1H +1N -1O"
    charged = kit.compare("CCO", "CC[O-]")
    assert charged["relation"] == "different_formula" and charged["charge_delta_b_minus_a"] == -1


def test_unparseable_input_raises_a_structure_error():
    with pytest.raises(kit.StructureError):
        kit.compare("C1CC", "CCO")


def test_stereo_audit_separates_open_centres_from_ring_fixed_ones():
    undefined = kit.stereo_audit("CC(O)CC", enumerate_isomers=False)
    assert undefined["n_undefined"] == 1 and not undefined["fully_defined"]
    assert undefined["n_isomers_consistent_with_drawing"] == 2
    # bridgeheads of norbornane and cage carbons are not open stereocentres
    for cage in ("C1CC2CCC1C2", "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3"):
        audit = kit.stereo_audit(cage, enumerate_isomers=False)
        assert audit["fully_defined"] and audit["n_undefined"] == 0, cage
    # one marked bridgehead of camphor fixes the other
    camphor = kit.stereo_audit("CC1(C)[C@@H]2CCC1(C)C(=O)C2", enumerate_isomers=False)
    assert camphor["fully_defined"]
    assert kit.compare("CC1(C)[C@@H]2CCC1(C)C(=O)C2", "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2")["relation"] == "identical"
    assert kit.stereo_audit("O[C@H](C(=O)O)[C@@H](O)C(=O)O")["meso"] is True
    assert kit.stereo_audit("C1CCC2CCCCC2C1", enumerate_isomers=False)["n_isomers_consistent_with_drawing"] == 2


@pytest.mark.parametrize(
    "smiles, total, pairs, achiral",
    [
        ("OC(C(=O)O)C(O)C(=O)O", 3, 1, 1),  # tartaric acid
        ("CC1CCC(C)CC1", 2, 0, 2),  # 1,4-dimethylcyclohexane
        ("CC1CCCC(C)C1", 3, 1, 1),  # 1,3-dimethylcyclohexane
        ("C1CC2CCC1C2", 1, 0, 1),  # norbornane: bridgeheads cannot be varied
        ("CC(O)C(C)Cl", 4, 2, 0),
    ],
)
def test_stereoisomer_counts(smiles, total, pairs, achiral):
    result = kit.enumerate_stereoisomers(smiles)
    assert (result["n_stereoisomers"], result["n_enantiomer_pairs"], result["n_achiral_isomers"]) == (
        total,
        pairs,
        achiral,
    )


@pytest.mark.parametrize(
    "smiles, substituent, n, constitutional, with_stereo",
    [
        ("CC(C)CC", "Cl", 1, 4, 6),
        ("CC1CCCCC1", "Cl", 1, 5, 12),
        ("CCCC", "Cl", 1, 2, 3),
        ("CCC", "Cl", 2, 4, 5),
        ("c1ccccc1", "Cl", 2, 3, 3),
        ("c1ccc2ccccc2c1", "Br", 2, 10, 10),
        ("C1CCCCC1", "Cl", 2, 4, 9),
        ("C12C3C4C1C5C2C3C45", "Cl", 2, 3, 3),  # cubane
        ("C1CC2CCC1C2", "Cl", 1, 3, 6),  # norbornane: 1-, 2-exo/endo (x2), 7-
    ],
)
def test_substitution_product_counts(smiles, substituent, n, constitutional, with_stereo):
    result = kit.substitution_products(smiles, substituent, n_substitutions=n)
    assert result["n_constitutional_isomers"] == constitutional
    assert result["n_including_stereoisomers"] == with_stereo


@pytest.mark.parametrize(
    "smiles, carbons, hydrogens",
    [
        ("CCO", 2, 3),
        ("Cc1ccccc1", 5, 4),
        ("c1ccc2ccccc2c1", 3, 2),
        ("CC(C)(C)C", 2, 1),
        ("CC(=O)OCC", 4, 3),
        ("CC1CCCCC1", 5, 8),  # ring CH2 protons cis/trans to methyl are diastereotopic
        ("OCC(O)CO", 2, 5),  # glycerol: diastereotopic CH2 protons
    ],
)
def test_symmetry_environments(smiles, carbons, hydrogens):
    result = kit.symmetry_environments(smiles)
    assert result["n_carbon_environments"] == carbons
    assert result["n_hydrogen_environments"] == hydrogens


def test_enantiotopic_and_diastereotopic_protons_are_distinguished():
    ethanol = kit.symmetry_environments("CCO")
    assert ethanol["n_hydrogen_environments_chiral_medium"] == 4
    assert ethanol["atoms_with_diastereotopic_h"] == []
    assert kit.symmetry_environments("C[C@H](O)CC")["atoms_with_diastereotopic_h"] == [3]
