"""Mass and charge accounting checks for written reactions.

These reactions are all parseable and all pass functional-group classification;
what distinguishes them is whether the equation as written conserves atoms and
charge. Before the balance check existed, every case below reported no warning.
"""

from codex_organic_chem.reaction import ReactionParts, mass_charge_balance
from codex_organic_chem.service import chem_reaction_analyze


def _balance(reaction_smiles: str) -> dict:
    return chem_reaction_analyze(reaction_smiles)["analysis"]["mass_charge_balance"]


def test_hydrogen_only_gain_is_reported_as_missing_reductant():
    result = _balance("CC(=O)C>>CC(O)C")
    assert result["status"] == "unbalanced"
    assert result["hydrogen_delta"] == 2
    assert result["element_delta_without_reagents"] == {"H": 2}
    assert any("net reduction or protonation" in message for message in result["messages"])


def test_hydrogen_only_loss_is_reported_as_missing_oxidant():
    result = _balance("CC(O)C>>CC(=O)C")
    assert result["status"] == "unbalanced"
    assert result["hydrogen_delta"] == -2
    assert any("net oxidation or deprotonation" in message for message in result["messages"])


def test_balanced_substitution_has_no_balance_message():
    result = _balance("CCBr.[OH-]>>CCO.[Br-]")
    assert result["status"] == "balanced"
    assert result["messages"] == []
    assert result["reactant_charge"] == -1
    assert result["product_charge"] == -1


def test_omitted_water_in_esterification_is_flagged():
    result = _balance("CC(=O)O.CCO>>CC(=O)OCC")
    assert result["status"] == "unbalanced"
    assert result["element_delta_without_reagents"] == {"H": -2, "O": -1}
    assert any("Atom balance does not close" in message for message in result["messages"])


def test_esterification_with_water_written_balances():
    assert _balance("CC(=O)O.CCO>>CC(=O)OCC.O")["status"] == "balanced"


def test_lost_charge_is_reported_separately_from_atom_imbalance():
    result = _balance("CC(=O)C.[OH-]>>CC(O)C")
    assert result["charge_delta"] == 1
    assert any("Formal charge is not conserved" in message for message in result["messages"])
    assert any("Atom balance does not close" in message for message in result["messages"])


def test_catalytic_reagent_slot_does_not_create_a_false_charge_warning():
    """A catalyst in the reagent slot must not be counted against the balance.

    Charge has to follow the same include/exclude decision as the atom count;
    summing charge over reactants+reagents unconditionally made this proton
    look like it vanished.
    """
    result = _balance("CC(=O)O.CCO>[H+]>CC(=O)OCC.O")
    assert result["status"] == "balanced"
    assert result["reagents_counted_in_balance"] is False
    assert result["messages"] == []


def test_stoichiometric_reagent_slot_is_counted_when_it_closes_the_balance():
    result = _balance("CCBr>[OH-]>CCO.[Br-]")
    assert result["status"] == "balanced"
    assert result["reagents_counted_in_balance"] is True
    assert result["messages"] == []


def test_unbalanced_reaction_is_scored_below_a_merely_unmapped_one():
    unbalanced = chem_reaction_analyze("CC(=O)C>>CC(O)C")
    balanced = chem_reaction_analyze("CCBr.[OH-]>>CCO.[Br-]")
    assert unbalanced["confidence"] < balanced["confidence"]
    assert unbalanced["confidence"] <= 0.3


def test_balance_messages_surface_in_top_level_warnings():
    result = chem_reaction_analyze("CC(=O)C>>CC(O)C")
    assert any("hydride/proton source is missing" in warning for warning in result["warnings"])
    assert result["analysis"]["tool_facts"]["mass_charge_balance_status"] == "unbalanced"


def test_balance_reports_formulas_for_each_role():
    result = mass_charge_balance(ReactionParts(reactants=["CCBr"], reagents=["[OH-]"], products=["CCO", "[Br-]"]))
    assert result["reactant_formula"] == "C2H5Br"
    assert result["reagent_formula"] == "HO"
    assert result["product_formula"] == "C2H6BrO"
