from codex_organic_chem.service import chem_input_review, chem_structure_review_result
from codex_organic_chem.structure_review import submit_review_item


def test_input_review_attaches_interactive_editor_session():
    result = chem_input_review(smiles="CC(=O)Oc1ccccc1C(=O)O")

    assert result["status"] == "awaiting_user_confirmation"
    assert result["session_id"]
    assert result["review_token"]
    assert result["review_url"].startswith("http://127.0.0.1:")
    assert result["items"][0]["status"] == "pending"
    # The static preview is still present for transcript context.
    assert "svg" in result["preview"]
    assert "chem_structure_review_result" in result["next_action"]


def test_input_review_session_round_trip_confirms_structure():
    started = chem_input_review(smiles="CCO")
    item_id = started["items"][0]["id"]

    submit_review_item(
        started["session_id"],
        item_id,
        {"status": "confirmed", "reviewed_smiles": "CCO"},
        token=started["review_token"],
    )
    final = chem_structure_review_result(started["session_id"], review_token=started["review_token"])

    assert final["status"] == "completed"
    assert final["items"][0]["status"] == "confirmed"
    assert final["items"][0]["canonical_smiles"] == "CCO"


def test_input_review_reaction_splits_into_component_items():
    result = chem_input_review(reaction_smiles="CBr.[OH-]>>CO.[Br-]")

    ids = [item["id"] for item in result["items"]]
    assert "reactant-1" in ids and "reactant-2" in ids
    assert "product-1" in ids and "product-2" in ids


def test_input_review_interactive_false_keeps_static_payload_only():
    result = chem_input_review(smiles="CCO", interactive=False)

    assert result["status"] == "awaiting_user_confirmation"
    assert "session_id" not in result
    assert "review_url" not in result
