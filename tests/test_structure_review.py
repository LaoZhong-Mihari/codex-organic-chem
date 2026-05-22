import json
from urllib.request import Request, urlopen

from codex_organic_chem.service import chem_input_review, chem_structure_review_batch, chem_structure_review_result
from codex_organic_chem.structure_review import cancel_review_session, complete_review_session, submit_review_item


def _read_json(url: str) -> dict:
    with urlopen(url, timeout=5) as response:
        return json.loads(response.read().decode("utf-8"))


def _post_json(url: str, payload: dict) -> dict:
    request = Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urlopen(request, timeout=5) as response:
        return json.loads(response.read().decode("utf-8"))


def test_structure_review_batch_starts_local_session_api():
    result = chem_structure_review_batch(
        items=[{"id": "ethanol", "label": "Ethanol", "smiles": "CCO"}],
        wait=False,
        timeout_s=60,
    )

    assert result["status"] == "pending"
    assert result["confirmation_required"] is True
    assert result["review_url"].startswith("http://127.0.0.1:")
    assert result["api_url"].startswith("http://127.0.0.1:")
    assert result["items"][0]["id"] == "ethanol"
    assert result["items"][0]["status"] == "pending"

    fetched = _read_json(f"{result['api_url']}/api/review-sessions/{result['session_id']}?token={result['review_token']}")
    assert fetched["session_id"] == result["session_id"]
    assert fetched["items"][0]["input_smiles"] == "CCO"


def test_structure_review_http_item_submission_returns_normalized_result():
    started = chem_structure_review_batch(
        items=[{"id": "ethyl-ether", "smiles": "CCO"}],
        wait=False,
        timeout_s=60,
    )

    posted = _post_json(
        (
            f"{started['api_url']}/api/review-sessions/{started['session_id']}"
            f"/items/ethyl-ether?token={started['review_token']}"
        ),
        {"status": "modified", "reviewed_smiles": "CCOC"},
    )

    assert posted["status"] == "completed"
    assert posted["items"][0]["status"] == "modified"
    assert posted["items"][0]["canonical_smiles"] == "CCOC"


def test_structure_review_confirmed_item_is_normalized():
    started = chem_structure_review_batch(
        items=[{"id": "ethanol", "smiles": "CCO"}],
        wait=False,
        timeout_s=60,
    )

    result = submit_review_item(
        started["session_id"],
        "ethanol",
        {"status": "confirmed", "reviewed_smiles": "CCO"},
        token=started["review_token"],
    )

    item = result["items"][0]
    assert result["status"] == "completed"
    assert item["status"] == "confirmed"
    assert item["canonical_smiles"] == "CCO"
    assert item["isomeric_smiles"] == "CCO"
    assert item["molfile"]
    assert "metadata" in item


def test_structure_review_result_fetches_existing_session_state():
    started = chem_structure_review_batch(
        items=[{"id": "ethanol", "smiles": "CCO"}],
        wait=False,
        timeout_s=60,
    )
    submit_review_item(
        started["session_id"],
        "ethanol",
        {"status": "confirmed", "reviewed_smiles": "CCO"},
        token=started["review_token"],
    )

    result = chem_structure_review_result(
        session_id=started["session_id"],
        review_token=started["review_token"],
        wait=False,
    )

    assert result["status"] == "completed"
    assert result["items"][0]["status"] == "confirmed"
    assert result["items"][0]["canonical_smiles"] == "CCO"


def test_structure_review_modified_and_invalid_items_are_isolated():
    started = chem_structure_review_batch(
        items=[
            {"id": "first", "smiles": "CCO"},
            {"id": "second", "smiles": "CCN"},
        ],
        wait=False,
        timeout_s=60,
    )

    first = submit_review_item(
        started["session_id"],
        "first",
        {"status": "confirmed", "reviewed_smiles": "CCOC"},
        token=started["review_token"],
    )
    assert first["status"] == "pending"
    assert first["items"][0]["status"] == "modified"
    assert first["items"][0]["canonical_smiles"] == "CCOC"
    assert first["items"][1]["status"] == "pending"

    final = submit_review_item(
        started["session_id"],
        "second",
        {"status": "modified", "reviewed_smiles": "not-a-smiles"},
        token=started["review_token"],
    )

    assert final["status"] == "completed"
    assert final["items"][0]["status"] == "modified"
    assert final["items"][0]["canonical_smiles"] == "CCOC"
    assert final["items"][1]["status"] == "error"
    assert any("Invalid SMILES" in warning for warning in final["items"][1]["warnings"])


def test_structure_review_cancel_and_complete_pending_items():
    cancelled = chem_structure_review_batch(
        items=[{"id": "one", "smiles": "CCO"}],
        wait=False,
        timeout_s=60,
    )

    cancel_result = cancel_review_session(cancelled["session_id"], token=cancelled["review_token"])
    assert cancel_result["status"] == "cancelled"
    assert cancel_result["confirmation_required"] is False

    completed = chem_structure_review_batch(
        items=[{"id": "one", "smiles": "CCO"}, {"id": "two", "smiles": "CCN"}],
        wait=False,
        timeout_s=60,
    )
    complete_result = complete_review_session(completed["session_id"], token=completed["review_token"])
    assert complete_result["status"] == "completed"
    assert [item["status"] for item in complete_result["items"]] == ["skipped", "skipped"]


def test_structure_review_wait_timeout_expires_session():
    result = chem_structure_review_batch(
        items=[{"id": "slow", "smiles": "CCO"}],
        wait=True,
        timeout_s=0.01,
    )

    assert result["status"] == "expired"
    assert result["confirmation_required"] is False
    assert result["warnings"]


def test_legacy_input_review_still_returns_static_confirmation_payload():
    result = chem_input_review(smiles="CCO")

    assert result["status"] == "awaiting_user_confirmation"
    assert result["confirmation_required"] is True
    assert "svg" in result["preview"]
