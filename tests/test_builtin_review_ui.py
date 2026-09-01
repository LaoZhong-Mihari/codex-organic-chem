import json
from urllib.request import Request, urlopen

from codex_organic_chem.service import chem_structure_review_batch
from codex_organic_chem.structure_review import _server_base_url


def _post(url: str, payload: dict) -> dict:
    request = Request(url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    with urlopen(request) as response:
        return json.loads(response.read())


def test_builtin_editor_serves_when_ketcher_dist_missing(monkeypatch):
    monkeypatch.setenv("CODEX_CHEM_KETCHER_DIST", "/nonexistent/ketcher/dist")
    monkeypatch.delenv("CODEX_CHEM_KETCHER_URL", raising=False)

    result = chem_structure_review_batch(items=[{"id": "e", "smiles": "CCO"}])

    assert result["editor"] == "builtin_fallback"
    # Auto-open is suppressed under pytest, so the URL is surfaced instead.
    assert result["browser_opened"] is False
    assert any(result["review_url"] in warning for warning in result["warnings"])

    with urlopen(result["review_url"]) as response:
        page = response.read().decode()
    assert response.status == 200
    assert "Structure review" in page


def test_preview_endpoint_renders_and_rejects(monkeypatch):
    monkeypatch.setenv("CODEX_CHEM_KETCHER_DIST", "/nonexistent/ketcher/dist")
    chem_structure_review_batch(items=[{"id": "e", "smiles": "CCO"}])

    ok = _post(f"{_server_base_url()}/api/preview", {"smiles": "c1ccccc1"})
    assert ok["status"] == "ok"
    assert "<svg" in ok["svg"]

    bad = _post(f"{_server_base_url()}/api/preview", {"smiles": "not-a-smiles"})
    assert bad["status"] == "invalid"


def test_builtin_flow_modified_smiles_round_trip(monkeypatch):
    monkeypatch.setenv("CODEX_CHEM_KETCHER_DIST", "/nonexistent/ketcher/dist")
    session = chem_structure_review_batch(items=[{"id": "asp", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}])

    final = _post(
        f"{_server_base_url()}/api/review-sessions/{session['session_id']}/items/asp"
        f"?token={session['review_token']}",
        {"status": "modified", "reviewed_smiles": "CC(=O)Oc1ccccc1C(=O)OC"},
    )
    assert final["status"] == "completed"
    assert final["items"][0]["status"] == "modified"
    assert final["items"][0]["canonical_smiles"] == "COC(=O)c1ccccc1OC(C)=O"
