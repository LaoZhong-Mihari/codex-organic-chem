import pytest


@pytest.fixture(autouse=True)
def _isolated_codex_chem_home(tmp_path_factory, monkeypatch):
    """Keep review sessions, the saved review mode and the detached review
    server of a test run out of the developer's own ~/.codex-chem."""
    monkeypatch.setenv("CODEX_CHEM_HOME", str(tmp_path_factory.mktemp("codex-chem-home")))
    monkeypatch.delenv("CODEX_CHEM_REVIEW", raising=False)
    monkeypatch.setenv("CODEX_CHEM_REVIEW_IDLE_S", "2")
    yield
