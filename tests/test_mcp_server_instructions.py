from codex_organic_chem import mcp_server


def test_server_instructions_point_to_the_skill(monkeypatch):
    monkeypatch.delenv("CODEX_CHEM_MCP_INSTRUCTIONS", raising=False)
    text = mcp_server.server_instructions()
    assert text is not None
    assert "organic-chemistry-assistant" in text
    assert "before your first call" in text


def test_server_instructions_can_be_switched_off(monkeypatch):
    for value in ("off", "0", "false", "OFF"):
        monkeypatch.setenv("CODEX_CHEM_MCP_INSTRUCTIONS", value)
        assert mcp_server.server_instructions() is None
    monkeypatch.setenv("CODEX_CHEM_MCP_INSTRUCTIONS", "on")
    assert mcp_server.server_instructions() == mcp_server.SERVER_INSTRUCTIONS
