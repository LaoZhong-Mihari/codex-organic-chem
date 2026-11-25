from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from codex_organic_chem import service
from codex_organic_chem import structure_review as sr
from codex_organic_chem.scheme_ocsr import load_gold_smiles


class _MockSocket:
    def __init__(self, raw_request: bytes) -> None:
        self._rfile = io.BytesIO(raw_request)
        self._wfile = io.BytesIO()

    def makefile(self, mode: str, *args, **kwargs):
        return self._rfile if "r" in mode else self._wfile

    def sendall(self, data: bytes) -> None:
        self._wfile.write(data)

    def close(self) -> None:
        pass


def _raw_http(method: str, path: str, body: dict | None = None, headers: dict[str, str] | None = None) -> tuple[int, dict[str, str], bytes]:
    hdrs = dict(headers or {})
    raw_body = b""
    if body is not None:
        raw_body = json.dumps(body).encode("utf-8")
        hdrs.setdefault("Content-Type", "application/json")
        hdrs.setdefault("Content-Length", str(len(raw_body)))
    hdrs.setdefault("Host", "127.0.0.1:8000")
    lines = [f"{method} {path} HTTP/1.1"] + [f"{k}: {v}" for k, v in hdrs.items()] + ["", ""]
    sock = _MockSocket("\r\n".join(lines).encode("utf-8") + raw_body)
    sr.ReviewRequestHandler(sock, ("127.0.0.1", 55555), None)
    raw_resp = sock._wfile.getvalue()
    head, _, payload = raw_resp.partition(b"\r\n\r\n")
    head_lines = head.decode("iso-8859-1").split("\r\n")
    status = int(head_lines[0].split()[1])
    resp_headers: dict[str, str] = {}
    for line in head_lines[1:]:
        if ":" in line:
            k, v = line.split(":", 1)
            resp_headers[k.strip().lower()] = v.strip()
    return status, resp_headers, payload


def test_mutating_post_endpoints_reject_missing_and_wrong_tokens(monkeypatch):
    monkeypatch.setattr(sr, "ensure_review_server", lambda *_a, **_kw: "http://127.0.0.1:8000")
    monkeypatch.setattr(sr, "_SERVER_BASE_URL", "http://127.0.0.1:8000")
    session = sr.create_review_session(items=[{"id": "1", "smiles": "CCO"}])
    sid = session["session_id"]
    token = session["review_token"]

    for endpoint, body in [
        (f"/api/review-sessions/{sid}/items/1", {"status": "confirmed", "reviewed_smiles": "CCO"}),
        (f"/api/review-sessions/{sid}/complete", {"verdict": "approved"}),
        (f"/api/review-sessions/{sid}/cancel", {}),
        (f"/api/review-sessions/{sid}/comment", {"comment": "hijacked"}),
        (f"/api/review-sessions/{sid}/route/step", {"key": "1", "decision": "approved"}),
    ]:
        status_missing, _, _ = _raw_http("POST", endpoint, body)
        assert status_missing == 403, f"{endpoint} allowed missing token"
        status_wrong, _, _ = _raw_http("POST", endpoint, body, {"X-Review-Token": "wrong-token"})
        assert status_wrong == 403, f"{endpoint} allowed wrong token"

    status_ok, _, raw = _raw_http(
        "POST",
        f"/api/review-sessions/{sid}/items/1",
        {"status": "confirmed", "reviewed_smiles": "CCO"},
        {"X-Review-Token": token},
    )
    assert status_ok == 200
    assert json.loads(raw.decode("utf-8"))["status"] == "completed"


def test_ketcher_static_blocks_path_traversal_when_dist_missing_or_present(tmp_path, monkeypatch):
    secret_file = tmp_path / "secret.txt"
    secret_file.write_text("TOP_SECRET_MARKER", encoding="utf-8")

    # Case 1: dist directory does not exist (default unbuilt state)
    missing_dist = tmp_path / "nested" / "ketcher" / "dist"
    monkeypatch.setenv("CODEX_CHEM_KETCHER_DIST", str(missing_dist))
    status, _, body = _raw_http("GET", "/../../../secret.txt")
    assert status == 404
    assert b"TOP_SECRET_MARKER" not in body

    # Case 2: dist directory exists; traversal via raw or URL-encoded .. is blocked
    missing_dist.mkdir(parents=True)
    (missing_dist / "index.html").write_text("<html>ok</html>", encoding="utf-8")
    for attack in ("/../../../secret.txt", "/%2e%2e/%2e%2e/%2e%2e/secret.txt"):
        status_attack, _, attack_body = _raw_http("GET", attack)
        assert status_attack in {403, 404}
        assert b"TOP_SECRET_MARKER" not in attack_body


def test_cors_blocks_external_origins_and_allows_loopback():
    status_ext, hdrs_ext, _ = _raw_http("GET", "/api/health", headers={"Origin": "https://evil.example"})
    assert status_ext == 200
    assert "access-control-allow-origin" not in hdrs_ext
    assert hdrs_ext.get("referrer-policy") == "no-referrer"

    opt_ext, _, _ = _raw_http("OPTIONS", "/api/health", headers={"Origin": "https://evil.example"})
    assert opt_ext == 403

    status_local, hdrs_local, _ = _raw_http("GET", "/api/health", headers={"Origin": "http://localhost:5173"})
    assert status_local == 200
    assert hdrs_local.get("access-control-allow-origin") == "http://localhost:5173"


def test_load_gold_smiles_parses_literals_without_executing_python_code(tmp_path):
    malicious = tmp_path / "evil_gold.py"
    marker = tmp_path / "pwned.txt"
    malicious.write_text(
        f"import pathlib\npathlib.Path({str(marker)!r}).write_text('owned')\nA_TO_I = [('A', '1', 'CCO')]\n",
        encoding="utf-8",
    )
    rows = load_gold_smiles(malicious)
    assert not marker.exists(), "load_gold_smiles executed arbitrary top-level Python code!"
    assert rows == [{"label": "A", "compound": "1", "smiles": "CCO"}]


def test_chem_python_scrubs_secrets_and_kills_child_process_groups(tmp_path, monkeypatch):
    monkeypatch.setenv("CODEX_CHEM_ENABLE_PYTHON", "1")
    monkeypatch.setenv("OPENAI_API_KEY", "sk-super-secret")
    monkeypatch.setenv("ANTHROPIC_AUTH_TOKEN", "ant-secret")

    res = service.chem_python("import os; print('OPENAI_API_KEY' in os.environ, 'ANTHROPIC_AUTH_TOKEN' in os.environ)")
    assert res["status"] == "ok"
    assert res["stdout"].strip() == "False False"

    pid_file = tmp_path / "child.pid"
    code = (
        "import subprocess, sys, time, pathlib\n"
        "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
        f"pathlib.Path({str(pid_file)!r}).write_text(str(p.pid))\n"
        "time.sleep(60)\n"
    )
    timed_out = service.chem_python(code, timeout_s=1)
    assert timed_out["status"] == "timeout"
    if pid_file.exists():
        import os

        child_pid = int(pid_file.read_text().strip())
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
