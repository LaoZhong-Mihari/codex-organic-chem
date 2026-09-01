import json
import os
import subprocess
import sys
from urllib.parse import parse_qs, quote, urlparse
from urllib.request import Request, urlopen


def _review_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["CODEX_CHEM_REVIEW_OPEN_BROWSER"] = "0"
    return environment


def test_cli_input_review_without_wait_returns_static_preview_only():
    completed = subprocess.run(
        [sys.executable, "-m", "codex_organic_chem", "input-review", "--smiles", "CCO"],
        check=True,
        capture_output=True,
        text=True,
        env=_review_environment(),
    )

    payload = json.loads(completed.stdout)
    assert payload["status"] == "awaiting_user_confirmation"
    assert "session_id" not in payload
    assert "review_url" not in payload
    assert any("--wait" in warning for warning in payload["warnings"])


def test_cli_input_review_wait_keeps_server_alive_until_confirmation():
    process = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "codex_organic_chem",
            "input-review",
            "--smiles",
            "CCO",
            "--wait",
            "--timeout-s",
            "15",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=_review_environment(),
    )
    try:
        assert process.stderr is not None
        review_line = process.stderr.readline().strip()
        assert review_line.startswith("Review URL: ")
        assert process.poll() is None

        review_url = review_line.removeprefix("Review URL: ")
        parsed = urlparse(review_url)
        query = parse_qs(parsed.query)
        session_id = query["session"][0]
        review_token = query["token"][0]
        api_url = query["api"][0]

        with urlopen(review_url, timeout=3) as response:
            assert response.status == 200

        submit_url = (
            f"{api_url}/api/review-sessions/{quote(session_id, safe='')}/items/input"
            f"?token={quote(review_token, safe='')}"
        )
        request = Request(
            submit_url,
            data=json.dumps({"status": "confirmed", "reviewed_smiles": "CCO"}).encode(),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urlopen(request, timeout=3) as response:
            assert json.load(response)["status"] == "completed"

        stdout, trailing_stderr = process.communicate(timeout=10)
        assert process.returncode == 0, trailing_stderr
        result = json.loads(stdout)
        assert result["status"] == "completed"
        assert result["items"][0]["status"] == "confirmed"
        assert result["items"][0]["canonical_smiles"] == "CCO"
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=3)
