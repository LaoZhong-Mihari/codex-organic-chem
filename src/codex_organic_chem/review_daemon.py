"""Detached review server, so a review page outlives the command that opened it.

A CLI command exits in seconds, but a person may take minutes to check a
structure. ``start_detached`` hands the session to a small background server
(started on demand, one per ``CODEX_CHEM_HOME``) and returns at once. The
server stops by itself after a stretch with no open session. Results are read
back from the session file on disk, so collecting them needs no connection.

    python -m codex_organic_chem.review_daemon      # what start_detached spawns
"""

from __future__ import annotations

import base64
import json
import os
import secrets
import subprocess
import sys
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from . import structure_review as sr

IDLE_SHUTDOWN_S = 900
STARTUP_TIMEOUT_S = 20


def _state_path() -> Path:
    return sr.home_dir() / "review-server.json"


def _code_version() -> str:
    """Changes whenever the review code or page changes, so an old server is
    replaced instead of serving a stale page."""
    here = Path(__file__).resolve().parent
    names = ["structure_review.py", "review_chem.py", "review_route.py", "review_daemon.py", "route_check.py",
             "structure_compare.py", "review_ui/workbench.html"]
    stamps = []
    for name in names:
        try:
            stamps.append(f"{name}:{int((here / name).stat().st_mtime)}")
        except OSError:
            stamps.append(f"{name}:missing")
    return "|".join(stamps)


def _read_state() -> dict[str, Any] | None:
    try:
        state = json.loads(_state_path().read_text(encoding="utf-8"))
        return state if isinstance(state, dict) else None
    except (OSError, ValueError):
        return None


def _request(url: str, payload: dict[str, Any] | None = None, headers: dict[str, str] | None = None, timeout: float = 30) -> dict[str, Any]:
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    request = Request(url, data=data, headers={"Content-Type": "application/json", **(headers or {})})
    try:
        with urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        try:
            body = json.loads(exc.read().decode("utf-8"))
        except Exception:
            body = {}
        raise RuntimeError("; ".join(body.get("warnings") or []) or f"Review server answered HTTP {exc.code}.") from exc


def _alive(state: dict[str, Any] | None) -> bool:
    if not state or state.get("code_version") != _code_version():
        return False
    try:
        return _request(f"{state['base_url']}/api/health", timeout=2).get("status") == "ok"
    except (RuntimeError, URLError, OSError, ValueError, KeyError):
        return False


def ensure_daemon() -> dict[str, Any]:
    """State of a running review server, starting one if there is none."""
    state = _read_state()
    if _alive(state):
        return state  # type: ignore[return-value]
    home = sr.home_dir()
    home.mkdir(parents=True, exist_ok=True)
    try:
        _state_path().unlink()
    except OSError:
        pass
    log = open(home / "review-server.log", "ab")  # noqa: SIM115 - handed to the child
    environment = os.environ.copy()
    environment["CODEX_CHEM_HOME"] = str(home)
    subprocess.Popen(
        [sys.executable, "-m", "codex_organic_chem.review_daemon"],
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=log,
        start_new_session=True,
        env=environment,
        cwd=str(home),
    )
    deadline = time.time() + STARTUP_TIMEOUT_S
    while time.time() < deadline:
        state = _read_state()
        if _alive(state):
            return state  # type: ignore[return-value]
        time.sleep(0.15)
    raise RuntimeError(f"The review server did not start; see {home / 'review-server.log'}.")


def _attach_images(payload: dict[str, Any]) -> dict[str, Any]:
    """Replace local image paths by inline data the server can hold."""
    payload = json.loads(json.dumps(payload))  # deep copy of plain data
    assets: dict[str, dict[str, str]] = {}
    warnings: list[str] = []

    def add(path: str) -> str | None:
        source = Path(path).expanduser()
        mime = sr.IMAGE_TYPES.get(source.suffix.lower())
        try:
            if mime is None:
                raise ValueError("use PNG, JPEG, GIF or WebP")
            data = source.read_bytes()
            if len(data) > sr.MAX_ASSET_BYTES:
                raise ValueError("file is too large")
        except (OSError, ValueError) as exc:
            warnings.append(f"Source image not attached: {source.name} ({exc}).")
            return None
        name = f"image-{len(assets) + 1}{source.suffix.lower()}"
        assets[name] = {"mime": mime, "data_b64": base64.b64encode(data).decode("ascii")}
        return name

    for item in payload.get("items") or []:
        if isinstance(item, dict) and item.get("image_path"):
            name = add(str(item.pop("image_path")))
            if name:
                item["image"] = name
    payload["images"] = [name for name in (add(str(p)) for p in payload.get("images") or []) if name]
    payload["assets"] = assets
    payload["_client_warnings"] = warnings
    return payload


def start_detached(payload: dict[str, Any], open_browser: bool = True) -> dict[str, Any]:
    """Create a session on the detached server and open it for the user."""
    state = ensure_daemon()
    prepared = _attach_images(payload)
    client_warnings = prepared.pop("_client_warnings")
    session = _request(
        f"{state['base_url']}/api/review-sessions", prepared, headers={"X-Admin-Token": state["admin_token"]}
    )
    session["warnings"] = [*session.get("warnings", []), *client_warnings]
    session["browser_opened"] = sr.open_review_in_browser(session["review_url"]) if open_browser else False
    if not session["browser_opened"]:
        session["warnings"].append(f"Open the review page manually: {session['review_url']}")
    session["server"] = "detached"
    return session


def stop_daemon() -> bool:
    state = _read_state()
    if not state:
        return False
    try:
        os.kill(int(state["pid"]), 15)
    except (OSError, KeyError, ValueError):
        return False
    return True


def main() -> None:
    admin_token = secrets.token_urlsafe(24)
    base_url = sr.ensure_review_server(port=int(os.environ.get("CODEX_CHEM_REVIEW_PORT") or 0), admin_token=admin_token)
    state_path = _state_path()
    state_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = state_path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(
            {
                "pid": os.getpid(),
                "base_url": base_url,
                "admin_token": admin_token,
                "code_version": _code_version(),
                "started_at": time.time(),
            }
        ),
        encoding="utf-8",
    )
    try:
        temporary.chmod(0o600)
    except OSError:
        pass
    temporary.replace(state_path)
    idle_s = float(os.environ.get("CODEX_CHEM_REVIEW_IDLE_S") or IDLE_SHUTDOWN_S)
    try:
        while True:
            time.sleep(1.0)
            current = _read_state()
            if not current or current.get("pid") != os.getpid():
                break  # a newer server has taken over the state file
            if sr.pending_session_count() == 0 and time.time() - sr.last_activity() > idle_s:
                break
    finally:
        current = _read_state()
        if current and current.get("pid") == os.getpid():
            try:
                state_path.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    main()
