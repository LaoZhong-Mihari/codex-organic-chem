from __future__ import annotations

import json
import mimetypes
import os
import secrets
import threading
import time
import webbrowser
from dataclasses import dataclass, field
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlencode, urlparse

from .rdkit_tools import normalize_structure


TERMINAL_ITEM_STATUSES = {"confirmed", "modified", "skipped", "error"}
DEFAULT_REVIEW_TIMEOUT_S = 1800


@dataclass
class ReviewItem:
    id: str
    label: str | None
    input_smiles: str
    initial_record: dict[str, Any]
    initial_warnings: list[str] = field(default_factory=list)
    result: dict[str, Any] | None = None

    def public_dict(self) -> dict[str, Any]:
        if self.result:
            return dict(self.result)
        return {
            "id": self.id,
            "label": self.label,
            "status": "pending",
            "input_smiles": self.input_smiles,
            "smiles": self.input_smiles,
            "warnings": list(self.initial_warnings),
            "initial_record": self.initial_record,
        }


@dataclass
class ReviewSession:
    id: str
    token: str
    items: list[ReviewItem]
    created_at: float
    expires_at: float
    status: str = "pending"
    completed_at: float | None = None
    warnings: list[str] = field(default_factory=list)
    condition: threading.Condition = field(default_factory=threading.Condition)


_SESSIONS: dict[str, ReviewSession] = {}
_SESSIONS_LOCK = threading.RLock()
_SERVER: ThreadingHTTPServer | None = None
_SERVER_THREAD: threading.Thread | None = None
_SERVER_BASE_URL: str | None = None


class ReviewHTTPServer(ThreadingHTTPServer):
    daemon_threads = True


def _utc_iso(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, timezone.utc).isoformat().replace("+00:00", "Z")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _ketcher_dist_path() -> Path:
    configured = os.environ.get("CODEX_CHEM_KETCHER_DIST")
    if configured:
        return Path(configured).expanduser().resolve()
    return (_repo_root() / "integrations" / "ketcher" / "dist").resolve()


def _read_json_body(handler: BaseHTTPRequestHandler) -> dict[str, Any]:
    length = int(handler.headers.get("Content-Length") or 0)
    if length <= 0:
        return {}
    raw = handler.rfile.read(length)
    if not raw:
        return {}
    payload = json.loads(raw.decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("JSON body must be an object.")
    return payload


def _session_by_id(session_id: str) -> ReviewSession | None:
    with _SESSIONS_LOCK:
        session = _SESSIONS.get(session_id)
    if session:
        _expire_session_if_needed(session)
    return session


def _expire_session_if_needed(session: ReviewSession) -> None:
    if session.status != "pending":
        return
    if time.time() <= session.expires_at:
        return
    with session.condition:
        if session.status == "pending":
            session.status = "expired"
            session.completed_at = time.time()
            session.warnings.append("Review session expired before all structures were confirmed.")
            session.condition.notify_all()


def _all_items_terminal(session: ReviewSession) -> bool:
    return all(item.result and item.result.get("status") in TERMINAL_ITEM_STATUSES for item in session.items)


def _item_by_id(session: ReviewSession, item_id: str) -> ReviewItem | None:
    for item in session.items:
        if item.id == item_id:
            return item
    return None


def _input_identity(item: ReviewItem) -> str | None:
    return item.initial_record.get("isomeric_smiles") or item.initial_record.get("canonical_smiles")


def _result_from_record(
    item: ReviewItem,
    submitted_status: str,
    reviewed_smiles: str | None,
    submitted_molfile: str | None,
    submitted_warnings: list[str] | None = None,
) -> dict[str, Any]:
    warnings = list(item.initial_warnings)
    warnings.extend(str(warning) for warning in (submitted_warnings or []))
    if submitted_status == "skipped":
        return {
            "id": item.id,
            "label": item.label,
            "status": "skipped",
            "input_smiles": item.input_smiles,
            "reviewed_smiles": reviewed_smiles or item.input_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": None,
            "inchi_key": None,
            "warnings": warnings + ["Structure was skipped during human review."],
        }
    if submitted_status == "error":
        return {
            "id": item.id,
            "label": item.label,
            "status": "error",
            "input_smiles": item.input_smiles,
            "reviewed_smiles": reviewed_smiles or item.input_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": submitted_molfile,
            "inchi_key": None,
            "warnings": warnings + ["Ketcher reported an error for this structure."],
        }

    if reviewed_smiles:
        record = normalize_structure(smiles=reviewed_smiles, source="structure_review")
    elif submitted_molfile:
        record = normalize_structure(molfile=submitted_molfile, source="structure_review")
    else:
        return {
            "id": item.id,
            "label": item.label,
            "status": "error",
            "input_smiles": item.input_smiles,
            "reviewed_smiles": None,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": None,
            "inchi_key": None,
            "warnings": warnings + ["Review submission did not include SMILES or Molfile."],
        }

    record_dict = record.to_dict()
    warnings.extend(record.warnings)
    reviewed_identity = record_dict.get("isomeric_smiles") or record_dict.get("canonical_smiles")
    if not reviewed_identity:
        return {
            "id": item.id,
            "label": item.label,
            "status": "error",
            "input_smiles": item.input_smiles,
            "reviewed_smiles": reviewed_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": submitted_molfile,
            "inchi_key": None,
            "warnings": warnings or ["Reviewed structure could not be normalized by RDKit."],
        }

    final_status = "modified" if submitted_status == "modified" or _input_identity(item) != reviewed_identity else "confirmed"
    return {
        "id": item.id,
        "label": item.label,
        "status": final_status,
        "input_smiles": item.input_smiles,
        "reviewed_smiles": reviewed_smiles or reviewed_identity,
        "canonical_smiles": record_dict.get("canonical_smiles"),
        "isomeric_smiles": record_dict.get("isomeric_smiles"),
        "molfile": record_dict.get("molblock"),
        "inchi_key": record_dict.get("inchi_key"),
        "warnings": warnings,
        "metadata": record_dict.get("metadata", {}),
    }


def _coerce_items(items: list[dict[str, Any]]) -> tuple[list[ReviewItem], list[str]]:
    coerced: list[ReviewItem] = []
    warnings: list[str] = []
    if not isinstance(items, list) or not items:
        raise ValueError("items must be a non-empty list of review item objects.")
    seen_ids: set[str] = set()
    for index, item in enumerate(items, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"Review item {index} must be an object.")
        smiles = str(item.get("smiles") or "").strip()
        if not smiles:
            raise ValueError(f"Review item {index} is missing a non-empty smiles field.")
        item_id = str(item.get("id") or f"item-{index}").strip() or f"item-{index}"
        if item_id in seen_ids:
            raise ValueError(f"Duplicate review item id: {item_id}")
        seen_ids.add(item_id)
        label = item.get("label")
        label = str(label) if label is not None else None
        initial = normalize_structure(smiles=smiles, source="structure_review_input")
        initial_dict = initial.to_dict()
        if initial.warnings:
            warnings.extend(f"{item_id}: {warning}" for warning in initial.warnings)
        coerced.append(
            ReviewItem(
                id=item_id,
                label=label,
                input_smiles=smiles,
                initial_record=initial_dict,
                initial_warnings=list(initial.warnings),
            )
        )
    return coerced, warnings


def _frontend_base_url() -> str:
    if os.environ.get("CODEX_CHEM_KETCHER_URL"):
        return os.environ["CODEX_CHEM_KETCHER_URL"].rstrip("/")
    return _server_base_url()


def _review_url(session: ReviewSession) -> str:
    query = urlencode({"session": session.id, "token": session.token, "api": _server_base_url()})
    return f"{_frontend_base_url()}/?{query}"


def _server_base_url() -> str:
    ensure_review_server()
    assert _SERVER_BASE_URL is not None
    return _SERVER_BASE_URL


def _session_payload(session: ReviewSession, include_token: bool = False) -> dict[str, Any]:
    _expire_session_if_needed(session)
    payload = {
        "status": session.status,
        "session_id": session.id,
        "api_url": _server_base_url(),
        "review_url": _review_url(session),
        "created_at": _utc_iso(session.created_at),
        "expires_at": _utc_iso(session.expires_at),
        "completed_at": _utc_iso(session.completed_at) if session.completed_at else None,
        "confirmation_required": session.status == "pending",
        "items": [item.public_dict() for item in session.items],
        "results": [item.result for item in session.items if item.result],
        "warnings": list(session.warnings),
    }
    dist_path = _ketcher_dist_path()
    if not os.environ.get("CODEX_CHEM_KETCHER_URL") and not (dist_path / "index.html").exists():
        payload["editor"] = "builtin_fallback"
        payload["warnings"].append(
            "Serving the built-in review editor (no build step required). For the full Ketcher "
            "drawing UI, run `cd integrations/ketcher && npm run build` or set CODEX_CHEM_KETCHER_URL."
        )
    else:
        payload["editor"] = "ketcher"
    if include_token:
        payload["review_token"] = session.token
    return payload


def _should_open_browser() -> bool:
    """Auto-open is on by default, but never during tests or when disabled."""
    if os.environ.get("CODEX_CHEM_REVIEW_OPEN_BROWSER", "").lower() in {"0", "false", "no"}:
        return False
    if "PYTEST_CURRENT_TEST" in os.environ:
        return False
    return True


def _open_review_in_browser(review_url: str) -> bool:
    if not _should_open_browser():
        return False
    try:
        return bool(webbrowser.open(review_url, new=1))
    except Exception:
        return False


def create_review_session(items: list[dict[str, Any]], timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S) -> dict[str, Any]:
    ensure_review_server()
    timeout = max(float(timeout_s), 0.001)
    coerced_items, warnings = _coerce_items(items)
    session = ReviewSession(
        id=secrets.token_urlsafe(18),
        token=secrets.token_urlsafe(24),
        items=coerced_items,
        created_at=time.time(),
        expires_at=time.time() + timeout,
        warnings=warnings,
    )
    with _SESSIONS_LOCK:
        _SESSIONS[session.id] = session
    return _session_payload(session, include_token=True)


def submit_review_item(
    session_id: str,
    item_id: str,
    payload: dict[str, Any],
    token: str | None = None,
) -> dict[str, Any]:
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and token != session.token:
        raise PermissionError("Invalid review session token.")
    with session.condition:
        if session.status != "pending":
            raise RuntimeError(f"Review session is already {session.status}.")
        item = _item_by_id(session, item_id)
        if not item:
            raise KeyError(f"Unknown review item id: {item_id}")
        submitted_status = str(payload.get("status") or "confirmed")
        if submitted_status not in TERMINAL_ITEM_STATUSES:
            raise ValueError(f"Unsupported item status: {submitted_status}")
        submitted_warnings = payload.get("warnings")
        if submitted_warnings is not None and not isinstance(submitted_warnings, list):
            submitted_warnings = [str(submitted_warnings)]
        item.result = _result_from_record(
            item,
            submitted_status=submitted_status,
            reviewed_smiles=payload.get("reviewed_smiles") or payload.get("smiles"),
            submitted_molfile=payload.get("molfile"),
            submitted_warnings=submitted_warnings,
        )
        if _all_items_terminal(session):
            session.status = "completed"
            session.completed_at = time.time()
        session.condition.notify_all()
        return _session_payload(session, include_token=False)


def complete_review_session(session_id: str, token: str | None = None) -> dict[str, Any]:
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and token != session.token:
        raise PermissionError("Invalid review session token.")
    with session.condition:
        if session.status == "pending":
            for item in session.items:
                if not item.result:
                    item.result = _result_from_record(item, "skipped", item.input_smiles, None)
            session.status = "completed"
            session.completed_at = time.time()
        session.condition.notify_all()
        return _session_payload(session, include_token=False)


def cancel_review_session(session_id: str, token: str | None = None) -> dict[str, Any]:
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and token != session.token:
        raise PermissionError("Invalid review session token.")
    with session.condition:
        if session.status == "pending":
            session.status = "cancelled"
            session.completed_at = time.time()
            session.warnings.append("Review session was cancelled by the user.")
        session.condition.notify_all()
        return _session_payload(session, include_token=False)


def review_session_result(
    session_id: str,
    token: str | None = None,
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
) -> dict[str, Any]:
    if wait:
        return wait_for_review_session(session_id, token=token, timeout_s=timeout_s)
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and token != session.token:
        raise PermissionError("Invalid review session token.")
    return _session_payload(session, include_token=False)


def wait_for_review_session(
    session_id: str,
    token: str | None = None,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
) -> dict[str, Any]:
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and token != session.token:
        raise PermissionError("Invalid review session token.")
    deadline = time.time() + max(float(timeout_s), 0.001)
    with session.condition:
        while session.status == "pending":
            _expire_session_if_needed(session)
            if session.status != "pending":
                break
            remaining = min(deadline, session.expires_at) - time.time()
            if remaining <= 0:
                session.status = "expired"
                session.completed_at = time.time()
                session.warnings.append("Review wait timed out before all structures were confirmed.")
                session.condition.notify_all()
                break
            session.condition.wait(timeout=remaining)
        return _session_payload(session, include_token=False)


def start_structure_review_batch(
    items: list[dict[str, Any]],
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
) -> dict[str, Any]:
    session = create_review_session(items, timeout_s=timeout_s)
    session["browser_opened"] = _open_review_in_browser(session["review_url"])
    if not session["browser_opened"]:
        session["warnings"] = [
            *session.get("warnings", []),
            f"Open the review page manually: {session['review_url']}",
        ]
    if wait:
        result = wait_for_review_session(
            session["session_id"], token=session["review_token"], timeout_s=timeout_s
        )
        result["browser_opened"] = session["browser_opened"]
        return result
    return session


def ensure_review_server() -> str:
    global _SERVER, _SERVER_BASE_URL, _SERVER_THREAD
    if _SERVER and _SERVER_BASE_URL:
        return _SERVER_BASE_URL
    with _SESSIONS_LOCK:
        if _SERVER and _SERVER_BASE_URL:
            return _SERVER_BASE_URL
        _SERVER = ReviewHTTPServer(("127.0.0.1", 0), ReviewRequestHandler)
        host, port = _SERVER.server_address[:2]
        _SERVER_BASE_URL = f"http://{host}:{port}"
        _SERVER_THREAD = threading.Thread(target=_SERVER.serve_forever, name="codex-chem-review-server", daemon=True)
        _SERVER_THREAD.start()
        return _SERVER_BASE_URL


# Dependency-free review UI served when the Ketcher build is absent. Structures
# are rendered server-side by RDKit (the /api/preview endpoint), so the page
# needs no JavaScript chemistry libraries: the user sees each depiction, edits
# the SMILES if needed with a live re-render, and confirms/skips each item.
_BUILTIN_REVIEW_PAGE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Structure Review</title>
<style>
 body{font:15px/1.45 Arial,Helvetica,sans-serif;margin:0;background:#f5f6f7;color:#111}
 header{background:#fff;border-bottom:1px solid #ddd;padding:14px 22px}
 h1{font-size:19px;margin:0}
 .sub{color:#666;font-size:13px;margin-top:3px}
 main{max-width:880px;margin:18px auto;padding:0 16px}
 .card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:16px;margin-bottom:14px;display:flex;gap:18px}
 .card.done{opacity:.62}
 .depiction{flex:0 0 240px;min-height:170px;display:flex;align-items:center;justify-content:center;border:1px solid #eee;border-radius:6px;background:#fff}
 .depiction svg{max-width:100%;max-height:220px;height:auto}
 .fields{flex:1;min-width:0}
 .label{font-weight:700;margin-bottom:2px}
 .warn{color:#9a3412;font-size:12.5px;margin:6px 0;white-space:pre-wrap}
 textarea{width:100%;box-sizing:border-box;font:13px Menlo,monospace;padding:7px;border:1px solid #ccc;border-radius:5px;min-height:44px}
 .row{margin-top:9px;display:flex;gap:8px;flex-wrap:wrap;align-items:center}
 button{font:600 13.5px Arial;padding:7px 14px;border-radius:6px;border:1px solid #bbb;background:#fff;cursor:pointer}
 button.primary{background:#0b57d0;border-color:#0b57d0;color:#fff}
 button:disabled{opacity:.45;cursor:default}
 .status{font-size:13px;font-weight:700}
 .status.confirmed{color:#0a7d33}.status.modified{color:#0b57d0}.status.skipped{color:#888}.status.invalid{color:#b00020}
 #banner{margin:14px 0;padding:11px 14px;border-radius:6px;background:#e7f0fe;border:1px solid #b9d2f8;display:none}
</style>
</head>
<body>
<header>
 <h1>Structure review</h1>
 <div class="sub">Check each rendered structure against the source. Edit the SMILES to correct it (preview updates), then confirm.</div>
</header>
<main><div id="banner"></div><div id="items"></div></main>
<script>
const params = new URLSearchParams(location.search);
const sessionId = params.get('session'), token = params.get('token');
const api = (params.get('api') || location.origin).replace(/\\/$/, '');
const hdrs = {'Content-Type':'application/json','X-Review-Token':token};
let sessionData = null;

async function fetchSession(){
  const r = await fetch(`${api}/api/review-sessions/${sessionId}?token=${encodeURIComponent(token)}`);
  if(!r.ok){banner(`Could not load review session (HTTP ${r.status}).`);throw new Error('load');}
  sessionData = await r.json();
}
function banner(text){const b=document.getElementById('banner');b.textContent=text;b.style.display='block';}
async function preview(smiles){
  const r = await fetch(`${api}/api/preview`,{method:'POST',headers:hdrs,body:JSON.stringify({smiles})});
  return r.json();
}
async function submit(itemId, body){
  const r = await fetch(`${api}/api/review-sessions/${sessionId}/items/${encodeURIComponent(itemId)}?token=${encodeURIComponent(token)}`,
    {method:'POST',headers:hdrs,body:JSON.stringify(body)});
  return r.json();
}
function render(){
  const host=document.getElementById('items');host.innerHTML='';
  for(const item of sessionData.items){
    const done = item.status && item.status !== 'pending';
    const card=document.createElement('div');card.className='card'+(done?' done':'');
    card.innerHTML=`
      <div class="depiction" id="dep-${item.id}">(rendering...)</div>
      <div class="fields">
        <div class="label">${item.label || item.id}</div>
        <div class="status ${item.status||'pending'}" id="st-${item.id}">${done?item.status:'awaiting review'}</div>
        <div class="warn" id="warn-${item.id}">${(item.warnings||[]).join('\\n')}</div>
        <textarea id="smi-${item.id}" ${done?'disabled':''}>${item.smiles || item.input_smiles || ''}</textarea>
        <div class="row">
          <button class="primary" id="ok-${item.id}" ${done?'disabled':''}>Confirm</button>
          <button id="re-${item.id}" ${done?'disabled':''}>Re-render preview</button>
          <button id="skip-${item.id}" ${done?'disabled':''}>Skip</button>
        </div>
      </div>`;
    host.appendChild(card);
    refreshPreview(item.id, item.smiles || item.input_smiles || '');
    if(!done){
      document.getElementById(`re-${item.id}`).onclick=()=>refreshPreview(item.id, val(item.id));
      document.getElementById(`ok-${item.id}`).onclick=()=>finish(item, 'confirmed');
      document.getElementById(`skip-${item.id}`).onclick=()=>finish(item, 'skipped');
      document.getElementById(`smi-${item.id}`).addEventListener('change',()=>refreshPreview(item.id, val(item.id)));
    }
  }
}
const val = id => document.getElementById(`smi-${id}`).value.trim();
async function refreshPreview(id, smiles){
  const dep=document.getElementById(`dep-${id}`);
  if(!smiles){dep.textContent='(empty SMILES)';return;}
  try{
    const p=await preview(smiles);
    if(p.svg){dep.innerHTML=p.svg;}
    else{dep.textContent='(invalid SMILES)';}
    const st=document.getElementById(`st-${id}`);
    if(st && !st.classList.contains('confirmed') && !st.classList.contains('skipped')){
      st.textContent = p.status==='ok' ? 'awaiting review' : 'invalid SMILES';
      st.className = 'status ' + (p.status==='ok' ? 'pending' : 'invalid');
    }
  }catch(e){dep.textContent='(preview failed)';}
}
async function finish(item, mode){
  const smiles = val(item.id);
  const body = mode==='skipped' ? {status:'skipped'} :
    {status: smiles===(item.input_smiles||'') ? 'confirmed' : 'modified', reviewed_smiles: smiles};
  const result = await submit(item.id, body);
  sessionData = result;
  render();
  if(result.status && result.status !== 'pending'){
    banner(`Review ${result.status}. You can close this tab; results are already available to the assistant.`);
  }
}
fetchSession().then(render).catch(()=>{});
</script>
</body>
</html>
"""


class ReviewRequestHandler(BaseHTTPRequestHandler):
    server_version = "CodexChemReview/0.1"

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        return

    def end_headers(self) -> None:
        origin = self.headers.get("Origin")
        self.send_header("Access-Control-Allow-Origin", origin or "*")
        self.send_header("Vary", "Origin")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, X-Review-Token")
        super().end_headers()

    def do_OPTIONS(self) -> None:
        self.send_response(204)
        self.end_headers()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        parts = [unquote(part) for part in parsed.path.strip("/").split("/") if part]
        if len(parts) == 3 and parts[:2] == ["api", "review-sessions"]:
            self._handle_get_session(parts[2], parsed.query)
            return
        self._serve_ketcher_static(parsed.path)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        parts = [unquote(part) for part in parsed.path.strip("/").split("/") if part]
        try:
            if parts == ["api", "preview"]:
                payload = _read_json_body(self)
                record = normalize_structure(smiles=str(payload.get("smiles") or ""), source="review_preview")
                record_dict = record.to_dict()
                self._write_json(
                    {
                        "status": "ok" if record_dict.get("canonical_smiles") else "invalid",
                        "svg": record_dict.get("svg"),
                        "canonical_smiles": record_dict.get("canonical_smiles"),
                        "warnings": record_dict.get("warnings", []),
                    }
                )
                return
            if len(parts) == 5 and parts[:2] == ["api", "review-sessions"] and parts[3] == "items":
                payload = _read_json_body(self)
                self._write_json(submit_review_item(parts[2], parts[4], payload, token=self._token(parsed.query)))
                return
            if len(parts) == 4 and parts[:2] == ["api", "review-sessions"] and parts[3] == "complete":
                self._write_json(complete_review_session(parts[2], token=self._token(parsed.query)))
                return
            if len(parts) == 4 and parts[:2] == ["api", "review-sessions"] and parts[3] == "cancel":
                self._write_json(cancel_review_session(parts[2], token=self._token(parsed.query)))
                return
            self._write_json({"status": "error", "warnings": ["Unknown review API endpoint."]}, status=404)
        except PermissionError as exc:
            self._write_json({"status": "error", "warnings": [str(exc)]}, status=403)
        except KeyError as exc:
            self._write_json({"status": "error", "warnings": [str(exc)]}, status=404)
        except Exception as exc:
            self._write_json({"status": "error", "warnings": [str(exc)]}, status=400)

    def _token(self, query: str) -> str | None:
        return self.headers.get("X-Review-Token") or parse_qs(query).get("token", [None])[0]

    def _authorize(self, session: ReviewSession, query: str) -> bool:
        return self._token(query) == session.token

    def _handle_get_session(self, session_id: str, query: str) -> None:
        session = _session_by_id(session_id)
        if not session:
            self._write_json({"status": "error", "warnings": ["Unknown review session."]}, status=404)
            return
        if not self._authorize(session, query):
            self._write_json({"status": "error", "warnings": ["Invalid review session token."]}, status=403)
            return
        self._write_json(_session_payload(session, include_token=False))

    def _serve_ketcher_static(self, request_path: str) -> None:
        dist = _ketcher_dist_path()
        if request_path in {"", "/"}:
            candidate = dist / "index.html"
        else:
            relative = request_path.lstrip("/")
            candidate = (dist / relative).resolve()
            if dist.exists():
                try:
                    candidate.relative_to(dist.resolve())
                except ValueError:
                    self.send_error(403)
                    return
            if not candidate.exists() and "." not in Path(relative).name:
                candidate = dist / "index.html"
        if candidate.exists() and candidate.is_file():
            content_type = mimetypes.guess_type(str(candidate))[0] or "application/octet-stream"
            data = candidate.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self._safe_write(data)
            return
        # No Ketcher build available: serve the dependency-free built-in
        # editor so human review still works out of the box.
        self._write_html(_BUILTIN_REVIEW_PAGE)

    def _write_json(self, payload: dict[str, Any], status: int = 200) -> None:
        data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self._safe_write(data)

    def _write_html(self, html: str, status: int = 200) -> None:
        data = html.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self._safe_write(data)

    def _safe_write(self, data: bytes) -> None:
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionAbortedError, ConnectionResetError):
            return
