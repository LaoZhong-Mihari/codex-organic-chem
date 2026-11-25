"""Review sessions: a local page where a person checks and corrects what the
assistant is about to rely on, and the HTTP server behind it.

Two kinds of session share one page (``review_ui/workbench.html``):

``structures``  a queue of molecules, each next to its source image when there
                is one. The reviewer confirms, corrects (click a stereocentre to
                flip it, pick an alternative, edit SMILES, draw in Ketcher),
                rejects or skips each one, and can leave a note.
``route``       a synthetic route with the ``route_check`` audit laid out step
                by step. Structures and conditions can be corrected in place;
                every edit re-runs the audit. Each step and each audit prompt
                can be answered, and the route is approved or sent back.

Every change is written to ``<CODEX_CHEM_HOME>/reviews/<session>.json``, so a
result survives the process that started the session and can be collected by
any later command.
"""

from __future__ import annotations

import base64
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

from . import review_chem
from .rdkit_tools import normalize_structure
from .review_route import RouteReview, strip_drawings


TERMINAL_ITEM_STATUSES = {"confirmed", "modified", "skipped", "rejected", "error"}
DEFAULT_REVIEW_TIMEOUT_S = 1800
REVIEW_MODES = ("auto", "ask", "off")
IMAGE_TYPES = {
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".gif": "image/gif",
    ".webp": "image/webp",
}
MAX_ASSET_BYTES = 25 * 1024 * 1024
MAX_JSON_BODY_BYTES = 32 * 1024 * 1024
LANGUAGES = ("en", "zh")


# --------------------------------------------------------------------------
# settings and storage


def home_dir() -> Path:
    configured = os.environ.get("CODEX_CHEM_HOME")
    return Path(configured).expanduser() if configured else Path.home() / ".codex-chem"


def reviews_dir() -> Path:
    path = home_dir() / "reviews"
    path.mkdir(parents=True, exist_ok=True)
    try:
        path.chmod(0o700)
    except OSError:
        pass
    return path


def _config_path() -> Path:
    return home_dir() / "config.json"


def _read_config() -> dict[str, Any]:
    try:
        data = json.loads(_config_path().read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def review_mode() -> str:
    """``auto`` (open the page at every review gate), ``ask`` (ask the user in
    chat first) or ``off``. ``CODEX_CHEM_REVIEW`` overrides the saved setting."""
    for candidate in (os.environ.get("CODEX_CHEM_REVIEW"), _read_config().get("review")):
        mode = str(candidate or "").strip().lower()
        if mode in REVIEW_MODES:
            return mode
    return "auto"


def set_review_mode(mode: str) -> dict[str, Any]:
    mode = str(mode).strip().lower()
    if mode not in REVIEW_MODES:
        raise ValueError(f"Review mode must be one of {', '.join(REVIEW_MODES)}.")
    config = _read_config()
    config["review"] = mode
    path = _config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return review_mode_status()


def review_mode_status() -> dict[str, Any]:
    env = str(os.environ.get("CODEX_CHEM_REVIEW") or "").strip().lower()
    return {
        "status": "ok",
        "review_mode": review_mode(),
        "saved_setting": _read_config().get("review"),
        "environment_override": env if env in REVIEW_MODES else None,
        "config_file": str(_config_path()),
        "meaning": {
            "auto": "open the review page at every gate without asking",
            "ask": "ask the user one line in chat before opening the page",
            "off": "never open the page unless the user asks for it",
        },
    }


def review_gate(confirmed: bool = False, what: str = "these structures") -> dict[str, Any] | None:
    """None when the page may be opened now, otherwise what to do instead.

    ``confirmed`` means the user asked for the page or agreed to it in chat.
    """
    mode = review_mode()
    if confirmed or mode == "auto":
        return None
    if mode == "off":
        return {
            "status": "review_disabled",
            "review_mode": "off",
            "confirmation_required": False,
            "next_action": (
                "Interactive review is switched off (CODEX_CHEM_REVIEW=off). Do not open a page. State your "
                "reading and its assumptions, label the result as not confirmed by a person, and continue. "
                "If the user asks for the page, call again with confirmed=true (CLI: --yes)."
            ),
        }
    return {
        "status": "ask_user_first",
        "review_mode": "ask",
        "confirmation_required": True,
        "user_prompt": (
            f"I can open a local page in your browser where you can check and correct {what} "
            "before I continue. Open it?"
        ),
        "next_action": (
            "Ask the user this one question. If they agree, call again with confirmed=true (CLI: --yes). "
            "If they decline, continue and label the result as not confirmed by a person."
        ),
    }


# --------------------------------------------------------------------------
# data model


@dataclass
class ReviewItem:
    id: str
    label: str | None
    input_smiles: str
    initial_record: dict[str, Any]
    initial_warnings: list[str] = field(default_factory=list)
    note: str | None = None
    question: str | None = None
    source: str | None = None
    image: str | None = None
    candidates: list[dict[str, str]] = field(default_factory=list)
    result: dict[str, Any] | None = None

    def context(self) -> dict[str, Any]:
        return {
            "note": self.note,
            "question": self.question,
            "source": self.source,
            "image": self.image,
            "candidates": list(self.candidates),
        }

    def public_dict(self) -> dict[str, Any]:
        if self.result:
            return {**self.context(), **self.result}
        return {
            "id": self.id,
            "label": self.label,
            "status": "pending",
            "input_smiles": self.input_smiles,
            "smiles": self.input_smiles,
            "warnings": list(self.initial_warnings),
            "initial_record": self.initial_record,
            **self.context(),
        }


@dataclass
class ReviewSession:
    id: str
    token: str
    items: list[ReviewItem]
    created_at: float
    expires_at: float
    kind: str = "structures"
    title: str | None = None
    intro: str | None = None
    language: str | None = None
    images: list[str] = field(default_factory=list)
    assets: dict[str, tuple[str, bytes]] = field(default_factory=dict)
    route: RouteReview | None = None
    comment: str = ""
    status: str = "pending"
    completed_at: float | None = None
    warnings: list[str] = field(default_factory=list)
    condition: threading.Condition = field(default_factory=threading.Condition)


_SESSIONS: dict[str, ReviewSession] = {}
_SESSIONS_LOCK = threading.RLock()
_SERVER: ThreadingHTTPServer | None = None
_SERVER_THREAD: threading.Thread | None = None
_SERVER_BASE_URL: str | None = None
_ADMIN_TOKEN: str | None = None  # set only in the detached review server
_LAST_ACTIVITY = time.time()


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


def _ketcher_available() -> bool:
    return bool(os.environ.get("CODEX_CHEM_KETCHER_URL")) or (_ketcher_dist_path() / "index.html").exists()


def _workbench_html() -> str:
    return (Path(__file__).resolve().parent / "review_ui" / "workbench.html").read_text(encoding="utf-8")


def _read_json_body(handler: BaseHTTPRequestHandler) -> dict[str, Any]:
    raw_len = handler.headers.get("Content-Length")
    try:
        length = int(raw_len or 0)
    except ValueError as exc:
        raise ValueError("Invalid Content-Length header.") from exc
    if length <= 0:
        return {}
    if length > MAX_JSON_BODY_BYTES:
        raise ValueError(f"Request body exceeds {MAX_JSON_BODY_BYTES // (1024 * 1024)} MB limit.")
    raw = handler.rfile.read(length)
    if not raw:
        return {}
    payload = json.loads(raw.decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("JSON body must be an object.")
    return payload


MAX_IN_MEMORY_TERMINAL_SESSIONS = 64
ASSET_RETENTION_AFTER_COMPLETE_S = 600


def _prune_sessions_locked() -> None:
    now = time.time()
    terminal: list[tuple[float, str]] = []
    for sid, sess in list(_SESSIONS.items()):
        if sess.status == "pending" and now > sess.expires_at:
            _expire_session_if_needed(sess)
        if sess.status != "pending":
            done_at = sess.completed_at or sess.expires_at or sess.created_at
            if sess.assets and now - done_at > ASSET_RETENTION_AFTER_COMPLETE_S:
                sess.assets.clear()
            terminal.append((done_at, sid))
    if len(terminal) > MAX_IN_MEMORY_TERMINAL_SESSIONS:
        terminal.sort()
        for _, sid in terminal[: len(terminal) - MAX_IN_MEMORY_TERMINAL_SESSIONS]:
            _SESSIONS.pop(sid, None)
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
            session.warnings.append("Review session expired before the review was finished.")
            session.condition.notify_all()
            _persist(session)


def pending_session_count() -> int:
    with _SESSIONS_LOCK:
        sessions = list(_SESSIONS.values())
    for session in sessions:
        _expire_session_if_needed(session)
    return sum(1 for session in sessions if session.status == "pending")


def last_activity() -> float:
    return _LAST_ACTIVITY


def _all_items_terminal(session: ReviewSession) -> bool:
    return all(item.result and item.result.get("status") in TERMINAL_ITEM_STATUSES for item in session.items)


def _item_by_id(session: ReviewSession, item_id: str) -> ReviewItem | None:
    for item in session.items:
        if item.id == item_id:
            return item
    return None


def _input_identity(item: ReviewItem) -> str | None:
    return item.initial_record.get("isomeric_smiles") or item.initial_record.get("canonical_smiles")


def _difference_from_input(item: ReviewItem, reviewed: str) -> dict[str, Any] | None:
    """What the reviewer changed, in words the assistant can act on."""
    try:
        described = review_chem.describe(reviewed, reference=item.input_smiles)
    except Exception:
        return None
    versus = described.get("versus_reference")
    if not versus or versus.get("relation") == "identical":
        return None
    return {k: v for k, v in versus.items() if k != "changed_atoms" and v not in (None, [], "")}


def _result_from_record(
    item: ReviewItem,
    submitted_status: str,
    reviewed_smiles: str | None,
    submitted_molfile: str | None,
    submitted_warnings: list[str] | None = None,
    comment: str | None = None,
) -> dict[str, Any]:
    warnings = list(item.initial_warnings)
    warnings.extend(str(warning) for warning in (submitted_warnings or []))
    base = {
        "id": item.id,
        "label": item.label,
        "input_smiles": item.input_smiles,
        "comment": (comment or "").strip(),
    }

    def without_structure(status: str, message: str, molfile: str | None = None) -> dict[str, Any]:
        return {
            **base,
            "status": status,
            "reviewed_smiles": reviewed_smiles or item.input_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": molfile,
            "inchi_key": None,
            "warnings": warnings + [message],
        }

    if submitted_status == "skipped":
        return without_structure("skipped", "Structure was skipped during human review.")
    if submitted_status == "rejected":
        return without_structure(
            "rejected",
            "The reviewer marked this structure as wrong without correcting it; do not use it. "
            "Read the comment, re-derive the structure and ask again.",
        )
    if submitted_status == "error":
        return without_structure("error", "Ketcher reported an error for this structure.", submitted_molfile)

    if reviewed_smiles:
        record = normalize_structure(smiles=reviewed_smiles, source="structure_review")
    elif submitted_molfile:
        record = normalize_structure(molfile=submitted_molfile, source="structure_review")
    else:
        result = without_structure("error", "Review submission did not include SMILES or Molfile.")
        result["reviewed_smiles"] = None
        return result

    record_dict = record.to_dict()
    warnings.extend(record.warnings)
    reviewed_identity = record_dict.get("isomeric_smiles") or record_dict.get("canonical_smiles")
    if not reviewed_identity:
        return {
            **base,
            "status": "error",
            "reviewed_smiles": reviewed_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "molfile": submitted_molfile,
            "inchi_key": None,
            "warnings": warnings or ["Reviewed structure could not be normalized by RDKit."],
        }

    final_status = "modified" if submitted_status == "modified" or _input_identity(item) != reviewed_identity else "confirmed"
    result = {
        **base,
        "status": final_status,
        "reviewed_smiles": reviewed_smiles or reviewed_identity,
        "canonical_smiles": record_dict.get("canonical_smiles"),
        "isomeric_smiles": record_dict.get("isomeric_smiles"),
        "molfile": record_dict.get("molblock"),
        "inchi_key": record_dict.get("inchi_key"),
        "warnings": warnings,
        "metadata": record_dict.get("metadata", {}),
    }
    if final_status == "modified":
        difference = _difference_from_input(item, reviewed_identity)
        if difference:
            result["changed_from_input"] = difference
    return result


# --------------------------------------------------------------------------
# building a session


def _add_asset(assets: dict[str, tuple[str, bytes]], path: str | os.PathLike[str]) -> str:
    source = Path(path).expanduser()
    mime = IMAGE_TYPES.get(source.suffix.lower())
    if mime is None:
        raise ValueError(f"Unsupported image type for review: {source.name} (use PNG, JPEG, GIF or WebP).")
    data = source.read_bytes()
    if len(data) > MAX_ASSET_BYTES:
        raise ValueError(f"Image is too large for review: {source.name}")
    name = f"image-{len(assets) + 1}{source.suffix.lower()}"
    assets[name] = (mime, data)
    return name


def _coerce_candidates(raw: Any, warnings: list[str], item_id: str) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for index, entry in enumerate(raw or [], start=1):
        smiles = entry.get("smiles") if isinstance(entry, dict) else entry
        label = (entry.get("label") if isinstance(entry, dict) else None) or f"Alternative {index}"
        try:
            candidates.append({"label": str(label), "smiles": review_chem.working_smiles(str(smiles or ""))})
        except Exception:
            warnings.append(f"{item_id}: alternative {index} could not be parsed and was dropped.")
    return candidates


def _coerce_items(
    items: list[dict[str, Any]], assets: dict[str, tuple[str, bytes]] | None = None
) -> tuple[list[ReviewItem], list[str]]:
    coerced: list[ReviewItem] = []
    warnings: list[str] = []
    assets = assets if assets is not None else {}
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
        image = item.get("image")
        if item.get("image_path"):
            try:
                image = _add_asset(assets, item["image_path"])
            except (OSError, ValueError) as exc:
                warnings.append(f"{item_id}: source image not attached ({exc}).")
                image = None
        elif image and image not in assets:
            warnings.append(f"{item_id}: source image {image!r} was not supplied.")
            image = None
        coerced.append(
            ReviewItem(
                id=item_id,
                label=label,
                input_smiles=smiles,
                initial_record=initial_dict,
                initial_warnings=list(initial.warnings),
                note=str(item["note"]) if item.get("note") else None,
                question=str(item["question"]) if item.get("question") else None,
                source=str(item["source"]) if item.get("source") else None,
                image=image,
                candidates=_coerce_candidates(item.get("candidates"), warnings, item_id),
            )
        )
    return coerced, warnings


def _frontend_base_url() -> str:
    return _server_base_url()


def _review_url(session: ReviewSession) -> str:
    query = {"session": session.id, "token": session.token, "api": _server_base_url()}
    if session.language in LANGUAGES:
        query["lang"] = session.language
    return f"{_frontend_base_url()}/?{urlencode(query)}"


def _server_base_url() -> str:
    ensure_review_server()
    assert _SERVER_BASE_URL is not None
    return _SERVER_BASE_URL


def _summary(session: ReviewSession) -> dict[str, Any]:
    if session.kind == "route" and session.route is not None:
        result = session.route.result()
        return {
            "verdict": result["verdict"],
            "needs_revision": result["needs_revision"],
            "edited_compounds": sorted(result["changed_compounds"]),
            "edited_conditions": sorted(result["changed_conditions"]),
            "flagged_steps": [d["step"] for d in result["step_decisions"] if d.get("decision") == "flagged"],
            "audit_status_after_review": result["audit_after_review"]["status"],
        }
    counts: dict[str, int] = {}
    for item in session.items:
        status = (item.result or {}).get("status", "pending")
        counts[status] = counts.get(status, 0) + 1
    return {
        "counts": counts,
        "all_usable": bool(session.items)
        and all((item.result or {}).get("status") in {"confirmed", "modified"} for item in session.items),
        "changed_items": [item.id for item in session.items if (item.result or {}).get("status") == "modified"],
        "unusable_items": [
            item.id for item in session.items if (item.result or {}).get("status") in {"rejected", "skipped", "error"}
        ],
    }


def _next_action(session: ReviewSession) -> str:
    if session.status == "pending":
        return (
            "The review page is open for the user. Tell them, then collect the result with "
            "chem_structure_review_result(session_id=..., wait=true) or "
            "`codex-chem review-result --session-id ... --wait`; a wait that times out leaves the session open."
        )
    if session.status != "completed":
        return f"The review ended as {session.status}; nothing was confirmed. Say so and ask how to proceed."
    if session.kind == "route":
        return (
            "Continue from route_result.route (the reviewer's version), never from your earlier copy. Address "
            "every flagged step, every audit prompt answered 'problem' and every comment, then re-run route-check."
        )
    return (
        "Use each item's canonical/isomeric SMILES from this result, never your pre-review SMILES. Items that "
        "are rejected, skipped or in error are not confirmed: read the comment and re-derive them."
    )


def _session_payload(session: ReviewSession, include_token: bool = False, for_page: bool = False) -> dict[str, Any]:
    _expire_session_if_needed(session)
    payload: dict[str, Any] = {
        "status": session.status,
        "kind": session.kind,
        "session_id": session.id,
        "title": session.title,
        "intro": session.intro,
        "language": session.language,
        "api_url": _server_base_url(),
        "review_url": _review_url(session),
        "created_at": _utc_iso(session.created_at),
        "expires_at": _utc_iso(session.expires_at),
        "completed_at": _utc_iso(session.completed_at) if session.completed_at else None,
        "confirmation_required": session.status == "pending",
        "items": [item.public_dict() for item in session.items],
        "results": [item.result for item in session.items if item.result],
        "comment": session.comment,
        "images": list(session.images),
        "summary": _summary(session),
        "warnings": list(session.warnings),
    }
    payload["editor"] = "ketcher" if _ketcher_available() else "builtin_fallback"
    if payload["editor"] == "builtin_fallback" and not for_page:
        payload["warnings"].append(
            "Free-hand drawing (Ketcher) is not built, so the page offers click-to-flip stereo edits, "
            "alternatives and SMILES/Molfile entry. For the drawing editor run "
            "`cd integrations/ketcher && npm install && npm run build`."
        )
    if session.kind == "route" and session.route is not None:
        payload["route_result"] = session.route.result()
        if for_page:
            payload["route"] = session.route.view()
    if for_page:
        payload["ketcher_url"] = os.environ.get("CODEX_CHEM_KETCHER_URL", "").rstrip("/") or None
    else:
        payload["next_action"] = _next_action(session)
        payload = strip_drawings(payload)
    if include_token:
        payload["review_token"] = session.token
    return payload


def _persist(session: ReviewSession) -> None:
    """Write the assistant-facing state to disk; never fail a request over it."""
    try:
        payload = _session_payload(session)
        path = reviews_dir() / f"{session.id}.json"
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
        try:
            temporary.chmod(0o600)
        except OSError:
            pass
        temporary.replace(path)
    except Exception:
        return


def _should_open_browser() -> bool:
    """Auto-open is on by default, but never during tests, CI, SSH, or headless Linux."""
    if os.environ.get("CODEX_CHEM_REVIEW_OPEN_BROWSER", "").lower() in {"0", "false", "no"}:
        return False
    if "PYTEST_CURRENT_TEST" in os.environ:
        return False
    if os.environ.get("CI", "").lower() in {"1", "true", "yes"} or os.environ.get("SSH_CONNECTION"):
        return False
    if os.name == "posix" and os.uname().sysname == "Linux":
        if not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY"):
            return False
    return True


def _open_review_in_browser(review_url: str) -> bool:
    if not _should_open_browser():
        return False
    try:
        return bool(webbrowser.open(review_url, new=1))
    except Exception:
        return False


def open_review_in_browser(review_url: str) -> bool:
    return _open_review_in_browser(review_url)


def create_review_session(
    items: list[dict[str, Any]] | None = None,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
    *,
    route: dict[str, Any] | None = None,
    title: str | None = None,
    intro: str | None = None,
    language: str | None = None,
    images: list[str] | None = None,
    assets: dict[str, tuple[str, bytes]] | None = None,
) -> dict[str, Any]:
    """Create a session. ``route`` makes it a route review; otherwise ``items``
    is the queue of structures. ``images`` are local paths of whole source
    pages shown above the queue; per-item crops go in ``image_path``."""
    ensure_review_server()
    timeout = max(float(timeout_s), 0.001)
    store: dict[str, tuple[str, bytes]] = dict(assets or {})
    warnings: list[str] = []
    session_images: list[str] = []
    for entry in images or []:
        if entry in store:
            session_images.append(entry)
            continue
        try:
            session_images.append(_add_asset(store, entry))
        except (OSError, ValueError) as exc:
            warnings.append(f"Source image not attached ({exc}).")
    route_review = None
    coerced_items: list[ReviewItem] = []
    if route is not None:
        route_review = RouteReview(route)
        title = title or route.get("title")
    else:
        coerced_items, item_warnings = _coerce_items(items or [], store)
        warnings.extend(item_warnings)
    session = ReviewSession(
        id=secrets.token_urlsafe(18),
        token=secrets.token_urlsafe(24),
        items=coerced_items,
        created_at=time.time(),
        expires_at=time.time() + timeout,
        kind="route" if route_review is not None else "structures",
        title=str(title) if title else None,
        intro=str(intro) if intro else None,
        language=language if language in LANGUAGES else None,
        images=session_images,
        assets=store,
        route=route_review,
        warnings=warnings,
    )
    with _SESSIONS_LOCK:
        _prune_sessions_locked()
        _SESSIONS[session.id] = session
    _persist(session)
    return _session_payload(session, include_token=True)


def create_review_session_from_request(payload: dict[str, Any]) -> dict[str, Any]:
    """Session creation over HTTP (detached server): images arrive as base64."""
    assets: dict[str, tuple[str, bytes]] = {}
    for name, entry in (payload.get("assets") or {}).items():
        mime = str(entry.get("mime") or "")
        if mime not in IMAGE_TYPES.values():
            raise ValueError(f"Unsupported image type: {mime}")
        data = base64.b64decode(entry.get("data_b64") or "")
        if len(data) > MAX_ASSET_BYTES:
            raise ValueError(f"Image is too large for review: {name}")
        assets[str(name)] = (mime, data)
    return create_review_session(
        items=payload.get("items"),
        timeout_s=payload.get("timeout_s") or DEFAULT_REVIEW_TIMEOUT_S,
        route=payload.get("route"),
        title=payload.get("title"),
        intro=payload.get("intro"),
        language=payload.get("language"),
        images=payload.get("images"),
        assets=assets,
    )


def _authorized_session(session_id: str, token: str | None) -> ReviewSession:
    session = _session_by_id(session_id)
    if not session:
        raise KeyError("Unknown review session.")
    if token is not None and not secrets.compare_digest(str(token), session.token):
        raise PermissionError("Invalid review session token.")
    return session


def _require_pending(session: ReviewSession) -> None:
    if session.status != "pending":
        raise RuntimeError(f"Review session is already {session.status}.")


def submit_review_item(
    session_id: str,
    item_id: str,
    payload: dict[str, Any],
    token: str | None = None,
) -> dict[str, Any]:
    session = _authorized_session(session_id, token)
    with session.condition:
        _require_pending(session)
        item = _item_by_id(session, item_id)
        if not item:
            raise KeyError(f"Unknown review item id: {item_id}")
        submitted_status = str(payload.get("status") or "confirmed")
        if submitted_status == "pending":
            item.result = None  # reopened by the reviewer
        else:
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
                comment=payload.get("comment"),
            )
        if "session_comment" in payload:
            session.comment = str(payload.get("session_comment") or "")
        if _all_items_terminal(session):
            session.status = "completed"
            session.completed_at = time.time()
        session.condition.notify_all()
        _persist(session)
        return _session_payload(session, include_token=False)


def set_session_comment(session_id: str, comment: str, token: str | None = None) -> dict[str, Any]:
    session = _authorized_session(session_id, token)
    with session.condition:
        _require_pending(session)
        session.comment = str(comment or "")
        if session.route is not None:
            session.route.comment = session.comment
        _persist(session)
        return _session_payload(session, include_token=False)


def update_route_review(
    session_id: str,
    action: str,
    payload: dict[str, Any],
    token: str | None = None,
) -> dict[str, Any]:
    """Apply one reviewer action to a route session and return the page state."""
    session = _authorized_session(session_id, token)
    with session.condition:
        _require_pending(session)
        if session.route is None:
            raise ValueError("This is not a route review session.")
        route = session.route
        if action == "compound":
            route.set_compound(str(payload.get("id")), smiles=payload.get("smiles"), molfile=payload.get("molfile"))
        elif action == "reset":
            route.reset_compound(str(payload["id"]) if payload.get("id") else None)
        elif action == "step":
            route.update_step(
                str(payload.get("key")),
                conditions=payload.get("conditions"),
                decision=payload.get("decision"),
                comment=payload.get("comment"),
            )
        elif action == "issue":
            route.answer_issue(str(payload.get("key")), payload.get("answer"), payload.get("comment"))
        else:
            raise ValueError(f"Unknown route review action: {action}")
        _persist(session)
        return _session_payload(session, for_page=True)


def complete_review_session(
    session_id: str,
    token: str | None = None,
    verdict: str | None = None,
    comment: str | None = None,
) -> dict[str, Any]:
    session = _authorized_session(session_id, token)
    with session.condition:
        if session.status == "pending":
            if comment is not None:
                session.comment = str(comment)
            if session.route is not None:
                session.route.finish(verdict, session.comment)
            for item in session.items:
                if not item.result:
                    item.result = _result_from_record(item, "skipped", item.input_smiles, None)
            session.status = "completed"
            session.completed_at = time.time()
        session.condition.notify_all()
        _persist(session)
        return _session_payload(session, include_token=False)


def cancel_review_session(session_id: str, token: str | None = None) -> dict[str, Any]:
    session = _authorized_session(session_id, token)
    with session.condition:
        if session.status == "pending":
            session.status = "cancelled"
            session.completed_at = time.time()
            session.warnings.append("Review session was cancelled by the user.")
        session.condition.notify_all()
        _persist(session)
        return _session_payload(session, include_token=False)


# --------------------------------------------------------------------------
# collecting a result


def _stored_result(session_id: str) -> dict[str, Any] | None:
    if not session_id or "/" in session_id or "\\" in session_id or session_id.startswith("."):
        return None
    path = reviews_dir() / f"{session_id}.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if payload.get("status") == "pending":
        try:
            expires = datetime.fromisoformat(str(payload.get("expires_at")).replace("Z", "+00:00")).timestamp()
        except ValueError:
            expires = None
        if expires is not None and time.time() > expires:
            payload["status"] = "expired"
            payload["confirmation_required"] = False
            payload.setdefault("warnings", []).append("Review session expired before the review was finished.")
            payload["next_action"] = "The review expired; nothing was confirmed. Say so and ask how to proceed."
    return payload


def _wait_for_stored_result(session_id: str, timeout_s: float) -> dict[str, Any] | None:
    deadline = time.time() + max(float(timeout_s), 0.0)
    payload = _stored_result(session_id)
    while payload is not None and payload.get("status") == "pending" and time.time() < deadline:
        time.sleep(0.4)
        payload = _stored_result(session_id) or payload
    if payload is not None and payload.get("status") == "pending":
        payload["wait_timed_out"] = True
    return payload


def review_session_result(
    session_id: str,
    token: str | None = None,
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
) -> dict[str, Any]:
    """State of a session held by this process, or by any other (from disk)."""
    if _session_by_id(session_id) is None:
        stored = _wait_for_stored_result(session_id, timeout_s) if wait else _stored_result(session_id)
        if stored is None:
            raise KeyError("Unknown review session.")
        return stored
    if wait:
        return wait_for_review_session(session_id, token=token, timeout_s=timeout_s)
    session = _authorized_session(session_id, token)
    return _session_payload(session, include_token=False)


def wait_for_review_session(
    session_id: str,
    token: str | None = None,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
) -> dict[str, Any]:
    """Block until the review ends or ``timeout_s`` passes.

    A wait that ends first leaves the session open (``wait_timed_out``); the
    session itself expires only at its own ``expires_at``.
    """
    session = _authorized_session(session_id, token)
    deadline = time.time() + max(float(timeout_s), 0.001)
    timed_out = False
    with session.condition:
        while session.status == "pending":
            _expire_session_if_needed(session)
            if session.status != "pending":
                break
            now = time.time()
            if now >= session.expires_at:
                session.status = "expired"
                session.completed_at = now
                session.warnings.append("Review wait timed out before the review was finished.")
                session.condition.notify_all()
                _persist(session)
                break
            if now >= deadline:
                timed_out = True
                break
            session.condition.wait(timeout=min(deadline, session.expires_at) - now)
        payload = _session_payload(session, include_token=False)
    if timed_out and payload["status"] == "pending":
        payload["wait_timed_out"] = True
    return payload


def _start(session: dict[str, Any], wait: bool, timeout_s: int | float) -> dict[str, Any]:
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


def start_structure_review_batch(
    items: list[dict[str, Any]],
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
    **options: Any,
) -> dict[str, Any]:
    return _start(create_review_session(items, timeout_s=timeout_s, **options), wait, timeout_s)


def start_route_review(
    route: dict[str, Any],
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
    **options: Any,
) -> dict[str, Any]:
    return _start(create_review_session(route=route, timeout_s=timeout_s, **options), wait, timeout_s)


def open_review(
    *,
    items: list[dict[str, Any]] | None = None,
    route: dict[str, Any] | None = None,
    wait: bool = False,
    timeout_s: int | float = DEFAULT_REVIEW_TIMEOUT_S,
    detached: bool = False,
    **options: Any,
) -> dict[str, Any]:
    """Open a review for the user: a queue of structures, or a route.

    ``detached`` hands the session to the background review server, so the page
    stays up after the calling process exits (what the CLI needs). Otherwise the
    session lives in this process (what a long-running MCP server needs).
    """
    options = {key: value for key, value in options.items() if value}
    if detached:
        from .review_daemon import start_detached

        payload = {"items": items, "route": route, "timeout_s": timeout_s, **options}
        session = start_detached({key: value for key, value in payload.items() if value is not None})
        if not wait:
            return session
        result = review_session_result(session["session_id"], wait=True, timeout_s=timeout_s)
        result["browser_opened"] = session["browser_opened"]
        return result
    if route is not None:
        return start_route_review(route, wait=wait, timeout_s=timeout_s, **options)
    return start_structure_review_batch(items or [], wait=wait, timeout_s=timeout_s, **options)


def ensure_review_server(port: int = 0, admin_token: str | None = None) -> str:
    global _SERVER, _SERVER_BASE_URL, _SERVER_THREAD, _ADMIN_TOKEN
    if _SERVER and _SERVER_BASE_URL:
        return _SERVER_BASE_URL
    with _SESSIONS_LOCK:
        if _SERVER and _SERVER_BASE_URL:
            return _SERVER_BASE_URL
        _SERVER = ReviewHTTPServer(("127.0.0.1", port), ReviewRequestHandler)
        if admin_token:
            _ADMIN_TOKEN = admin_token
        host, bound_port = _SERVER.server_address[:2]
        _SERVER_BASE_URL = f"http://{host}:{bound_port}"
        _SERVER_THREAD = threading.Thread(target=_SERVER.serve_forever, name="codex-chem-review-server", daemon=True)
        _SERVER_THREAD.start()
        return _SERVER_BASE_URL


# --------------------------------------------------------------------------
# HTTP


def _is_loopback_origin(origin: str | None) -> bool:
    if not origin or origin == "null":
        return False
    try:
        parsed = urlparse(origin)
    except ValueError:
        return False
    if parsed.scheme not in {"http", "https"}:
        return False
    return (parsed.hostname or "").lower() in {"127.0.0.1", "localhost", "::1"}


class ReviewRequestHandler(BaseHTTPRequestHandler):
    server_version = "CodexChemReview/2.0"

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        return

    def end_headers(self) -> None:
        origin = self.headers.get("Origin")
        if _is_loopback_origin(origin):
            self.send_header("Access-Control-Allow-Origin", str(origin))
            self.send_header("Vary", "Origin")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type, X-Review-Token, X-Admin-Token")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_OPTIONS(self) -> None:
        if self.headers.get("Origin") and not _is_loopback_origin(self.headers.get("Origin")):
            self.send_response(403)
            self.end_headers()
            return
        self.send_response(204)
        self.end_headers()

    def _touch(self) -> None:
        global _LAST_ACTIVITY
        _LAST_ACTIVITY = time.time()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        parts = [unquote(part) for part in parsed.path.strip("/").split("/") if part]
        if parts == ["api", "health"]:
            self._write_json({"status": "ok", "pending_sessions": pending_session_count()})
            return
        if len(parts) == 3 and parts[:2] == ["api", "review-sessions"]:
            self._handle_get_session(parts[2], parsed.query)
            return
        if len(parts) == 5 and parts[:2] == ["api", "review-sessions"] and parts[3] == "assets":
            self._handle_get_asset(parts[2], parts[4], parsed.query)
            return
        if parsed.path in {"", "/", "/index.html"}:
            self._touch()
            self._write_html(_workbench_html())
            return
        self._serve_ketcher_static(parsed.path)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        parts = [unquote(part) for part in parsed.path.strip("/").split("/") if part]
        try:
            payload = _read_json_body(self)
            if parts == ["api", "preview"]:
                self._touch()
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
            if parts == ["api", "review-sessions"]:
                admin_hdr = self.headers.get("X-Admin-Token") or ""
                if not _ADMIN_TOKEN or not secrets.compare_digest(admin_hdr, _ADMIN_TOKEN):
                    raise PermissionError("Session creation over HTTP is not enabled for this caller.")
                self._touch()
                self._write_json(create_review_session_from_request(payload))
                return
            if len(parts) >= 4 and parts[:2] == ["api", "review-sessions"]:
                result = self._session_action(parts[2], parts[3:], payload, self._token(parsed.query))
                self._touch()
                self._write_json(result)
                return
            self._write_json({"status": "error", "warnings": ["Unknown review API endpoint."]}, status=404)
        except PermissionError as exc:
            self._write_json({"status": "error", "warnings": [str(exc)]}, status=403)
        except KeyError as exc:
            self._write_json({"status": "error", "warnings": [str(exc.args[0]) if exc.args else str(exc)]}, status=404)
        except Exception as exc:
            self._write_json({"status": "error", "warnings": [str(exc)]}, status=400)

    def _session_action(
        self, session_id: str, rest: list[str], payload: dict[str, Any], token: str | None
    ) -> dict[str, Any]:
        req_token = token or ""
        session = _authorized_session(session_id, req_token)
        if rest[0] == "items" and len(rest) == 2:
            submit_review_item(session_id, rest[1], payload, token=req_token)
            with session.condition:
                return _session_payload(session, for_page=True)
        if rest == ["complete"]:
            complete_review_session(
                session_id, token=req_token, verdict=payload.get("verdict"), comment=payload.get("comment")
            )
            with session.condition:
                return _session_payload(session, for_page=True)
        if rest == ["cancel"]:
            cancel_review_session(session_id, token=req_token)
            with session.condition:
                return _session_payload(session, for_page=True)
        if rest == ["comment"]:
            set_session_comment(session_id, str(payload.get("comment") or ""), token=req_token)
            return {"status": "ok"}
        if rest[0] == "route" and len(rest) == 2:
            return update_route_review(session_id, rest[1], payload, token=req_token)
        # stateless chemistry for the page, still behind the session token
        del session
        options = {
            "reference": payload.get("reference"),
            "highlight_atoms": payload.get("highlight_atoms"),
            "show_indices": bool(payload.get("show_indices")),
            "width": int(payload.get("width") or 460),
            "height": int(payload.get("height") or 340),
        }
        if rest == ["describe"]:
            return review_chem.describe(payload.get("smiles"), payload.get("molfile"), **options)
        if rest == ["edit"]:
            try:
                edited = review_chem.apply_edit(
                    str(payload.get("smiles") or ""), str(payload.get("op")), payload.get("target"), payload.get("value")
                )
            except review_chem.StructureError as exc:
                return {"status": "invalid", "warnings": [str(exc)]}
            return review_chem.describe(edited, **options)
        raise KeyError("Unknown review API endpoint.")

    def _token(self, query: str) -> str | None:
        return self.headers.get("X-Review-Token") or parse_qs(query).get("token", [None])[0]

    def _authorize(self, session: ReviewSession, query: str) -> bool:
        supplied = self._token(query)
        return bool(supplied) and secrets.compare_digest(str(supplied), session.token)

    def _handle_get_session(self, session_id: str, query: str) -> None:
        session = _session_by_id(session_id)
        if not session:
            self._write_json({"status": "error", "warnings": ["Unknown review session."]}, status=404)
            return
        if not self._authorize(session, query):
            self._write_json({"status": "error", "warnings": ["Invalid review session token."]}, status=403)
            return
        self._touch()
        for_page = parse_qs(query).get("view", [""])[0] == "page"
        with session.condition:
            payload = _session_payload(session, include_token=False, for_page=for_page)
        self._write_json(payload)

    def _handle_get_asset(self, session_id: str, name: str, query: str) -> None:
        session = _session_by_id(session_id)
        if not session or not self._authorize(session, query) or name not in session.assets:
            self.send_error(404)
            return
        self._touch()
        mime, data = session.assets[name]
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self._safe_write(data)

    def _serve_ketcher_static(self, request_path: str) -> None:
        dist = _ketcher_dist_path().resolve()
        if not dist.is_dir():
            self.send_error(404)
            return
        clean_path = unquote(request_path)
        if clean_path.rstrip("/") == "/ketcher":
            candidate = (dist / "index.html").resolve()
        else:
            relative = clean_path.lstrip("/")
            candidate = (dist / relative).resolve()
        try:
            candidate.relative_to(dist)
        except ValueError:
            self.send_error(403)
            return
        if candidate.is_file():
            self._touch()
            content_type = mimetypes.guess_type(str(candidate))[0] or "application/octet-stream"
            data = candidate.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self._safe_write(data)
            return
        self.send_error(404)

    def _write_json(self, payload: dict[str, Any], status: int = 200) -> None:
        data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self._safe_write(data)

    def _write_html(self, html: str, status: int = 200) -> None:
        data = html.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header(
            "Content-Security-Policy",
            "default-src 'self' http://127.0.0.1:* http://localhost:*; "
            "img-src 'self' data: blob: http://127.0.0.1:* http://localhost:*; "
            "style-src 'self' 'unsafe-inline'; script-src 'self' 'unsafe-inline'; "
            "frame-src 'self' http://127.0.0.1:* http://localhost:*;",
        )
        self.end_headers()
        self._safe_write(data)

    def _safe_write(self, data: bytes) -> None:
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionAbortedError, ConnectionResetError):
            return
