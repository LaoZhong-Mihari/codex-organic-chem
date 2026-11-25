from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from collections.abc import Sequence
from typing import Any

from .models import ToolStatus

TOOL_ENV_VARS = {
    "obabel": "CODEX_CHEM_OBABEL_PATH",
    "osra": "CODEX_CHEM_OSRA_PATH",
    "xtb": "CODEX_CHEM_XTB_PATH",
    "crest": "CODEX_CHEM_CREST_PATH",
}


EXTRA_SEARCH_DIRS = (Path("/opt/homebrew/bin"), Path("/usr/local/bin"))


def bundled_tool_prefix() -> Path:
    configured = os.environ.get("CODEX_CHEM_TOOL_PREFIX")
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".local" / "share" / "codex-organic-chem" / "conda-tools"


_STATUS_CACHE: dict[tuple[str, tuple[str, ...]], ToolStatus] = {}


def _candidate_executable_paths(name: str) -> list[str]:
    candidates: list[str] = []
    env_name = TOOL_ENV_VARS.get(name, f"CODEX_CHEM_{name.upper()}_PATH")
    env_value = os.environ.get(env_name)
    if env_value:
        expanded = Path(env_value).expanduser()
        try:
            if expanded.exists():
                candidates.append(str(expanded))
        except OSError:
            pass
    bundled = bundled_tool_prefix() / "bin" / name
    bundled_ocsr = Path.home() / ".local" / "share" / "codex-organic-chem" / "ocsr-tools" / "osra-osx64" / "bin" / name
    if name == "osra":
        try:
            if bundled_ocsr.exists():
                candidates.append(str(bundled_ocsr))
        except OSError:
            pass
    if name in {"xtb", "crest"}:
        try:
            if bundled.exists():
                candidates.append(str(bundled))
        except OSError:
            pass
    found = shutil.which(name)
    if found:
        candidates.append(found)
    for directory in EXTRA_SEARCH_DIRS:
        extra = directory / name
        try:
            if extra.exists():
                candidates.append(str(extra))
        except OSError:
            continue
    try:
        if bundled.exists():
            candidates.append(str(bundled))
    except OSError:
        pass
    deduped: list[str] = []
    seen: set[str] = set()
    for item in candidates:
        if item not in seen:
            seen.add(item)
            deduped.append(item)
    return deduped


def executable_path(name: str) -> str | None:
    candidates = _candidate_executable_paths(name)
    for candidate in candidates:
        try:
            if os.access(candidate, os.X_OK):
                return candidate
        except OSError:
            continue
    return candidates[0] if candidates else None


def json_dumps(data: Any) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True)


def extract_version_line(text: str) -> str | None:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return None
    patterns = [
        re.compile(r"\bversion\b", re.IGNORECASE),
        re.compile(r"\bv?\d+\.\d+(?:\.\d+)?\b"),
    ]
    for pattern in patterns:
        for line in lines:
            if pattern.search(line):
                return line
    return lines[0]


def _xtb_opt_healthy(candidate: str, timeout_s: int = 5) -> tuple[bool, str | None]:
    import tempfile

    with tempfile.TemporaryDirectory(prefix="codex-chem-xtb-probe-") as tmp:
        xyz = Path(tmp) / "h2.xyz"
        xyz.write_text("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
        proc = subprocess.run(
            [candidate, "h2.xyz", "--opt", "--gfn", "2"],
            cwd=tmp,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        if proc.returncode == 0 and (Path(tmp) / "xtbopt.xyz").exists():
            return True, None
        detail = (proc.stderr or proc.stdout).strip().splitlines()
        summary = detail[0] if detail else f"exit code {proc.returncode}"
        return False, f"xtb binary at {candidate} failed --opt self-check ({summary}); install the conda-forge xtb build."


def executable_status(
    name: str,
    version_args: Sequence[str] = ("--version",),
    timeout_s: int = 5,
) -> ToolStatus:
    candidates = _candidate_executable_paths(name)
    if not candidates:
        env_name = TOOL_ENV_VARS.get(name, f"CODEX_CHEM_{name.upper()}_PATH")
        return ToolStatus(
            name=name,
            status="unavailable",
            message=f"{name} not found on PATH, Homebrew, bundled tool env, or {env_name}.",
        )
    cache_key = ( "::".join(candidates), tuple(version_args) )
    cached = _STATUS_CACHE.get(cache_key)
    if cached is not None:
        return cached
    last_error: str | None = None
    for candidate in candidates:
        try:
            proc = subprocess.run(
                [candidate, *version_args],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout_s,
            )
            version = extract_version_line(proc.stdout or proc.stderr)
            if name == "xtb":
                ok, reason = _xtb_opt_healthy(candidate, timeout_s=timeout_s)
                if not ok:
                    last_error = reason
                    continue
            result = ToolStatus(name=name, status="available", version=version, path=candidate, message=None)
            _STATUS_CACHE[cache_key] = result
            return result
        except subprocess.TimeoutExpired as exc:  # pragma: no cover - slow binary
            result = ToolStatus(name=name, status="available", version=None, path=candidate, message=str(exc))
            _STATUS_CACHE[cache_key] = result
            return result
        except OSError as exc:  # PermissionError, Bad CPU type, missing interpreter, etc.
            last_error = str(exc)
            continue
        except Exception as exc:  # pragma: no cover
            last_error = str(exc)
            continue
    result = ToolStatus(
        name=name,
        status="unavailable",
        version=None,
        path=candidates[0],
        message=last_error or f"{name} is present at {candidates[0]} but could not be executed.",
    )
    _STATUS_CACHE[cache_key] = result
    return result


def run_command(
    args: Sequence[str],
    timeout_s: int = 60,
    cwd: str | None = None,
) -> tuple[int, str, str]:
    resolved_args = list(args)
    if resolved_args:
        path = executable_path(resolved_args[0])
        if path:
            resolved_args[0] = path
    proc = subprocess.run(
        resolved_args,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_s,
        cwd=cwd,
    )
    return proc.returncode, proc.stdout, proc.stderr


def split_nonempty(value: str | None, sep: str = ".") -> list[str]:
    if not value:
        return []
    return [part for part in value.split(sep) if part]
