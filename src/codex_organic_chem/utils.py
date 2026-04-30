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


def bundled_tool_prefix() -> Path:
    return Path.home() / ".local" / "share" / "codex-organic-chem" / "conda-tools"


def executable_path(name: str) -> str | None:
    env_name = TOOL_ENV_VARS.get(name, f"CODEX_CHEM_{name.upper()}_PATH")
    env_value = os.environ.get(env_name)
    if env_value and Path(env_value).expanduser().exists():
        return str(Path(env_value).expanduser())
    bundled = bundled_tool_prefix() / "bin" / name
    if name in {"xtb", "crest"} and bundled.exists():
        return str(bundled)
    found = shutil.which(name)
    if found:
        return found
    if bundled.exists():
        return str(bundled)
    return None


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


def executable_status(
    name: str,
    version_args: Sequence[str] = ("--version",),
    timeout_s: int = 5,
) -> ToolStatus:
    path = executable_path(name)
    if not path:
        env_name = TOOL_ENV_VARS.get(name, f"CODEX_CHEM_{name.upper()}_PATH")
        return ToolStatus(
            name=name,
            status="unavailable",
            message=f"{name} not found on PATH, bundled tool env, or {env_name}.",
        )
    version = None
    message = None
    try:
        proc = subprocess.run(
            [path, *version_args],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        version = extract_version_line(proc.stdout or proc.stderr)
    except Exception as exc:  # pragma: no cover - depends on local tools
        message = str(exc)
    return ToolStatus(name=name, status="available", version=version, path=path, message=message)


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
