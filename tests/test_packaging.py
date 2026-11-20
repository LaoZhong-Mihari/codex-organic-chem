"""Install-time guarantees: declared dependencies, the MCP extra, the project
MCP config and the environment report of `codex-chem doctor`."""

import ast
import json
import re
import sys
import tomllib
from pathlib import Path

from codex_organic_chem import external

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
IMPORT_TO_DIST = {"PIL": "pillow"}


def _dist(requirement: str) -> str:
    return re.split(r"[<>=!~\[; ]", requirement, maxsplit=1)[0].strip().lower()


def test_third_party_imports_are_declared():
    declared = {_dist(r) for r in PYPROJECT["project"]["dependencies"]}
    extras = {_dist(r) for reqs in PYPROJECT["project"]["optional-dependencies"].values() for r in reqs}
    used = set()
    for path in (ROOT / "src" / "codex_organic_chem").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                used |= {alias.name.split(".")[0] for alias in node.names}
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                used.add(node.module.split(".")[0])
    third_party = {name for name in used if name not in sys.stdlib_module_names and name != "codex_organic_chem"}
    missing = {name for name in third_party if IMPORT_TO_DIST.get(name, name).lower() not in declared | extras}
    assert not missing, f"imported but not declared in pyproject.toml: {sorted(missing)}"


def test_mcp_extra_stays_on_v1():
    (requirement,) = [r for r in PYPROJECT["project"]["optional-dependencies"]["mcp"] if _dist(r) == "mcp"]
    assert "<2" in requirement.replace(" ", ""), "codex-chem-mcp uses FastMCP, which mcp 2.x renamed"


def test_project_mcp_config_installs_the_mcp_extra():
    server = json.loads((ROOT / ".mcp.json").read_text(encoding="utf-8"))["mcpServers"]["codex-organic-chem"]
    args = server["args"]
    assert server["command"] == "uv" and args[0] == "run" and args[-1] == "codex-chem-mcp"
    assert args[args.index("--extra") + 1] == "mcp", "a plain `uv sync` would leave the server without mcp"


def test_doctor_reports_the_python_environment():
    report = external.doctor_report()
    env = report["python_environment"]
    assert {"python", "rdkit", "numpy", "pillow", "mcp", "review_ui"} <= set(env)
    assert report["summary"]["python_env_ready"] is True
    assert env["review_ui"]["status"] == "available"


def test_doctor_names_the_fix_when_mcp_is_missing(monkeypatch):
    real = external._dist_version
    monkeypatch.setattr(external, "_dist_version", lambda dist: None if dist == "mcp" else real(dist))
    report = external.doctor_report()
    assert report["summary"]["mcp_server_ready"] is False
    assert "mcp: uv sync --extra mcp" in report["summary"]["next_steps"]
