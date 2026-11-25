"""Check that codex-organic-chem is installed and wired up correctly.

    uv run python scripts/check_install.py

Runs `codex-chem doctor`, one CLI command whose answer is known, and a real MCP
handshake with the server command from .mcp.json (the one Claude Code starts).
Exits 0 when the core install and the MCP server work. Optional engines (xTB,
CREST, Open Babel, OCSR models) are reported but not required.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.executable).parent
results: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, ok, detail))
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{': ' + detail if detail else ''}", flush=True)


def console_script(name: str) -> str | None:
    return shutil.which(name, path=str(BIN)) or shutil.which(name)


def run_json(args: list[str], timeout: int = 300) -> dict:
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or proc.stdout).strip()[-400:])
    return json.loads(proc.stdout)


def check_doctor(cli: str) -> None:
    report = run_json([cli, "doctor"])
    summary, env = report["summary"], report["python_environment"]
    record("Python 3.12 and required packages", summary["python_env_ready"],
           ", ".join(f"{k} {v.get('version') or v['status']}" for k, v in env.items() if k in ("python", "rdkit", "numpy", "pillow")))
    record("review page files", env["review_ui"]["status"] == "available")
    record("mcp package for the MCP server", summary["mcp_server_ready"], env["mcp"].get("version") or env["mcp"]["fix"])
    optional = [("xTB + CREST", summary["calculation_stack_ready"]), ("OCSR engine", summary["image_ocr_stack_ready"]),
                ("Open Babel", report["statuses"]["obabel"]["status"] == "available")]
    print("  optional engines: " + ", ".join(f"{name} {'yes' if ok else 'no'}" for name, ok in optional))
    for step in summary["next_steps"]:
        print(f"    next: {step}")


def check_cli(cli: str) -> None:
    out = run_json([cli, "compare", "--a", "C[C@H](O)CC", "--b", "C[C@@H](O)CC"])
    record("CLI answers a known question", out.get("relation") == "enantiomers", f"relation = {out.get('relation')}")


def check_mcp() -> None:
    import os

    server = json.loads((ROOT / ".mcp.json").read_text(encoding="utf-8"))["mcpServers"]["codex-organic-chem"]
    search_path = os.pathsep.join(
        [os.environ.get("PATH", ""), "/opt/homebrew/bin", "/usr/local/bin", str(Path.home() / ".local" / "bin")]
    )
    resolved_cmd = shutil.which(server["command"]) or shutil.which(server["command"], path=search_path) or server["command"]
    command = [resolved_cmd, *server.get("args", [])]
    proc = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1)
    lines: list[str] = []
    errors: list[str] = []
    threading.Thread(target=lambda: lines.extend(proc.stdout), daemon=True).start()
    threading.Thread(target=lambda: errors.extend(proc.stderr), daemon=True).start()

    def send(message: dict) -> None:
        proc.stdin.write(json.dumps(message) + "\n")
        proc.stdin.flush()

    def reply(message_id: int, timeout: float = 180) -> dict | None:
        end = time.time() + timeout
        while time.time() < end and proc.poll() is None:
            for line in list(lines):
                try:
                    message = json.loads(line)
                except ValueError:
                    continue
                if message.get("id") == message_id:
                    return message
            time.sleep(0.1)
        return None

    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "check_install", "version": "1"}}})
        init = reply(1)
        if init is None:
            record("MCP server starts (as in .mcp.json)", False, "".join(errors).strip()[-400:] or "no answer")
            return
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        tools = (reply(2) or {}).get("result", {}).get("tools", [])
        record("MCP server starts (as in .mcp.json)", len(tools) >= 26, f"{len(tools)} tools")
    finally:
        proc.kill()


def main() -> int:
    print(f"codex-organic-chem install check in {ROOT}")
    cli = console_script("codex-chem")
    if cli is None:
        record("codex-chem command", False, "not found; run `uv sync --locked --extra mcp` in the repository")
    else:
        for check in (check_doctor, check_cli):
            try:
                check(cli)
            except Exception as exc:  # report and keep going
                record(check.__name__.replace("check_", ""), False, str(exc))
    try:
        check_mcp()
    except Exception as exc:
        record("MCP server starts (as in .mcp.json)", False, str(exc))
    failed = [name for name, ok, _ in results if not ok]
    print("install OK" if not failed else f"{len(failed)} problem(s): " + "; ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
