# Installing codex-organic-chem

Written so that a coding agent (Claude Code, Codex) or a person can set the
repository up by following it in order. Every step ends with a check. Run all
commands from the repository root.

## 1. uv

The only prerequisite is [uv](https://docs.astral.sh/uv/). It provides Python
3.12 (the project requires `>=3.12,<3.13`) and installs every Python package at
the version pinned in `uv.lock`.

```bash
uv --version
```

If uv is missing, install it with `brew install uv` (macOS) or
`curl -LsSf https://astral.sh/uv/install.sh | sh`. An agent should ask the user
before installing system software.

## 2. Python environment

```bash
uv sync --locked --extra mcp
```

- `--locked` installs exactly what `uv.lock` pins. If uv says the lock file is
  out of date, `pyproject.toml` was changed without re-locking: run `uv lock`,
  then this command again.
- `--extra mcp` adds the MCP server. Add `--extra dev` to run the tests.
- Do not install into another interpreter with pip. If pip is unavoidable, use
  Python 3.12 and `pip install -e ".[mcp]"`; the extra keeps `mcp` below 2.0.

## 3. Check the install

```bash
uv run python scripts/check_install.py
```

The last line must be `install OK`. The script runs `codex-chem doctor`, asks
the CLI a question with a known answer, and starts the MCP server exactly as
`.mcp.json` does. For every `FAIL`, run the command shown on its `next:` line
and check again.

## 4. Connect the MCP server and the skill

**Claude Code.** `.mcp.json` in the repository root starts the server with
`uv run --extra mcp codex-chem-mcp`. Approve the server when Claude Code asks.
If the session started before step 2, reconnect the server (`/mcp`) or start a
new session. The skill in `.claude/skills/organic-chemistry-assistant` is picked
up automatically.

**Codex.** Add the server to `~/.codex/config.toml`, with the absolute path of
this repository:

```toml
[mcp_servers.codex-organic-chem]
command = "uv"
args = ["run", "--directory", "/absolute/path/to/this/repository", "--extra", "mcp", "codex-chem-mcp"]
```

The skill is in `.agents/skills/organic-chemistry-assistant`.

**Other hosts.** Use the same command and arguments as for Codex, and copy the
skill folder into the host's skills directory. Keep one copy of the skill per
host: an older copy under the same name is listed next to this one, and the
model may load either.

If the host reports that `uv` cannot be found (desktop apps do not always see
the shell's `PATH`), replace `"uv"` with the absolute path printed by
`which uv`.

## 5. Optional engines

None of these is needed for the core. `codex-chem doctor` lists which are
missing, and the tools that need them say so instead of guessing.

| Engine | Used for | Install |
|---|---|---|
| xTB, CREST | energies, Fukui indices, conformer ensembles | macOS: `scripts/install_external_tools_macos.sh` (Homebrew; also installs Open Babel). An agent should ask the user first. |
| Open Babel | file-format conversion | same script, or `brew install open-babel` |
| OCSR models (MolScribe, DECIMER, MolGrapher, OpenChemIE, OSRA) | reading structures from images | separate environments; see [docs/external-tools.md](docs/external-tools.md) |
| Ketcher | drawing editor inside the review page | `cd integrations/ketcher && npm install && npm run build`; without it the review page uses its built-in editor |

## 6. Tests (for development)

```bash
uv sync --locked --extra mcp --extra dev
uv run pytest
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| MCP server fails with `No module named 'mcp'` | environment synced without the `mcp` extra | `uv sync --locked --extra mcp` |
| MCP server fails with `No module named 'mcp.server.fastmcp'` | mcp 2.x installed, usually by pip without the lock file | `uv sync --locked --extra mcp`, or `pip install "mcp>=1.2,<2"` |
| pip refuses the package because of `Requires-Python` | an interpreter other than 3.12 | use uv, or a Python 3.12 environment |
| The MCP server times out the first time it starts | its first `uv run` is still downloading packages | run step 2 first, then reconnect |
| `codex-chem: command not found` | the project environment is not activated | prefix the command with `uv run`; from another directory use `uv run --directory /path/to/repository codex-chem ...` |
| No browser opens for a review page | headless machine or no default browser | the command prints the URL; `CODEX_CHEM_REVIEW_OPEN_BROWSER=0` stops it from trying |
| xTB, CREST or OCSR reported unavailable | optional engines not installed | step 5 |
