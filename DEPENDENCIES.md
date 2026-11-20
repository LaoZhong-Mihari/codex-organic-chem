# Dependencies

This repository contains four related parts:

- the `organic-chemistry-assistant` Codex skill in `.agents/skills/organic-chemistry-assistant`,
- its synchronized Claude Code mirror in `.claude/skills/organic-chemistry-assistant`,
- the `codex-organic-chem` Python CLI/MCP package in `src/codex_organic_chem`,
- optional local editor integrations under `integrations/`.

The pinned dependency graphs live in `uv.lock` for Python and
`integrations/ketcher/package-lock.json` for the Ketcher frontend.

## Required runtime

- Python `>=3.12,<3.13` (uv provides it)
- `uv` for environment creation, locking, and running commands
- RDKit `>=2024.9.1`, NumPy `>=1.26`, Pillow `>=10.0`

Install and check the toolchain (step by step, with troubleshooting, in
[INSTALL.md](INSTALL.md)):

```bash
uv sync --locked --extra mcp
uv run python scripts/check_install.py
```

## Python extras

The project defines optional extras in `pyproject.toml`:

- `mcp`: installs `mcp>=1.2.0,<2` for the `codex-chem-mcp` server. The server
  uses `FastMCP`, which mcp 2.x renamed, so the extra stays on 1.x.
- `dev`: installs `pytest>=8.2.0` for the test suite.

For normal development and MCP use:

```bash
uv sync --extra dev --extra mcp
uv run pytest
```

The build backend is `hatchling>=1.25`.

## Codex and Claude Code skill files

The skill itself has no separate package manager. It depends on this repository's
CLI/MCP tools being available in the same environment:

- `.agents/skills/organic-chemistry-assistant/SKILL.md`
- `.agents/skills/organic-chemistry-assistant/agents/openai.yaml`
- `.claude/skills/organic-chemistry-assistant/SKILL.md`
- `.mcp.json`, which launches `uv run --extra mcp codex-chem-mcp` (the extra is
  installed on first start if the environment was synced without it)

Codex installs the `.agents` version into its skills directory. Claude Code
discovers the `.claude/skills` mirror while working in this repository and uses
the project-level `.mcp.json`. In either environment, make sure `codex-chem` /
`codex-chem-mcp` run from this project. A regression test prevents the two
`SKILL.md` files from drifting.

## Optional external chemistry tools

The CLI works without these tools. Keep them out of the core environment unless
the task needs them; `codex-chem doctor` reports what is configured locally.

- Open Babel: executable `obabel`; used for structure conversion and CDXML/SDF
  workflows. macOS: `brew install open-babel`.
- xTB: executable `xtb`; used by `chem_compute --task xtb_opt`,
  `--task xtb_reactivity`, and `--task xtb_thermo`. These run real GFN2-xTB
  optimization/energy, atom-mapped charge and Fukui calculations, and
  frequency/free-energy calculations. macOS:
  `brew tap grimme-lab/qc && brew install xtb`.
- CREST: executable `crest`; used by `chem_compute --task crest` for a real
  conformer search with ensemble energies and Boltzmann populations. macOS:
  `brew tap grimme-lab/qc && brew install crest`. The bundled installer falls
  back to a micromamba conda-forge environment at
  `~/.local/share/codex-organic-chem/conda-tools`.
- MolScribe: recommended molecule-crop OCSR adapter. Install in a separate ML
  environment and set `CODEX_CHEM_MOLSCRIBE_CMD`.
- OSRA: useful legacy fallback for molecule crops. On Apple Silicon macOS, the
  tested local fallback is
  `~/.local/share/codex-organic-chem/ocsr-tools/osra-osx64/bin/osra`; the code
  auto-detects that path.
- ChemSchematicResolver: installed/tested as a Rosetta x86_64 Python 3.6 helper
  for schematic labels, but off by default because the Q8 crop benchmark
  produced no structure candidates. Enable only with `CODEX_CHEM_CSR_CMD` or
  `CODEX_CHEM_ENABLE_CSR_DEFAULT=1`.
- DECIMER, MolGrapher, OpenChemIE, and RxnScribe: experimental/auxiliary OCSR
  adapters. Keep them in separate environments and configure the matching
  `CODEX_CHEM_*_CMD` only when a task benefits from them.

macOS helper:

```bash
scripts/install_external_tools_macos.sh
```

MolScribe can be wired to the bundled adapter without adding heavy ML packages
to the main RDKit environment:

```bash
uv venv --python 3.10 ~/.local/share/codex-organic-chem/molscribe-venv
uv pip install --python ~/.local/share/codex-organic-chem/molscribe-venv/bin/python MolScribe huggingface_hub
uv pip install --python ~/.local/share/codex-organic-chem/molscribe-venv/bin/python 'numpy<2' --reinstall
export CODEX_CHEM_MOLSCRIBE_CMD="$HOME/.local/share/codex-organic-chem/molscribe-venv/bin/python scripts/molscribe_adapter.py --input {input}"
```

The adapter downloads the upstream MolScribe checkpoint on first use unless
`CODEX_CHEM_MOLSCRIBE_CHECKPOINT` points to a local `.pth` checkpoint. Keep
molecule crops tight: full PDF pages that include text, arrows, conditions, and
check boxes are not reliable OCSR inputs.

## Optional editor and figure tools

These tools are used for human review and polishing of publication-oriented
chemistry figures. They are not required for the Python tests.

- ChemDraw: commercial Revvity Signals ChemDraw app. Set
  `CODEX_CHEM_CHEMDRAW_APP` if auto-detection does not find it.
- ChemDoodle Desktop or ChemDoodle Web Components: used for object-bound
  mechanism review. Set `CODEX_CHEM_CHEMDOODLE_APP` or
  `CODEX_CHEM_CHEMDOODLE_WEB_DIR`.
- ChemAxon Marvin / Marvin JS: proprietary ChemAxon software; not bundled. Set
  `CODEX_CHEM_MARVIN_APP` or `CODEX_CHEM_MARVIN_JS_DIR`.
- Ketcher: open-source structure editor, optional. The review page itself is a
  single dependency-free HTML file served by the Python package; when
  `integrations/ketcher/dist` exists (or `CODEX_CHEM_KETCHER_DIST` points at a
  build) the page gains a `Draw...` button that embeds Ketcher. Without it the
  page still offers click-to-flip stereo edits, alternatives and SMILES/Molfile
  entry.
- Inkscape or Adobe Illustrator: optional final vector-polish tools. They do
  not validate chemistry.

Useful environment overrides:

```bash
export CODEX_CHEM_OBABEL_PATH=/path/to/obabel
export CODEX_CHEM_XTB_PATH=/path/to/xtb
export CODEX_CHEM_CREST_PATH=/path/to/crest
export CODEX_CHEM_OSRA_PATH=/path/to/osra
export CODEX_CHEM_MOLSCRIBE_CMD='your-molscribe-command --input {input}'
export CODEX_CHEM_DECIMER_CMD='your-decimer-command --input {input}'
export CODEX_CHEM_MOLGRAPHER_CMD='your-molgrapher-command --input {input}'
export CODEX_CHEM_OPENCHEMIE_CMD='your-openchemie-command --input {input}'
export CODEX_CHEM_CSR_CMD='your-chemschematicresolver-command --input {input}'
export CODEX_CHEM_ENABLE_CSR_DEFAULT=1
export CODEX_CHEM_RXNSCRIBE_CMD='your-rxnscribe-command --input {input}'
export CODEX_CHEM_OCR_TIMEOUT_S=300
export CODEX_CHEM_CHEMDRAW_APP="/Applications/ChemDraw.app"
export CODEX_CHEM_CHEMDOODLE_APP="/Applications/ChemDoodle.app"
export CODEX_CHEM_CHEMDOODLE_WEB_DIR="/path/to/ChemDoodleWeb"
export CODEX_CHEM_MARVIN_APP="/Applications/MarvinSketch.app"
export CODEX_CHEM_MARVIN_JS_DIR="/path/to/marvin-js"
export CODEX_CHEM_KETCHER_URL="http://localhost:8080"
export CODEX_CHEM_KETCHER_DIST="/path/to/ketcher/build"
export CODEX_CHEM_REVIEW_OPEN_BROWSER=0   # return the review URL without opening a browser
export CODEX_CHEM_REVIEW=ask              # auto | ask | off; overrides `codex-chem review-mode`
export CODEX_CHEM_HOME="$HOME/.codex-chem" # review results, saved settings, review-server state
export CODEX_CHEM_ENABLE_PYTHON=1   # opt-in: MCP tool chem_python (arbitrary local code)
```

## Ketcher frontend integration

The optional Ketcher wrapper requires Node.js compatible with Vite 7:

- Node.js `^20.19.0 || >=22.12.0`
- npm
- direct npm dependencies from `integrations/ketcher/package.json`:
  `@vitejs/plugin-react`, `vite`, `typescript`, `react`, `react-dom`,
  `ketcher-react`, and `ketcher-standalone`

Build locally:

```bash
cd integrations/ketcher
npm install
npm run build
npm run preview -- --port 4173
```

`node_modules/` and `dist/` are intentionally ignored and should be rebuilt
from `package-lock.json`.

The review server uses only Python's standard library and listens on
`127.0.0.1`; every session has its own random token. CLI review commands hand
their session to a background copy of that server
(`python -m codex_organic_chem.review_daemon`, started on demand, stopped after
15 idle minutes, state in `$CODEX_CHEM_HOME/review-server.json`), so the page
outlives the command. The MCP server hosts its sessions in-process. Session
results are written to `$CODEX_CHEM_HOME/reviews/`.

## Network access

Most tools are local. `chem_literature_search` sends the query text to the
OpenAlex Works API (abstracts, citation counts, open-access links) and to the
Crossref Works API as fallback and supplement. No API key and no contact
address is sent. The verification kit (`compare`, `enumerate`, `route-check`,
`model-kit`) needs only RDKit and NumPy, which RDKit already depends on; xTB is
used for `model-kit` refinement when installed.

## Code execution through MCP (opt-in)

Hosts with a shell should run `uv run python` and import
`codex_organic_chem.kit`. For MCP-only hosts the server can expose
`chem_python`, which runs an arbitrary Python snippet as the local user with a
timeout. It is registered only when the server starts with
`CODEX_CHEM_ENABLE_PYTHON=1`:

```json
{"mcpServers": {"codex-organic-chem": {"command": "uv", "args": ["run", "codex-chem-mcp"],
  "env": {"CODEX_CHEM_ENABLE_PYTHON": "1"}}}}
```

Leave it off when the server may be driven by content you do not trust.

## Licenses and third-party files

- This project's original source is MIT licensed; see `LICENSE`.
- ChemDoodle Web Components files under `integrations/chemdoodle/` are
  third-party iChemLabs files with GPLv3/commercial licensing terms in their
  file headers. See `THIRD_PARTY_NOTICES.md`.
- Proprietary editor packages such as ChemDraw, ChemDoodle commercial builds,
  and ChemAxon Marvin are not redistributable here; configure local paths with
  environment variables.
