# Dependencies

This repository contains three related parts:

- the `organic-chemistry-assistant` Codex skill in `.agents/skills/organic-chemistry-assistant`,
- the `codex-organic-chem` Python CLI/MCP package in `src/codex_organic_chem`,
- optional local editor integrations under `integrations/`.

The pinned dependency graphs live in `uv.lock` for Python and
`integrations/ketcher/package-lock.json` for the Ketcher frontend.

## Required runtime

- Python `>=3.12,<3.13`
- `uv` for environment creation, locking, and running commands
- RDKit `>=2024.9.1`

Install and check the core toolchain:

```bash
uv venv --python 3.12
uv sync
uv run codex-chem doctor
```

## Python extras

The project defines optional extras in `pyproject.toml`:

- `mcp`: installs `mcp>=1.2.0` for the `codex-chem-mcp` server.
- `dev`: installs `pytest>=8.2.0` for the test suite.

For normal development and MCP use:

```bash
uv sync --extra dev --extra mcp
uv run pytest
```

The build backend is `hatchling>=1.25`.

## Codex skill files

The skill itself has no separate package manager. It depends on this repository's
CLI/MCP tools being available in the same environment:

- `.agents/skills/organic-chemistry-assistant/SKILL.md`
- `.agents/skills/organic-chemistry-assistant/agents/openai.yaml`
- `.mcp.json`, which launches `uv run codex-chem-mcp`

To install the skill into another Codex setup, copy or symlink the
`.agents/skills/organic-chemistry-assistant` folder into that setup's skills
directory and make sure `codex-chem` / `codex-chem-mcp` run from this project.

## Optional external chemistry tools

The CLI works without these tools. Keep them out of the core environment unless
the task needs them; `codex-chem doctor` reports what is configured locally.

- Open Babel: executable `obabel`; used for structure conversion and CDXML/SDF
  workflows. macOS: `brew install open-babel`.
- xTB: executable `xtb`; used by `chem_compute --task xtb_opt`. macOS:
  `brew tap grimme-lab/qc && brew install xtb`.
- CREST: executable `crest`; used by `chem_compute --task crest`. macOS:
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
- Ketcher: open-source structure/reaction correction gate. Use the bundled
  Vite integration or set `CODEX_CHEM_KETCHER_URL` /
  `CODEX_CHEM_KETCHER_DIST`.
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

## Network access

Most tools are local. `chem_literature_search` uses the Crossref Works API
(`https://api.crossref.org/works`) for literature metadata when requested.

## Licenses and third-party files

- This project's original source is MIT licensed; see `LICENSE`.
- ChemDoodle Web Components files under `integrations/chemdoodle/` are
  third-party iChemLabs files with GPLv3/commercial licensing terms in their
  file headers. See `THIRD_PARTY_NOTICES.md`.
- Proprietary editor packages such as ChemDraw, ChemDoodle commercial builds,
  and ChemAxon Marvin are not redistributable here; configure local paths with
  environment variables.
