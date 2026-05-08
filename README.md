# Codex Organic Chemistry Assistant

Local open-source toolchain for helping Codex reason about organic chemistry.
It gives Codex deterministic tools for:

- molecule/reaction image parsing with OCSR adapters and human-review warnings,
- RDKit structure normalization, drawing, descriptors, conformers, and charges,
- optional Open Babel, xTB, and CREST integration when binaries are installed,
- rule-based reaction sanity checks, retrosynthesis hints, condition families, and mechanism drafts,
- a Codex Skill and MCP server exposing the same workflow.

This is a research assistant, not an experimental safety authority. Treat all
synthetic, mechanistic, and computational suggestions as hypotheses requiring
literature review and expert validation.

## Repository Contents

- `.agents/skills/organic-chemistry-assistant`: the Codex skill instructions and
  UI metadata.
- `src/codex_organic_chem`: the local Python CLI/MCP implementation used by the
  skill.
- `integrations/ketcher`: optional Ketcher review UI source; build artifacts are
  intentionally not committed.
- `integrations/chemdoodle`: optional ChemDoodle Web Components viewer files.

See [DEPENDENCIES.md](DEPENDENCIES.md) for the complete dependency map,
including Python, MCP, optional chemistry binaries, editor integrations, npm
packages, environment variables, and third-party licensing notes.

## Quick Start

```bash
uv venv --python 3.12
uv sync --extra dev --extra mcp
uv run codex-chem doctor
uv run codex-chem figure-tools
uv run codex-chem normalize --smiles "CC(=O)Oc1ccccc1C(=O)O"
uv run codex-chem draw --smiles "c1ccccc1C(=O)O" --output svg
uv run codex-chem compute --smiles "CCO" --task descriptors --task conformers
uv run codex-chem reaction-analyze --reaction "CBr.[OH-]>>CO.[Br-]" --mode sanity_check
uv run codex-chem mechanism-draft --reaction "CBr.[OH-]>>CO.[Br-]"
```

## Codex Skill Setup

The skill folder is packaged with the project:

```text
.agents/skills/organic-chemistry-assistant/
├── SKILL.md
└── agents/openai.yaml
```

For another Codex environment, copy or symlink that folder into the target
skills directory and keep this repository available so `.mcp.json` can launch:

```bash
uv run codex-chem-mcp
```

## MCP

The repository includes `.mcp.json`:

```bash
uv run codex-chem-mcp
```

The MCP tools have the planned names:

- `chem_parse_image`
- `chem_normalize_structure`
- `chem_draw`
- `chem_compute`
- `chem_reaction_analyze`
- `chem_mechanism_draft`
- `chem_mechanism_render`
- `chem_mechanism_spec_example`
- `chem_figure_tool_status`
- `chem_input_review`
- `chem_literature_search`
- `chem_synthesis_suggest`

## Required Review Gate

For any molecule or reaction recognized from an image, scan, screenshot, or
ambiguous text, the assistant must render the parsed candidate and wait for
explicit user confirmation before continuing.

```bash
codex-chem input-review --image-path scheme.png --kind reaction
codex-chem input-review --smiles "CC(=O)Oc1ccccc1C(=O)O"
```

Downstream publication mechanisms, synthesis suggestions, expensive
calculations, and literature-backed route expansion should only run after the
user confirms the rendered structure or supplies corrected machine-readable
input.

## Mechanism Figures

Publication-oriented mechanism work now has two levels. `mechanism-draft` gives
a mechanistic hypothesis and explanation, but it is not treated as final
manuscript artwork:

```bash
codex-chem mechanism-draft \
  --reaction "CBr.[OH-]>>CO.[Br-]" \
  --quality publication
```

This returns `awaiting_user_confirmation` with a reaction preview. After the
rendered reaction is confirmed:

```bash
codex-chem mechanism-draft \
  --reaction "CBr.[OH-]>>CO.[Br-]" \
  --quality publication \
  --confirmed
```

For accurate figure output, use an explicit atom-mapped mechanism canvas:

```bash
codex-chem mechanism-render --example
codex-chem mechanism-render --spec mechanism.spec.json --output-dir mechanism_out
```

The spec must provide every intermediate and every curved arrow anchor using
atom-map numbers or atom-map bond pairs. It supports:

- `spec_version: "2.0"` semantic mechanism specs,
- intermediate molecules as mapped SMILES/Molfile,
- curved electron-pair and radical arrows,
- lone-pair dots, with publication defaults that show only mechanistically
  relevant electron-source lone pairs unless `lone_pair_display="all"` is set,
- formal charges and custom charge labels,
- partial charges such as `δ+` and `δ-`,
- graph edits for bond formation/breaking, charge changes, proton transfers,
  resonance, radicals, and counterion bookkeeping,
- panel titles, notes, captions, SVG output, and ChemDraw-like CDXML output.

When `--output-dir` is used, the result includes `editor_followup` with local
ChemDraw/ChemDoodle/Marvin/Inkscape/Illustrator open-command templates for the
generated files when those tools are detected.

`mechanism-render` now returns and writes a multi-format publication/review
package:

- `mechanism.svg`: RDKit-based vector layout with atom-map anchored electron
  arrows, shortened endpoints, restrained journal-style arrow weights, optional
  source-only lone-pair display, and cross-panel alignment.
- `mechanism.cdxml`: ChemDraw-like CDXML containing fragments, nodes, bonds,
  text, reaction-step records, and electron-flow arrow objects.
- `mechanism.chemdoodle.json`: ChemDoodle JSON with `Pusher` shapes whose
  source/target ids point to atom or bond objects and preserve electron count;
  atoms also carry lone-pair counts such as `numLonePair`.
- `mechanism.chemdoodle.html`: a local ChemDoodle Web Components viewer when
  `integrations/chemdoodle/ChemDoodleWeb.js` is installed.
- `mechanism.trace.json`: an AI-readable step trace with starting state,
  electron moves, graph edits, bond/charge changes, expected state, validation
  checks, and review boundaries.
- `mechanism.ketcher.json`: an adapter payload for using Ketcher as the
  open-source correction gate for structures, reactions, and atom maps.
- `mechanism.marvin.json`: an adapter payload for ChemAxon Marvin/Marvin JS
  review when a licensed Marvin package is installed locally.

`mechanism-render` also runs a publication-figure gate. If an arrow has no
electron source/sink, starts from an atom center instead of a lone pair or bond,
references missing atom-map numbers, or depends on a hidden lone pair, the result
is marked `blocked_for_publication` or `ok_with_warnings` and includes
`publication_checks`.

This mirrors how manuscript mechanisms are usually made: ChemDraw/CDXML or
similar chemical drawing software for the chemical canvas, with final manual
polishing in ChemDraw, ChemDoodle, Marvin, ChemSketch, Illustrator, or Inkscape.
ChemDraw/CDXML, ChemDoodle, and Marvin-style tools support reaction/mechanism
arrows and chemical annotations in ways RDKit alone does not. The generated
`mechanism_trace` is intended to help an AI reason step by step through the
public, auditable mechanism state rather than guessing from a single overall
reaction equation.

Implementation rule learned from ChemDoodle/Marvin-style editors: electron
pushers must be attached to chemical objects, not arbitrary SVG curves.
Therefore `mechanism-render` treats `from_lone_pair_atom_map`, `from_bond`,
`to_atom_map`, and `to_bond` as chemical anchors. It warns when disconnected
reactants are packed into one molecule entry, because that gives poor RDKit
layout and unreadable arrow geometry; use separate molecule entries for each
reactant, intermediate, counterion, or byproduct.

For manuscript-like figures, keep the visual electron annotations sparse:
`presentation_mode="publication"` defaults to source-only lone-pair display, so
an anionic nucleophile can carry one visible lone-pair pair at the arrow tail
while nonreacting heteroatom lone pairs and product lone pairs stay hidden.
Use full lone-pair counts only for teaching, debugging, or when the chemistry
requires them.

### Optional editor integrations

Ketcher can be built locally from the bundled integration:

```bash
cd integrations/ketcher
npm install
npm run build
npm run preview -- --port 4173
```

`codex-chem figure-tools` auto-detects `integrations/ketcher/dist` after a
local build. The repository commits `package.json` and `package-lock.json`, but
ignores `node_modules/` and `dist/`.

ChemDoodle Web Components are detected from `integrations/chemdoodle` or
`CODEX_CHEM_CHEMDOODLE_WEB_DIR`; the renderer writes
`mechanism.chemdoodle.html` so atom-bound lone pairs and object-bound pusher
arrows can be reviewed in a ChemDoodle-style renderer. The bundled ChemDoodle
Web Components files are third-party GPLv3/commercial-license files; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Marvin JS/Desktop requires a ChemAxon license. Place a licensed Marvin JS
package in `integrations/marvin` or set `CODEX_CHEM_MARVIN_JS_DIR`; the renderer
still writes `mechanism.marvin.json` and CDXML/MOL/RXN-facing artifacts for
manual Marvin review.

## Stepwise Synthesis Suggestions

Synthesis suggestions are intentionally one decision layer at a time. The tool
does not dump a full route in one pass.

```bash
codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O"
codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O" --confirmed
```

Each option includes local tool data and, when enabled, Crossref literature
metadata. The next step is to choose or reject one disconnection before any
route expansion.

## Optional External Tools

The code detects these tools if present and reports `unavailable` otherwise. Run:

```bash
codex-chem doctor
```

On macOS, install the supported local tools with:

```bash
scripts/install_external_tools_macos.sh
```

- RDKit: core parser/drawer/computation, installed by `uv sync`.
- Open Babel: install with `brew install open-babel`; executable is `obabel`.
- xTB: install with `brew tap grimme-lab/qc && brew install xtb`, or use the conda-forge fallback at `~/.local/share/codex-organic-chem/conda-tools`; used by `chem_compute task=xtb_opt`.
- CREST: install with `brew install crest` from the Grimme tap, or use the conda-forge fallback created at `~/.local/share/codex-organic-chem/conda-tools`.
- For `xtb`/`crest`, the software prefers `CODEX_CHEM_*_PATH`, then the local conda-forge tool env, then PATH. This avoids broken Homebrew builds when the fallback exists.
- OSRA: optional OCSR fallback. macOS has no stable Homebrew formula in this setup; use Docker/source OSRA and set `CODEX_CHEM_OSRA_PATH` if needed.
- MolScribe/RxnScribe: optional modern image OCR adapters. Install them in separate ML environments and set `CODEX_CHEM_MOLSCRIBE_CMD` / `CODEX_CHEM_RXNSCRIBE_CMD`.
- `scripts/molscribe_adapter.py` is a ready-to-use MolScribe command adapter for cropped molecule images. It prints JSON/SMILES for `codex-chem parse-image`, expands common textbook labels such as Me/OMe/NO2/CO2Et/t-Bu/Ph/OMs/OTs/NHTs, preserves MolScribe chiral atom tags when possible, and downloads the upstream checkpoint on first use when `CODEX_CHEM_MOLSCRIBE_CHECKPOINT` is not set.
- Publication mechanism polish: ChemDraw, ChemDoodle, or ChemAxon Marvin are the
  chemically aware editors to learn from; Inkscape (`brew install --cask
  inkscape`) or Illustrator can polish the exported SVG but do not validate
  chemistry. `codex-chem doctor` lists this in `publication_figure_guide`.

Check local editor support with:

```bash
codex-chem figure-tools
```

This reports whether ChemDraw, ChemDoodle, Marvin, Ketcher, Inkscape,
Illustrator, and Open Babel CDXML conversion are available. ChemDraw support is
file/workflow integration: the tool prepares SVG/spec/CDXML packages and gives
open-command templates, but does not assume unattended ChemDraw GUI automation.

Custom OCSR command variables may contain `{input}` and should print JSON with a
`smiles` field, or plain SMILES lines.

Useful environment overrides:

```bash
export CODEX_CHEM_OBABEL_PATH=/path/to/obabel
export CODEX_CHEM_XTB_PATH=/path/to/xtb
export CODEX_CHEM_CREST_PATH=/path/to/crest
export CODEX_CHEM_OSRA_PATH=/path/to/osra
export CODEX_CHEM_MOLSCRIBE_CMD='your-molscribe-command --input {input}'
export CODEX_CHEM_RXNSCRIBE_CMD='your-rxnscribe-command --input {input}'
export CODEX_CHEM_CHEMDRAW_APP="/Applications/ChemDraw.app"
export CODEX_CHEM_CHEMDOODLE_APP="/Applications/ChemDoodle.app"
export CODEX_CHEM_MARVIN_APP="/Applications/MarvinSketch.app"
export CODEX_CHEM_KETCHER_URL="http://localhost:8080"
export CODEX_CHEM_KETCHER_DIST="/path/to/ketcher/build"
```
