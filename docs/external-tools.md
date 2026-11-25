# Optional external tools

[Back to README](../README.md)

## Optional External Tools

The code detects these tools if present and reports `unavailable` otherwise. Run:

```bash
uv run codex-chem doctor
```

On macOS, install the supported local tools with:

```bash
scripts/install_external_tools_macos.sh
```

- RDKit: core parser/drawer/computation, installed by `uv sync`.
- Open Babel: install with `brew install open-babel`; executable is `obabel`.
- xTB: install with `brew tap grimme-lab/qc && brew install xtb`, or use the
  conda-forge fallback at `~/.local/share/codex-organic-chem/conda-tools`; used
  by `chem_compute` tasks `xtb_opt`, `xtb_reactivity`, and `xtb_thermo` for
  optimized geometry/energy, atom-mapped charges and Fukui indices, and
  frequencies/free-energy corrections.
- CREST: install with `brew install crest` from the Grimme tap, or use the
  conda-forge fallback created at
  `~/.local/share/codex-organic-chem/conda-tools`; `chem_compute task=crest`
  runs the conformer search and returns ensemble energies and populations.
- For `xtb`/`crest`, the software prefers `CODEX_CHEM_*_PATH`, then the local conda-forge tool env, then PATH. This avoids broken Homebrew builds when the fallback exists.
- OSRA: optional OCSR fallback. On Apple Silicon macOS, the tool auto-detects the Rosetta/x86_64 conda install at `~/.local/share/codex-organic-chem/ocsr-tools/osra-osx64/bin/osra`; set `CODEX_CHEM_OSRA_PATH` only to override it.
- ChemSchematicResolver: optional legacy schematic/label resolver. It is installed locally but kept off the default crop OCSR path because the Q8 benchmark produced no structure candidates; set `CODEX_CHEM_CSR_CMD` or `CODEX_CHEM_ENABLE_CSR_DEFAULT=1` only when label-heavy schemes need it.
- MolScribe/DECIMER/MolGrapher/OpenChemIE/RxnScribe: optional modern image OCR adapters. Install them in separate ML environments and set the matching `CODEX_CHEM_*_CMD` variable.
- `scripts/molscribe_adapter.py` is a ready-to-use MolScribe command adapter for cropped molecule images. It prints JSON/SMILES for `codex-chem parse-image`, expands common textbook labels such as Me/OMe/NO2/CO2Et/t-Bu/Ph/OMs/OTs/NHTs/MEMO/SO2Tol, preserves R/X placeholders as dummy atoms for manual or context-aware resolution, preserves MolScribe chiral atom tags when possible, and downloads the upstream checkpoint on first use when `CODEX_CHEM_MOLSCRIBE_CHECKPOINT` is not set.
- `scripts/decimer_adapter.py`, `scripts/molgrapher_adapter.py`, `scripts/openchemie_adapter.py`, and `scripts/csr_adapter.py` provide JSON-line wrappers for local DECIMER, MolGrapher, OpenChemIE, and ChemSchematicResolver environments.
- Publication mechanism polish: ChemDraw, ChemDoodle, or ChemAxon Marvin are the
  chemically aware editors to learn from; Inkscape (`brew install --cask
  inkscape`) or Illustrator can polish the exported SVG but do not validate
  chemistry. `codex-chem doctor` lists this in `publication_figure_guide`.

Check local editor support with:

```bash
uv run codex-chem figure-tools
```

This reports whether ChemDraw, ChemDoodle, Marvin, Ketcher, Inkscape,
Illustrator, and Open Babel CDXML conversion are available. ChemDraw support is
file/workflow integration: the tool prepares SVG/spec/CDXML packages and gives
open-command templates, but does not assume unattended ChemDraw GUI automation.

Custom OCSR command variables may contain `{input}` and should print JSON with a
`smiles`, `canonical_smiles`, `reaction_smiles`, or `molfile` field, a JSON
`candidates` list, or plain SMILES lines. For line-angle/skeletal molecule
crops, `codex-chem parse-image` also tries white-background, padded,
high-contrast, and thickened line-art image variants before ranking candidates.
Configured adapters may therefore be called more than once per crop and should
be idempotent. The winning candidate records `metadata.image_variant` so poor
raw-crop recognition can be traced back to a preprocessing variant. Candidate
ranking still penalizes collapsed tiny outputs and low-confidence label-artifact
atoms, and preserves raw tool metadata/warnings instead of forcing a single
merged answer.

Scheme-level OCSR helpers:

```bash
uv run codex-chem crop-image scheme.png --grid --boxes '[{"id": "1", "box": [0.05, 0.10, 0.32, 0.45]}]' --output-dir ocsr_crops
uv run codex-chem ocsr-benchmark --gold-smiles render_q8_sample_style.py --image-dir ocsr_crops --out tsv
uv run codex-chem parse-scheme --image scheme.png --crops ocsr_crops --gold-map legend.json
```

`crop-image` writes 2× upscaled sub-crops and an optional normalized `0.1`-step
`grid_overlay.png` ruler image for dense pages. `ocsr-benchmark` reports exact
graph match, main-fragment match, dummy count, fragment count, stereo warnings,
and tool warnings for golden crop sets. `parse-scheme` accepts either a crop
directory or an inline JSON list of bounding boxes, keeps raw candidates,
applies optional `R`/`X`/`R1` definitions with stereo-preserving RDKit graph
replacement, and emits route-consistency warnings for likely OCR jumps such as
arene ortho/meta/para changes.

Local Q8 benchmark note, 2026-05-10/2026-05-11: DECIMER, MolGrapher, OpenChemIE,
OSRA, and ChemSchematicResolver were installed or tested against the 9 crop
golden set from `render_q8_sample_style.py`. None of the tested open-source
tools produced an exact graph or main-fragment match on this difficult scheme.
MolScribe was the most useful crop-level default because it preserved larger
scaffolds and dummy placeholders for `R`/`X`/`R1`/`R2`; DECIMER sometimes
recognized the core but emitted non-RDKit tokens such as `[R]`/`[X]`;
MolGrapher was useful for OCR label metadata but often collapsed structures or
read labels as atoms; the OpenChemIE crop adapter usually duplicated MolScribe
behavior or returned `<invalid>`. OpenChemIE reaction mode did segment one full
scheme into reactants/products/conditions with bbox/molfile data, so it remains
useful as a scheme-level auxiliary source, but the output still needs conversion
and manual review before reaction SMILES use. OSRA is available as a Rosetta
x86_64 fallback and produced RDKit-sanitizable dummy-atom skeletons on all 9
crops, but no exact/main-fragment matches. ChemSchematicResolver is available
through a local Python 3.6 x86_64 environment with a compiled `osra_rgroup`
extension; on these isolated crops it produced no structure candidates, so it
remains best treated as a label/schematic auxiliary rather than the primary crop
recognizer.

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
export CODEX_CHEM_DISABLE_IMAGE_PREPROCESSING=1
export CODEX_CHEM_KEEP_PREPROCESSED_IMAGES=1
export CODEX_CHEM_CHEMDRAW_APP="/Applications/ChemDraw.app"
export CODEX_CHEM_CHEMDOODLE_APP="/Applications/ChemDoodle.app"
export CODEX_CHEM_MARVIN_APP="/Applications/MarvinSketch.app"
export CODEX_CHEM_KETCHER_URL="http://localhost:8080"
export CODEX_CHEM_KETCHER_DIST="/path/to/ketcher/build"
```
