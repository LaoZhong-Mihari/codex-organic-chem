#!/usr/bin/env bash
set -euo pipefail

TOOL_PREFIX="${CODEX_CHEM_TOOL_PREFIX:-$HOME/.local/share/codex-organic-chem/conda-tools}"

if ! command -v brew >/dev/null 2>&1; then
  echo "Homebrew is required on macOS: https://brew.sh" >&2
  exit 1
fi

echo "Installing Open Babel and xTB via Homebrew..."
brew tap grimme-lab/qc >/dev/null
brew install open-babel xtb micromamba

echo "Installing CREST. Homebrew is tried first; conda-forge is used as fallback."
if brew install crest; then
  echo "CREST installed via Homebrew."
else
  echo "Homebrew CREST failed; installing CREST into $TOOL_PREFIX with micromamba."
  micromamba create -y -p "$TOOL_PREFIX" -c conda-forge crest xtb
fi

cat <<'EOF'

Optional image OCR adapters still need separate setup:
- OSRA: on Apple Silicon macOS, the tested fallback is the osx-64 conda package at ~/.local/share/codex-organic-chem/ocsr-tools/osra-osx64; codex-chem auto-detects it.
- MolScribe: install in a separate ML environment and set CODEX_CHEM_MOLSCRIBE_CMD='your-command {input}'.
- DECIMER/MolGrapher/OpenChemIE: install each in a separate ML environment and set the matching CODEX_CHEM_*_CMD variable.
- ChemSchematicResolver: installed/tested as a Rosetta/x86_64 Python 3.6 env at ~/.local/share/codex-organic-chem/ocsr-tools/csr-osx64 plus scripts/csr_adapter.py, but off by default; enable with CODEX_CHEM_CSR_CMD or CODEX_CHEM_ENABLE_CSR_DEFAULT=1.
- RxnScribe: install in a separate ML environment and set CODEX_CHEM_RXNSCRIBE_CMD='your-command {input}'.

Run:
  codex-chem doctor
EOF
