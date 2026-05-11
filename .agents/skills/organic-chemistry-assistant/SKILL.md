---
name: organic-chemistry-assistant
description: Use this skill when Codex needs to validate, draw, compute, parse, or explain organic molecules, reaction schemes, retrosynthesis ideas, or mechanisms with the local codex-organic-chem toolchain.
metadata:
  short-description: Local organic chemistry reasoning tools for Codex
---

# Organic Chemistry Assistant

Use `codex-chem` or the `codex-organic-chem` MCP tools before making nontrivial
chemistry claims. Treat tool output as evidence, not as automatic truth: render
uncertain structures, state warnings, and label unverified mechanism or synthesis
ideas as hypotheses.

Check capabilities only when the task needs optional tools:

```bash
codex-chem doctor
```

## Dependency Posture

Keep the dependency chain small.

- Core: this repository plus RDKit. These are required for normalization,
  validation, drawing, and most reasoning gates.
- Useful optional tools: MolScribe for molecule crops, OSRA as a fallback,
  Open Babel for format conversion, xTB/CREST only when geometry, charge,
  conformers, or rough energetics matter.
- Auxiliary/off-by-default OCSR: ChemSchematicResolver, DECIMER, MolGrapher,
  OpenChemIE, and RxnScribe. Use them only when configured or explicitly helpful
  for a scheme-level question. Do not treat them as required dependencies.
- Figure editors: ChemDraw, ChemDoodle, Marvin, Ketcher, Inkscape, and
  Illustrator are human-review or polishing tools. They are not required for
  ordinary structure validation.

If an optional tool is unavailable, report the limitation and continue with the
best available local evidence. Do not invent missing OCSR, computation, or
literature results.

## Visual Chemistry

1. Convert images to machine-readable chemistry before reasoning.
   - Use `chem_parse_image(path, kind)` for a crop or `chem_parse_scheme(...)`
     for multi-compound schemes.
   - Keep ranked candidates and warnings. For complex schemes, resolve visible
     `R`, `X`, `R1`, `R2`, and abbreviations with graph replacement, not string
     substitution.
   - Multimodal vision may read compound numbers, legends, and abbreviations,
     but final structures must be RDKit-validated SMILES/Molfile/Rxnfile.
   - Run route-consistency checks when adjacent intermediates should preserve a
     scaffold, ring size, stereochemical anchor, or arene substitution pattern.
2. Before synthesis, mechanism, calculation, or final claims, render the selected
   candidate with `chem_input_review` and wait for explicit user confirmation.
3. Never infer stereochemistry from an image unless the parsed structure and
   rendered preview make it explicit.

## Structure And Reaction Flow

- Normalize every molecule with `chem_normalize_structure` before reasoning.
  Check sanitization, dummy atoms, metals, large molecules, and stereo warnings.
- Draw only from validated structures using `chem_draw` or another
  chemistry-aware renderer. SVG/PIL/canvas may be used for page layout, but not
  for hand-drawing molecular bonds.
- For reactions, use `chem_reaction_analyze` for sanity checks or high-level
  condition/retrosynthesis notes. Keep facts, tool output, literature evidence,
  rule-based inference, and LLM assumptions separate.
- For mechanisms, call `chem_input_review` first. After confirmation, use
  `chem_mechanism_draft` as an explanatory hypothesis and `chem_mechanism_render`
  for explicit atom-mapped figure packages.
- Do not pack disconnected reactants, counterions, or byproducts into one
  molecule entry when drawing mechanisms.

## Publication Figures

Use this section for manuscript-ready schemes or mechanisms.

- Prefer ChemDraw ACS-style defaults when available: black structures, restrained
  arrows, standard bond lengths/line widths, readable final lettering, vector
  master output, and no UI/card/slide styling.
- Put prose in the caption or answer, not inside the scheme. Inside the graphic,
  use only structure labels, compound numbers, conditions, short annotations,
  and essential mechanistic labels.
- Curved arrows must be object-bound: start at a visible lone pair or bond and
  end at the receiving atom or bond. Do not use arbitrary SVG curves as chemical
  evidence.
- Show only mechanistically relevant lone pairs by default. Use all lone pairs
  only for teaching/debug figures.
- Treat `publication_checks.blocking_issues` as hard stops. A manuscript-ready
  mechanism needs explicit intermediates, atom-map anchored arrows, charges or
  partial charges where relevant, vector output, and human chemical review.

## Synthesis Reasoning

- Work one decision layer at a time. Use `chem_synthesis_suggest(...,
  confirmed=false)` to render and request confirmation, then
  `confirmed=true` for the first disconnection layer.
- Offer the best route plus realistic alternatives. Explain why alternatives may
  fail: chemoselectivity, redox mismatch, regio/stereocontrol, leaving groups,
  protecting groups, isolation, solubility, or safety.
- Use `chem_literature_search` for unusual or route-defining steps. Cite
  DOI/title/journal/year when available; otherwise say the support is weak.
- Use computation only when it can answer the question. RDKit/xTB/CREST are
  screening evidence, not proof of a mechanism or transition state.
- For wet-lab advice, stay at high-level condition families unless the user has
  provided validated protocols and explicitly asks for operational detail.

## Quick Calls

```bash
uv run codex-chem normalize --smiles "CC(=O)Oc1ccccc1C(=O)O"
uv run codex-chem draw --smiles "c1ccccc1C(=O)O" --output svg
uv run codex-chem parse-image crop.png --kind molecule
uv run codex-chem parse-scheme --image scheme.png --crops ocsr_crops --gold-map legend.json
uv run codex-chem input-review --reaction-smiles "CBr.[OH-]>>CO.[Br-]"
uv run codex-chem reaction-analyze --reaction "CBr.[OH-]>>CO.[Br-]" --mode sanity_check
uv run codex-chem compute --smiles "CCO" --task descriptors --task conformers
uv run codex-chem mechanism-render --spec mechanism.spec.json --output-dir mechanism_out
uv run codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O" --confirmed
```

## Hard Rules

- Do not present raw OCSR as confirmed chemistry.
- Do not present retrosynthesis hints as validated routes.
- Do not call mechanism drafts publication-ready unless the input was confirmed
  and `publication_package.publication_ready` is true.
- Do not let optional tool absence become a hallucination source.
- Use `codex-chem doctor`, `codex-chem figure-tools`, and
  `scripts/install_external_tools_macos.sh` for local capability and setup
  guidance.
