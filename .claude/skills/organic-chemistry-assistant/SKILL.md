---
name: organic-chemistry-assistant
description: Use when validating, drawing, computing, parsing, or explaining organic molecules, reaction schemes, mechanisms, retrosynthesis, or literature-supported route screening with the local codex-organic-chem toolchain (CLI `codex-chem` or MCP `codex-organic-chem`).
metadata:
  short-description: Structure-grounded organic chemistry reasoning and publication figures
---

# Organic Chemistry Assistant

Use the `codex-chem` CLI or `codex-organic-chem` MCP tools before making
nontrivial chemistry claims. Tool output is evidence, not truth: keep warnings
visible, label unverified mechanism or synthesis ideas as hypotheses, and never
let a missing optional tool become a reason to invent results. Check optional
capabilities only when needed: `codex-chem doctor`.

## Reason from structure, not from names

- Normalize every molecule with `chem_normalize_structure` before reasoning;
  check sanitization, dummy atoms, metals, and stereo warnings.
- Inspect the rendered molecule, then build a structure map: functional groups,
  heteroatoms, pi systems, acidic/basic sites, nucleophilic/electrophilic
  sites, leaving/protecting groups, steric congestion, stereo anchors, and
  groups that must stay untouched. Reason from mapped reactive sites — never
  from a compound name, formula, or remembered reaction pattern without
  checking the required atoms and competing sites exist.
- `chem_reaction_analyze` for sanity checks and condition/retro notes. It closes
  mass/charge balance (hydrogens included): `analysis.mass_charge_balance`
  `unbalanced` means the equation as written is wrong, not just
  under-annotated — fix it or say what is missing before reasoning on it. Keep
  tool facts, literature evidence, rule inference, and assumptions separate.
- `chem_compute` when data can change the decision: `xtb_reactivity` ranks
  nucleophilic/electrophilic sites by Fukui index and GFN2 charge per RDKit atom
  index; `xtb_opt`/`xtb_thermo`/`crest` add geometry, free energy and imaginary
  modes, conformer populations. Ions and `--solvent` are handled; read warnings.

## Mechanistic reasoning

1. Establish substrate and reagent roles from validated structures.
2. Mark likely reactive atoms/bonds and the groups controlling chemo-, regio-,
   and stereoselectivity or acid/base state.
3. Enumerate plausible elementary steps before choosing one.
4. Balance each step with `chem_reaction_analyze`, then check what it cannot:
   conserved scaffolds, stereochemistry, which atoms need mapping.
5. State why the chosen path beats close alternatives, especially when another
   functional group could react under the same conditions.

## Images in: OCSR

- `chem_parse_image` (crops) / `chem_parse_scheme` (multi-compound schemes).
  Preprocessing variants and all configured engines run automatically; ranking
  rewards independent engines agreeing on one canonical structure.
- Read the top candidate's `metadata.consensus`: `multi_engine` is strong;
  `single_read` means compare the rendered candidate against the source image
  before trusting it. Resolve `R`/`X` groups and abbreviations by graph
  replacement, never string substitution.
- Raw OCSR is never confirmed chemistry. Confirmation is **interactive, in the
  editor — not "look at this picture"**: `chem_input_review` (single input or
  OCSR result) and `chem_structure_review_batch` (multiple SMILES) open the
  local structure editor in the browser automatically (Ketcher when built, a
  built-in fallback otherwise). Tell the user it is open (give `review_url` if
  `browser_opened` is false), block on
  `chem_structure_review_result(session_id=..., wait=true)`, and continue only
  with the returned structures — the user may have *corrected* them, so never
  reuse the pre-review SMILES. Never infer stereochemistry an image does not
  make explicit.
- From the CLI, use `codex-chem input-review ... --wait`; without `--wait` the
  command intentionally returns a static preview because its daemon review
  server cannot outlive the command. MCP sessions remain live in the MCP
  server process.

## Images out: figures

All renderers share one fixed-geometry engine (uniform bond length, ACS-style
defaults) and audit their own output — overlap, clipping, bond-length drift,
blank exports. Results carry `publication_checks`; blocking issues downgrade
`status` to `blocked_for_publication`.

- Molecules/reactions: `chem_draw` (SVG master; PNG derived after the audit via
  `output_file`).
- Routes/schemes: `chem_route_figure` (`style_preset`: `acs` default, `rsc`,
  `presentation`). Use it for any rendered scheme; RDKit one-call previews are
  for internal checks only.
- Mechanisms: `chem_input_review` first; after confirmation `chem_mechanism_draft`
  (hypothesis) and `chem_mechanism_render` (atom-mapped panels, curved arrows
  anchored to lone pairs/bonds/atoms, CDXML/ChemDoodle/Ketcher exports).
- Chemistry-aware renderers draw all bonds and arrows; hand-built SVG/PIL/canvas
  may only do page composition. Prose belongs in the caption, not the graphic.
- Report artifact status honestly: confirmed, draft, or blocked — with the
  audit's warnings. `blocked_for_publication` output must be fixed or reported
  as blocked, never shipped. A figure is "publication-ready" only when the
  input was confirmed, `publication_checks` passes, and a human reviewed it.

## Synthesis reasoning

- Work one decision layer at a time: `chem_synthesis_suggest(...,
  confirmed=false)` renders and requests confirmation; `confirmed=true`
  proceeds with the first disconnection layer.
- Screen broad, then narrow by structural fit, functional-group tolerance,
  step economy, selectivity risk, redox/protecting-group burden, precedent,
  and starting-material availability. Explain why alternatives fail.
- `chem_literature_search` for route-defining or disputed steps: exact
  transformation first, then analogous class/catalyst/substrate. Cite
  DOI/title/year, and grade support: exact precedent > close analog > same
  class > review > unsupported inference. Weak support must be called weak.
- Wet-lab advice stays at condition-family level unless the user supplies
  validated protocols and asks for operational detail.

## Quick calls

```bash
uv run codex-chem normalize --smiles "CC(=O)Oc1ccccc1C(=O)O"
uv run codex-chem draw --smiles "c1ccccc1C(=O)O" --output svg
uv run codex-chem parse-image crop.png --kind molecule
uv run codex-chem input-review --reaction-smiles "CBr.[OH-]>>CO.[Br-]" --wait
uv run codex-chem review-batch --input reviews.json --wait
uv run codex-chem route-figure --spec route.spec.json --output-dir out --format svg --format png
uv run codex-chem reaction-analyze --reaction "CBr.[OH-]>>CO.[Br-]" --mode sanity_check
uv run codex-chem compute --smiles "CC(=O)C" --task xtb_reactivity --solvent water
uv run codex-chem mechanism-render --spec mechanism.spec.json --output-dir out
uv run codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O" --confirmed
```
