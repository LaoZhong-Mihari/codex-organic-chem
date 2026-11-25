---
name: organic-chemistry-assistant
description: Use for organic chemistry whose answer must be right atom by atom - products and selectivity, stepwise reaction mechanisms and intermediate SMILES lists (including strict JSON or persona prompts), synthetic routes and retrosynthesis, stereochemistry (R/S, E/Z, cis/trans, axial, meso), counting isomers, stereoisomers or NMR signals, structures deduced from formulas, spectra or dense page images, and molecule, scheme or mechanism figures. Covers exam and olympiad problems (also as page images, in any language) as much as lab and paper work. Load it before you start, even when the chemistry looks familiar or tools are already available, because answers fail on bookkeeping (a flipped stereocentre, a lost carbon, a miscounted isomer, an unfinished cascade, a group that does not survive a reagent) and this skill says which check catches each slip with codex-organic-chem (MCP, CLI `codex-chem`, Python `codex_organic_chem.kit`). Not for chemistry without structures (pH, redox balancing) or toolchain code edits.
metadata:
  short-description: Structure-grounded organic chemistry with checks that can fail
---

# Organic Chemistry Assistant

You already know the chemistry. What goes wrong is bookkeeping and missing
facts: a flipped wedge, a lost carbon, a miscounted isomer, an unfinished
cascade, a group that does not survive the reagent, a 3D claim never looked at,
a catalyst whose known behaviour was not looked up. The tools exist to catch
exactly those. A check that cannot come back negative (parsing a single SMILES
you just wrote) proves nothing; pick checks that would fail if you were wrong.

## The loop

1. **Write structures, not names.** Every species you reason about gets a SMILES
   with all stereochemistry (wedges, E/Z, axial) defined. On dense multi-structure
   page images, first zoom/crop each structure with `chem_crop_image` (`crop-image`)
   and follow [structure-input](references/structure-input.md); interactive runs
   end at the user's review page unless review mode is `off`.
2. **Reason and commit** to an answer, intermediate list, or route.
3. **Try to break it** with the matching check below, and with your own code
   when available (`chem_python` or `uv run python` with RDKit and
   `from codex_organic_chem import kit`). Never count atoms, isomers or NMR
   signals, compare stereoisomers, audit a mechanism list, or judge a 3D
   relationship in your head.
4. **Fix what the checks expose**, re-run them, and say in one line what changed
   after checking (unless the user requires strict machine-only JSON output, in
   which case fix silently and emit the exact requested JSON). Report separately
   what a tool verified, what literature supports, and what is your inference.

| About to claim | Check that can fail |
| --- | --- |
| two structures are the same / which isomer this is | `compare --a --b` / `chem_compare_structures` (one SMILES: stereo audit) |
| a number: stereoisomers, substitution products, 13C/1H environments | `enumerate --task stereoisomers\|substitution\|environments` / `chem_enumerate` |
| a stepwise reaction mechanism (list of intermediates) | `mechanism-steps-check` / `chem_mechanism_steps_check` (valence, radical parity, proton-transfer skeleton, step deltas, target closure) |
| a multi-step synthetic sequence or route | `route-check` / `chem_route_check` (atoms, every stereocentre, bystander groups, target identity) |
| cis or trans on a ring, which face (Re/Si) is open, can two atoms reach, anti-periplanar, flat or not, can this ring fusion exist, which isomer is more stable | `model-kit --task cis_trans\|face\|distance\|dihedral\|planarity\|embed\|isomer_energies` / `chem_model_kit` |
| which site is most nucleophilic/electrophilic; conformer populations; free energies | `compute --task xtb_reactivity\|crest\|xtb_thermo` / `chem_compute` |
| a single written equation is balanced | `reaction-analyze --mode sanity_check` / `chem_reaction_analyze` |
| curved-arrow figure spec (MechanismSpec 2.1) is consistent | `mechanism-validate` / `chem_mechanism_validate` |
| this reaction is known / this catalyst gives this selectivity / this reagent is available | `literature-search` / `chem_literature_search` and web search when allowed |

API and formats: [kit-api](references/kit-api.md) (also exposed as MCP resources `chem://skill/references/...`).

Two things no check can tell you: what the user's drawing really shows, and
what they already know about a step. Structures you read from an image, and a
route they will act on, go to the review page (`chem_structure_review_batch` /
`chem_input_review` / `chem_route_review`) when a person is present; skip it in
batch/headless runs or when the tool answers `review_disabled`. See
[human-review](references/human-review.md).

## Task protocols (read the one you need)

- **Derive or check a reaction mechanism (list of intermediates):** Write every
  elementary step (`proton_transfer`, `nucleophilic_addition`,
  `leaving_group_elimination`, `1,2_shift`, etc.) with explicit charges and call
  `chem_mechanism_steps_check(steps=..., reactants=..., products=..., conditions=...)`.
  - Do **not** collapse protonation + C–X bond cleavage, or nucleophilic attack +
    proton transfer, into one step unless genuinely concerted.
  - **Cascade / isomer rule:** When a problem supplies a valid target product and
    `chem_mechanism_steps_check` / `chem_compare_structures` reports that your
    adduct is a `constitutional_isomer` or `tautomer` of the target, **do not
    claim the target SMILES is a typo**: look for a subsequent ring-opening /
    reclosure (retro-aldol, retro-Claisen/Dieckmann, retro-Mannich), 1,2-shift,
    Smiles/Dimroth rearrangement, or tautomerization that reaches the target.
- **Plan or audit a synthesis** -> [route-planning](references/route-planning.md).
  Commit to one full route, write it as data, run `route-check`, fix, then
  present with the three riskiest steps and tactical fallbacks.
- **Predict a product, selectivity or best conditions** ->
  [reaction-outcome](references/reaction-outcome.md). Enumerate competing
  pathways before choosing; catalyst- and reagent-specific behaviour is decided
  by precedent, not by structure checks.
- **Structures from images, dense exam pages, OCSR** ->
  [structure-input](references/structure-input.md). Crop/zoom individual
  molecules first (`chem_crop_image`) before reading wedges or running OCSR.
- **A person checks or corrects structures or a route** ->
  [human-review](references/human-review.md).
- **Mechanisms with curved arrows (artwork)** -> [mechanism-spec](references/mechanism-spec.md)
  (write a MechanismSpec 2.1, `mechanism-validate`, then render).
- **Any figure for a human** -> [figures](references/figures.md) and
  [mechanism-graphics](references/mechanism-graphics.md).

## What a passing check means

- `route-check` and `mechanism-steps-check` `blocking` issues are hard conflicts
  (unparseable structure, skeleton changed during proton transfer, unclosed
  equation, last species is not the target). `review` items come from heuristic
  mapping/delta rules: resolve or justify each. `consistent` proves chemical
  bookkeeping, not kinetic feasibility or regioselectivity.
- `model-kit` is sterics and ground-state geometry from conformers. It finds no
  transition state; torsional, electronic and chelation control must be argued.
- xTB/CREST numbers are approximate: gaps under ~1.5-2 kcal/mol decide nothing.
- `mechanism-validate` `semantic_checks.status=valid` proves supported electron
  bookkeeping only; `incomplete` or `unsupported` cannot back a validation claim.
- Literature hits are abstract-level. Cite DOI, grade as exact precedent >
  close analogue > same class > review, and never invent a reference.
- If a tool is missing (`codex-chem doctor`), say so; do not invent its output.

## Quick calls

```bash
uv run codex-chem compare --a "C[C@H](O)CC" --b "C[C@@H](O)CC"
uv run codex-chem enumerate --smiles "CC1CCCCC1" --task substitution --substituent Cl
uv run codex-chem mechanism-steps-check --reactants "C=CCC1C(=O)CCCC1=O.CC(=O)C=C" --products "CC(=O)CCC1(CC=C)C(=O)CCCC1=O" --steps '["C=CCC1=C([O-])CCCC1=O","C=CCC1(CC=C(C)[O-])C(=O)CCCC1=O","CC(=O)CCC1(CC=C)C(=O)CCCC1=O"]'
uv run codex-chem crop-image exam_page.jpg --grid --output-dir crops/
uv run codex-chem route-check --route route.json --figure-dir out --issues-only
uv run codex-chem model-kit --task face --smiles "CCC(=O)C" --atoms 2
uv run codex-chem compute --smiles "CC(=O)C" --task xtb_reactivity --solvent water
uv run codex-chem literature-search --query "cobalt alkyne migratory hydroboration" --rows 5
```
