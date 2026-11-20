---
name: organic-chemistry-assistant
description: Use for organic chemistry whose answer must be right atom by atom - products and selectivity, mechanisms, synthetic routes and retrosynthesis, stereochemistry (R/S, E/Z, cis/trans, meso), counting isomers, stereoisomers or NMR signals, structures deduced from formulas, spectra or images, and molecule, scheme or mechanism figures. Covers exam and olympiad problems (also as page images, in any language) as much as lab and paper work. Load it before you start, even when the chemistry looks familiar and even when chemistry tools are already available, because answers fail on bookkeeping (a flipped stereocentre, a lost carbon, a miscounted isomer, a group that does not survive a reagent) and this skill says which check catches each slip and how to run it with the codex-organic-chem tools (MCP, CLI `codex-chem`, Python `codex_organic_chem.kit`). Not for chemistry without structures (pH, stoichiometry, redox balancing) or for work on the toolchain's own code.
metadata:
  short-description: Structure-grounded organic chemistry with checks that can fail
---

# Organic Chemistry Assistant

You already know the chemistry. What goes wrong is bookkeeping and missing
facts: a flipped wedge, a lost carbon, a miscounted isomer, a group that does
not survive the reagent, a 3D claim never looked at, a catalyst whose known
behaviour was not looked up. The tools exist to catch exactly those. A check
that cannot come back negative (parsing a SMILES you just wrote) proves nothing;
pick checks that would fail if you were wrong.

## The loop

1. **Write structures, not names.** Every species you reason about gets a SMILES
   with all stereochemistry (wedges and E/Z) defined. From an image: follow
   [structure-input](references/structure-input.md) first; it ends at the
   user's review page, also when your reading and OCSR agree.
2. **Reason and commit** to an answer or a route.
3. **Try to break it** with the matching check below, and with your own code:
   `uv run python` with RDKit and `from codex_organic_chem import kit` is your
   model kit and calculator. Never count atoms, isomers or signals, compare
   stereoisomers, or judge a 3D relationship in your head.
4. **Fix what the checks expose**, re-run them, and say in one line what changed
   after checking. Report separately what a tool verified, what literature
   supports, and what is your inference.

| About to claim | Check that can fail |
| --- | --- |
| two structures are the same / which isomer this is | `compare --a --b` (one SMILES: stereo audit) |
| a number: stereoisomers, substitution products, 13C/1H environments | `enumerate --task stereoisomers\|substitution\|environments` |
| a multi-step sequence or route | `route-check` (atoms, every stereocentre, bystander groups, target identity) |
| cis or trans on a ring, which face is open, can two atoms reach, anti-periplanar, flat or not, can this ring fusion exist, which isomer is more stable | `model-kit --task cis_trans\|face\|distance\|dihedral\|planarity\|embed\|isomer_energies` |
| which site is most nucleophilic/electrophilic; conformer populations; free energies | `compute --task xtb_reactivity\|crest\|xtb_thermo` |
| a written equation is balanced | `reaction-analyze` (`mass_charge_balance`) |
| electron-pushing steps are consistent | `mechanism-validate` |
| this reaction is known / this catalyst gives this selectivity / this reagent is available | `literature-search` and web search when allowed |

MCP names: `chem_compare_structures`, `chem_enumerate`, `chem_route_check`,
`chem_model_kit`, `chem_compute`, `chem_reaction_analyze`,
`chem_mechanism_validate`, `chem_literature_search`. API and the route format:
[kit-api](references/kit-api.md).

Two things no check can tell you: what the user's drawing really shows, and
what they already know about a step. For those the user has a review page.
Structures you read from an image, and a route they will act on, go to it
(`review-batch` / `input-review` / `route-review`) as a matter of course; the
command answers `ask_user_first` or `review_disabled` when their setting says
otherwise. See [human-review](references/human-review.md).

## Task protocols (read the one you need)

- **Plan or audit a synthesis** -> [route-planning](references/route-planning.md).
  Commit to one full route, write it as data, run `route-check`, fix, then
  present with the three riskiest steps and tactical fallbacks.
- **Predict a product, selectivity or best conditions** ->
  [reaction-outcome](references/reaction-outcome.md). Enumerate competing
  pathways before choosing; catalyst- and reagent-specific behaviour is decided
  by precedent, not by structure checks.
- **Structures from images, OCSR** ->
  [structure-input](references/structure-input.md).
- **A person checks or corrects structures or a route** ->
  [human-review](references/human-review.md): when to open the review page,
  what to send, how to wait, how to use what comes back.
- **Mechanisms with curved arrows** -> [mechanism-spec](references/mechanism-spec.md)
  (write a MechanismSpec 2.1, `mechanism-validate`, then render).
- **Any figure for a human** -> [figures](references/figures.md) and, for
  mechanisms, [mechanism-graphics](references/mechanism-graphics.md).

## What a passing check means

- `route-check` `blocking` issues are facts (structure unparseable, declared
  balance does not close, last product is not the target). `review` items come
  from heuristic atom mapping and reagent rules: answer each in one line or fix
  it; some are false alarms. `consistent` means the bookkeeping holds, not that
  a step works.
- `model-kit` is sterics and ground-state geometry from conformers. It finds no
  transition state; torsional, electronic and chelation control must be argued.
- xTB/CREST numbers are approximate: gaps under ~1.5-2 kcal/mol decide nothing.
- `mechanism-validate` `semantic_checks.status=valid` proves the supported
  electron bookkeeping only; `incomplete` or `unsupported` cannot back a claim
  that a mechanism was validated, and no status establishes kinetics.
- A legal, balanced, stereo-defined product is a well-formed candidate, not the
  observed product. Selectivity claims need precedent or a stated assumption.
- Literature hits are abstract-level. Cite DOI, grade as exact precedent >
  close analogue > same class > review, and never invent a reference.
- Wet-lab advice stays at condition-family level unless the user supplies a
  validated protocol and asks for operational detail.
- If a tool is missing (`codex-chem doctor`), say so; do not invent its output.

## Quick calls

```bash
uv run codex-chem compare --a "C[C@H](O)CC" --b "C[C@@H](O)CC"
uv run codex-chem enumerate --smiles "CC1CCCCC1" --task substitution --substituent Cl
uv run codex-chem route-check --example > route.json
uv run codex-chem route-check --route route.json --figure-dir out --issues-only
uv run codex-chem route-review --route route.json    # the user checks it in the browser
uv run codex-chem review-result --session-id ID --wait
uv run codex-chem model-kit --task face --smiles "O=C1C[C@@H]2CC[C@H]1C2" --atoms 1
uv run codex-chem model-kit --task isomer_energies --smiles "C/C=C/C" --smiles "C/C=C\\C" --method xtb
uv run codex-chem compute --smiles "CC(=O)C" --task xtb_reactivity --solvent water
uv run codex-chem literature-search --query "cobalt alkyne migratory hydroboration" --rows 5
uv run codex-chem mechanism-validate --spec mechanism.spec.json
uv run python -c "from codex_organic_chem import kit; print(kit.compare('CCO','COC')['relation'])"
```
