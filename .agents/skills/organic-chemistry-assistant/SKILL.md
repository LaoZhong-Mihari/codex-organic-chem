---
name: organic-chemistry-assistant
description: Use this skill when Codex needs to understand, validate, draw, compute, or explain organic molecules, reaction schemes, retrosynthesis hints, or reaction mechanisms using the local codex-organic-chem CLI/MCP toolchain. Trigger for molecule images, SMILES/Molfile cleanup, reaction plausibility checks, lightweight RDKit/xTB calculations, and mechanism drafts.
metadata:
  short-description: Local organic chemistry reasoning tools for Codex
---

# Organic Chemistry Assistant

Use the local `codex-chem` CLI or `codex-organic-chem` MCP tools before making
nontrivial chemistry claims. Prefer tool facts over memory, and mark unverified
mechanistic or synthetic ideas as hypotheses.

First check local capability when a task depends on optional tools:

```bash
codex-chem doctor
```

## Core workflow

1. Convert visual chemistry into machine-readable structure first.
   - Use `chem_parse_image(path, kind)` for screenshots, scans, schemes, or hand drawings.
   - Immediately render the candidate with `chem_input_review`; stop until the user confirms the rendered molecule/reaction or provides corrected SMILES/Molfile/Rxnfile.
2. Normalize every molecule before reasoning.
   - Use `chem_normalize_structure(smiles=...)` or Molfile input.
   - Check warnings for stereochemistry, metals, wildcard atoms, large structures, and sanitization issues.
3. Draw or compute only from validated structures.
   - Use `chem_draw` for SVG/Molfile output.
   - Use `chem_compute` for descriptors, conformers, charges, and optional xTB/CREST status.
4. For reactions, separate facts from inference.
   - Use `chem_reaction_analyze(input, mode)` for sanity checks, forward notes, condition families, or retrosynthesis hints.
   - For manuscript-style mechanisms, first call `chem_input_review(reaction_smiles=...)`; after user confirmation, call `chem_mechanism_draft(reaction, style, quality="publication", structure_confirmed=true)` only as the hypothesis/explanation layer.
   - For actual figure output, use `chem_mechanism_render(spec)` with explicit intermediates, atom-map anchored arrows, mechanistically relevant lone pairs, formal charges, and partial charges. Do not call a template mechanism "publication-ready" unless the explicit canvas is provided and chemically reviewed.
   - For publication-style mechanism figures, be visually conservative: show only the lone pair(s) that are the actual electron source for an arrow unless the user asks for all electrons; omit product/intermediate lone-pair dots unless they explain the step; use restrained black or single-color curved arrows; keep curves short, smooth, and clear of atom labels; avoid decorative atom colors unless they carry meaning.
   - Treat `publication_checks.blocking_issues` as hard stops. Fix missing intermediates, missing atom maps, hidden lone pairs, or arrows that start from atom centers before presenting a mechanism as manuscript figure material.
   - Use `chem_figure_tool_status()` or `codex-chem figure-tools` before recommending ChemDraw/ChemDoodle/Marvin/Ketcher follow-up. With `chem_mechanism_render(..., output_dir=...)`, inspect `editor_followup` for local open-command templates. ChemDraw support is file/workflow integration, not unattended GUI automation.
   - Do not pack disconnected reactants, counterions, or byproducts into one molecule entry when drawing mechanisms. Use separate molecule entries so layout and electron-flow arrows stay readable.
   - Explain every elementary step: electron movement, bond changes, charge/proton changes, why the step is plausible, and what must be checked.
5. For synthesis suggestions, work one decision layer at a time.
   - Use `chem_synthesis_suggest(target_smiles, confirmed=false)` first to render and require confirmation.
   - After confirmation, use `chem_synthesis_suggest(..., confirmed=true)` to propose only the first disconnection layer.
   - Give several plausible disconnection choices at the current decision layer, rank the most feasible route, and keep alternatives visible. Do not expand a full operational multistep protocol until the user chooses/rejects the previous disconnection.
   - Include tool-derived data and literature metadata separately; use `chem_literature_search` when additional precedent evidence is needed.
6. Always state boundaries.
   - Tool facts: RDKit parsing, descriptors, drawings, atom-map deltas.
   - Literature facts: DOI/title/journal/year/search query from available literature search tools.
   - Rule inferences: functional-group, reaction-class, and disconnection heuristics.
   - LLM assumptions: any explanation not directly produced by a tool.

## ACS manuscript figure standard

Use this section whenever the user asks for a publishable mechanism, reaction
scheme, structure figure, ACS-style figure, or manuscript-ready chemistry art.
The goal is a compact chemical scheme, not an explanatory slide.

ACS/ChemDraw structure defaults to emulate:

- Use ChemDraw's ACS/ACS-1996 document settings when ChemDraw is available.
- Chain angle: 120 degrees.
- Bond spacing: 18% of width.
- Fixed bond length: 14.4 pt / 0.508 cm / 0.2 in.
- Bold width: 2.0 pt.
- Line width: 0.6 pt.
- Margin width: 1.6 pt.
- Hash spacing: 2.5 pt.
- Text font: Arial or Helvetica, 10 pt for chemical labels in ChemDraw.
- Preferences: units in points, tolerance 5 pixels, US Letter page, 100% scale.
- Scheme/structure block width: target final production size. Prefer one-column
  width when possible; ACS examples cite 240 pt / 3.33 in for many journal
  graphics and 11.3 cm for structure blocks in some author guidelines. Use
  double-column only when the chemistry cannot be read at one-column width;
  keep within 504 pt / 7 in, or 23.6 cm where the journal explicitly uses the
  structure-block rule.
- Maximum depth for normal graphics: 660 pt / 9.167 in including caption space.
- Minimum final lettering: never below 4.5 pt; prefer 6-8 pt or larger after
  reduction. Lines should not be thinner than 0.5 pt.
- File outputs: keep a vector master (`.cdxml`/`.cdx`, `.svg`, `.pdf`, or
  `.eps`). If a raster file is required, export at final size: 1200 dpi for
  black-and-white line art, 600 dpi for grayscale, 300 dpi for RGB color.
  Embed or outline fonts in EPS/PDF.

ACS-style mechanism/scheme design rules:

- Do not make final manuscript figures with PIL/PNG card layouts, rounded panels,
  large titles, explanatory paragraphs, colored callouts, or screenshots. Those
  are rough discussion sketches only. The file
  `/Users/wzzzz11/Documents/Codex/2026-04-28/files-mentioned-by-the-user-122668/render_stepwise_mechanism_detail.py`
  is an example of a non-ACS sketch style: useful for brainstorming, not for
  submission.
- Put prose in the manuscript text or caption, not inside the scheme. Inside the
  graphic, use only structure labels, compound numbers, conditions, short arrow
  annotations, and essential mechanistic labels.
- Use black structures and black reaction/mechanism arrows by default. Use color
  only sparingly to identify a reaction center or compare pathways; the figure
  must remain readable in grayscale.
- Avoid card borders, UI panels, background boxes, shadows, gradients, and
  decorative colors. If multiple steps are necessary, use a clean multi-panel
  scheme with panel letters or step numbers, aligned structures, and even
  whitespace.
- Show mechanistically relevant electron flow only. Use source-only lone pairs
  for nucleophiles (`layout.lone_pair_display="source_only"` and
  `visible_count: 1`); omit product/intermediate lone pairs unless they are
  needed to understand that elementary step.
- Curved arrows must be object-bound: start at a visible lone pair or bond and
  end at the receiving atom or bond. Keep arrows smooth, short, and clear of
  atom labels; do not use freehand arcs as chemical evidence.
- Use boldface only for compound numbers, not atom labels or prose captions.
- A final mechanism must have a vector master plus a chemistry audit trail:
  mapped structures, explicit intermediates, atom-map anchored arrows,
  `mechanism.trace.json`, and human chemical review.

Tool workflow for ACS-ready chemistry art:

1. Use Ketcher or ChemDraw to correct structures/reactions first. Ketcher is
   mainly the open-source input correction gate; do not rely on Ketcher SVG
   alone for final electron-pushing mechanisms.
2. Use `chem_mechanism_render(spec, output_dir=...)` to generate the checked
   `mechanism.svg`, `mechanism.cdxml`, `mechanism.trace.json`,
   `mechanism.chemdoodle.json`, `mechanism.ketcher.json`, and
   `mechanism.marvin.json` package. In the spec set
   `journal_style="acs"`, `presentation_mode="publication"`, hide titles/notes,
   and use ACS final-size dimensions rather than a large poster canvas.
3. Open `mechanism.cdxml` in ChemDraw when available and apply the ACS
   stylesheet/document settings. Use ChemDraw for final arrow geometry,
   typography, compound numbering, and journal-style polishing.
4. If ChemDraw is not available, use ChemDoodle or Marvin for object-bound
   mechanism review, then Inkscape/Illustrator only for final vector cleanup.
   Vector polish tools must not change the chemistry.
5. Export the vector master and, only if needed, a final-size high-resolution
   TIFF/PNG according to the ACS resolution rules above.
6. Before calling a figure manuscript-ready, verify: no UI/card layout, no
   paragraphs in the art, ACS dimensions, legible final font size, line width
   above minimum, vector master saved, atom maps/charges/arrows reviewed, and
   all speculative steps labeled as proposed in the caption or surrounding text.

## Synthesis reasoning standard

Use this section whenever proposing, comparing, validating, or explaining
organic synthesis routes.

Route exploration requirements:

- Generate multiple plausible retrosynthetic disconnections or forward routes
  when chemistry permits. Rank them as primary route, backup route(s), and
  lower-confidence/speculative route(s).
- For each route, separate what is known from what is inferred: direct
  literature precedent, close-analog precedent, local tool/computation result,
  first-principles organic chemistry reasoning, and LLM hypothesis.
- Prefer the route with the best combined balance of step economy, chemoselective
  control, functional-group tolerance, reagent availability, safety,
  purification realism, stereochemical/regiochemical control, and evidence.
- Do not hide alternatives. Include why each alternative may fail: competing
  reactions, overreaction, poor leaving group, incompatible acid/base conditions,
  redox mismatch, steric congestion, tautomer/protonation ambiguity,
  regioselectivity, stereoselectivity, solubility, or isolation risk.

Evidence and validation requirements:

- Use local tools before making strong claims: normalize structures with
  `chem_normalize_structure`, sanity-check reactions with
  `chem_reaction_analyze`, and use `chem_compute` for descriptors/conformers,
  charges, xTB, or CREST when the question depends on geometry, strain,
  conformers, charge distribution, or relative thermodynamics.
- Use `chem_literature_search` for named transformations, unusual steps,
  proposed key bond formations, protecting-group choices, rearrangements, and
  steps that determine the route's credibility. Cite DOI/title/journal/year or
  say explicitly that no close precedent was found.
- When computation is feasible, report the method and limits. Treat RDKit/xTB
  results as screening evidence, not proof. Do not overclaim transition states
  unless an actual transition-state calculation/literature precedent supports it.
- When literature is unavailable or weak, reason from first principles and label
  the conclusion as an inference.

First-principles analysis checklist for every important step:

- Electron flow: nucleophile/electrophile identity, orbital interaction, leaving
  group, bond made/broken, and charge/proton bookkeeping.
- Acid/base: approximate pKa logic, protonation state under conditions, buffer or
  counterion role, and whether proton transfer is intra- or intermolecular.
- Selectivity: chemo-, regio-, stereo-, and site-selectivity; steric and
  conformational effects; directing groups; hard/soft matching.
- Thermodynamics/kinetics: driving force, reversibility, strain relief/formation,
  aromaticity, precipitation/gas evolution, redox balance, and likely
  rate-determining issue.
- Practicality: solvent/reagent compatibility, moisture/air sensitivity,
  temperature sensitivity, purification/isolation plausibility, scale and safety
  concerns at a high level.

Output format for synthesis answers:

- Start with the recommended route and a confidence level. Then list alternative
  routes with pros/cons and what evidence would decide between them.
- For each proposed step, include: transformation, rationale, evidence level,
  key risks, validation needed, and whether the step is literature-backed,
  computation-supported, or hypothesis-only.
- Include mechanism or route figures whenever useful. Use the ACS manuscript
  figure standard above and `chem_mechanism_render` / `chem_draw` outputs rather
  than explanatory PNG card layouts. Prefer vector files plus a rendered preview.
- Keep the explanation teachable like a textbook: show the reaction logic and
  electron flow, but keep the final graphic clean and publishable.

## CLI quick calls

```bash
uv run codex-chem normalize --smiles "CC(=O)Oc1ccccc1C(=O)O"
uv run codex-chem draw --smiles "c1ccccc1C(=O)O" --output svg
uv run codex-chem compute --smiles "CCO" --task descriptors --task conformers
uv run codex-chem reaction-analyze --reaction "CBr.[OH-]>>CO.[Br-]" --mode sanity_check
uv run codex-chem input-review --reaction-smiles "CBr.[OH-]>>CO.[Br-]"
uv run codex-chem mechanism-draft --reaction "CBr.[OH-]>>CO.[Br-]" --quality publication --confirmed
uv run codex-chem mechanism-render --example
uv run codex-chem mechanism-render --spec mechanism.spec.json --output-dir mechanism_out
uv run codex-chem figure-tools
uv run codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O" --confirmed
```

## Safety and quality rules

- Do not present retrosynthesis hints as validated routes.
- Do not present mechanism drafts as publication-ready unless the input was rendered/confirmed and `publication_package.publication_ready` is true.
- If the user wants publication-level mechanism art, require explicit intermediates and atom-map anchors. Literature-quality mechanism figures need visible electron-source lone pairs or bond anchors, charges/partial charges only where they clarify the mechanism, and arrows that start at electron sources and end at electron sinks.
- In `presentation_mode="publication"`, prefer `layout.lone_pair_display="source_only"` and `visible_count: 1` for nucleophile lone-pair sources. Use `lone_pair_display="all"` only for teaching/debug figures where showing every heteroatom lone pair is the point.
- Literature-style mechanism drawings are usually made in ChemDraw/CDXML, ChemDoodle, or ChemAxon Marvin, then polished in ChemDraw, Illustrator, or Inkscape. RDKit-only mechanism overlays are schematic and must be labeled as such.
- Prefer Ketcher as the open-source correction/editor option for molecule or reaction input review. Treat it as an input correction tool unless the local deployment supports all needed mechanism annotations.
- Follow the ChemDoodle/Marvin design pattern: electron-pushing arrows attach to chemical objects (atom, bond, lone pair) and are then laid out; arbitrary SVG curves are not sufficient evidence of a valid mechanism.
- Do not continue after OCSR/image recognition until the user has confirmed the rendered structure.
- Do not infer stereochemistry from images unless the parsed structure and preview make it explicit.
- For wet-lab advice, include safety/literature-review caveats and avoid operational instructions beyond high-level condition families unless the user explicitly provides validated protocols.
- If Open Babel, xTB, CREST, OSRA, MolScribe, or RxnScribe are unavailable, report that as a tool limitation instead of filling the gap from memory.
- For installation guidance, use `codex-chem doctor` or `scripts/install_external_tools_macos.sh` in the project. CREST may live in `~/.local/share/codex-organic-chem/conda-tools/bin` and is discovered automatically.
- `codex-chem doctor` includes `publication_figure_guide`; `codex-chem figure-tools` includes local availability and open-command templates for ChemDraw/ChemDoodle/Marvin/Ketcher/Inkscape/Illustrator.
