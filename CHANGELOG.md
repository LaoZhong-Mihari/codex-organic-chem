# Changelog

All notable changes to this project are documented here.

## [Unreleased]

### Added

- `mechanism-steps-check` / `chem_mechanism_steps_check` / `kit.mechanism_steps_check`:
  audits stepwise reaction intermediate sequences (plain SMILES lists or
  oMeBench-style step objects) without requiring 2D curved-arrow coordinates,
  checking valence, radical parity, heavy-atom skeleton invariance on proton
  transfers, element provenance, consecutive formula/charge deltas, and target
  closure (flagging `incomplete_cascade_constitutional_isomer` when a sequence
  stops at an un-rearranged constitutional isomer of the product).
- `crop-image` / `chem_crop_image` / `kit.crop_image` and inline bounding-box
  support in `chem_parse_scheme`: renders a 0.1-step normalized coordinate grid
  overlay and 2× Lanczos upscaled sub-crops for dense multi-structure pages.
- MCP tool docstrings across all 28 tools, `Literal` parameter enums,
  `chem://skill/...` MCP resources, `CODEX_CHEM_MCP_COMPACT=1` compact profile,
  and `.github/workflows/ci.yml` for Ubuntu + macOS across Python 3.12+.

### Changed

- `route-check` discrimination benchmark reaches **282/282 (100%)** with zero
  extra false alarms on the clean reference routes: added coordinate-free
  stereospecificity checks for concerted `E`/`Z`-alkene `[2+1]` ring formations
  (epoxidation, cyclopropanation, aziridination), *syn*-dihydroxylation and
  *anti*-dihalogenation; tightened deprotection and `-CH4O`/`-C2H6O` delta rules
  against conditions and substrate groups; prevented missing single-letter IDs
  (`"C"`, `"B"`, `"I"`) from parsing as inline methane/borane/HI; eliminated
  `TBDPS`/`Tr`/`Ph` automorphism exhaustion and moved CIP assignment outside the
  MCS loop (>50× faster on protected polyols); and refined `enone` SMARTS and
  multi-site (2×–6×) polyol protection recognition.
- `model-kit --task face` (`kit.face`) assigns IUPAC 2013 `Re`/`Si` faces via
  topological `rdCIPLabeler` probes, resolving `Re/Si undefined` on achiral
  prochiral carbonyls/alkenes and inverted `Re`/`Si` labels on chiral esters.
- `enumerate --task symmetry` (`kit.symmetry_environments`) resolves
  diastereotopic $^{13}\text{C}$/heavy-atom NMR signals in 3D, and `compare`
  handles OpenSMILES axial allenes (`=[C@]=`/`=[C@@]=`), isotopologues,
  isotopomers, and `tautomers_different_stereo`.
- `mechanism-validate` and `mechanism-render` now support spectator aromatic
  rings, localized Kekulé replay on active aromatic rings ($S_E\text{Ar}$),
  spectator alkali/halide counterions, intermolecular $\pi$-bond attacks
  (`from_bond` to cross-molecule `to_bond` or `to_atom_map`), teaching-mode
  lone-pair/partial-charge rendering, and spec-compliant integer `Element` and
  geometry attributes in CDXML.
- Route figures (`route-figure` / `chem_route_figure`) keep a shared core in
  place (`align_structures`), avoid subscripting precatalyst tags (`XPhos-Pd-G3`),
  measure and wrap wide CJK glyphs accurately, render `row["note"]` headers, and
  preserve duplicate reactants in homodimerizations (`route_to_figure_spec`).
- README section on route and mechanism figures, with the spec of its route
  figure in `docs/examples/hydroxyproline-route-figure.json`.

### Fixed

- Hardened local review server (`structure_review.py`, `workbench.html`,
  `review_daemon.py`): enforced constant-time `X-Review-Token` checks on all
  HTTP `POST` endpoints, blocked path traversal when `integrations/ketcher/dist`
  is absent, restricted CORS and `?api=` to loopback origins, capped request
  body size, replaced `exec_module` in `load_gold_smiles` with AST literal
  parsing, scrubbed secret env vars and killed process groups on timeout in
  `chem_python`, added `flock` daemon startup locking, and preserved unfocused
  drafts in `workbench.html`.
- Fixed stereocentre inversion during `scheme_ocsr.resolve_placeholders`,
  corrected `MEM` (`*COCCOC`) vs `MEMO` (`*OCOCCOC`) attachment and MolScribe
  `@`/`@@` chiral direction, added UFF fallback for boron/silicon in
  `semiempirical.py`, blocked Hessian runs after failed `--opt`, and fixed
  `CODEX_CHEM_REVIEW=off|ask` state transitions.

## [1.1.0] - 2026-11-19

### Added

- Add a verification kit built around checks whose result differs between a
  right and a wrong answer (`codex_organic_chem.kit`, CLI, MCP):
  `route-check` / `chem_route_check` audits a route step by step (formula change
  read against the stated conditions, exact closure with declared reagents and
  byproducts, every stereocentre followed by spatial arrangement under an atom
  mapping, bystander groups versus 33 reagent classes, untraceable atoms, last
  product versus target centre by centre) and can draw every numbered compound;
  `compare` / `chem_compare_structures` classifies two structures and locates
  differing centres; `enumerate` / `chem_enumerate` counts stereoisomers,
  substitution products and symmetry-distinct C/H environments;
  `model-kit` / `chem_model_kit` answers 3D questions from conformer ensembles
  (facial shielding, reach, dihedral populations, planarity, embeddability,
  isomer energies, optional GFN2-xTB refinement).

- `scripts/benchmark_route_check.py`: on 282 single-compound corruptions of three
  curated eight-step routes the audit output changes at blocking or review level
  for 280 (99%); the earlier balance and normalization tools change for 67 (24%)
  and flag all 24 correct steps as unbalanced. Clean routes give no blocking
  issue and 0-3 review prompts. This measures tool discrimination, not model
  accuracy.

- The MCP server sends instructions on connection: what the tools check, and
  to load the skill before the first tool call (`CODEX_CHEM_MCP_INSTRUCTIONS=off`
  leaves them out). In the trigger test they added nothing measurable on top of
  the new description.

- Add an opt-in `chem_python` MCP tool (`CODEX_CHEM_ENABLE_PYTHON=1`) for hosts
  without a shell, and a `chem-route-critic` Claude Code subagent for an
  independent second reading of a route.

- Add MechanismSpec 2.1 and `mechanism-validate` / `chem_mechanism_validate`:
  resolved states, atom/isotope/H/charge conservation, simultaneous electron
  replay, source-electron budgets, and graph-edit checks for localized
  main-group pair transfers and paired homolysis.

- Add five curated mechanism classes, negative controls and exact SVG/JSON
  round-trip checks. Expand and synchronize the mechanism skill contract and
  reusable example. These checks are not a model accuracy comparison.

- Add skeletal/ring drawing cases, geometric regressions and a skill workflow
  grounded in inspected published Schemes and ChemDraw drawing controls.
  Export 2× PNG previews from the vector master. Passing tests and geometry
  checks is explicitly separate from visual review or publication readiness.


### Changed

- Rewrite the skill around a reason / try-to-break-it / fix loop with a
  claim-to-check table. Route planning, reaction-outcome prediction, structure
  input, figures and the kit API move to on-demand references, so the always-
  loaded text is about a third shorter (825 words, was 1285). Route planning now commits to
  one full route and audits it instead of asking for confirmation at each
  disconnection; `chem_synthesis_suggest` remains for interactive brainstorming.

- The skill now gets loaded on exam-style problems. With the old description
  Opus 5.5 loaded it in 3 of 18 sessions on the OrgReason-Real olympiad items
  (as in the benchmark run of 2026-11-17) and 0 of 8 on olympiad-style pages
  written for the test; with the new one 25 of 27 and 6 of 8. Other organic
  questions went from 21/32 to 25/32, and 24 near-miss sessions (pH,
  translating a procedure, an email, code ...) never loaded it (internal
  trigger test). In the benchmark itself (one seed,
  2026-11-18) the skill was loaded in 7 of 9 runs in both D and E, against 1
  and 2 of 9 on 2026-11-17. The description now says that exam and
  olympiad problems, also as page images, are in scope and why to load the
  skill when the chemistry looks familiar or the tools are already listed. It
  stays a plain YAML scalar under 1024 characters (new test).

- Replace the two review pages (a bare Ketcher queue and a SMILES text box)
  with one review workbench served by the package. Structures: source image
  beside the structure, the assistant's reading, question and alternative
  readings, click a stereocentre to flip it, R/S and E/Z buttons for open
  elements, mirror, SMILES/Molfile entry, embedded Ketcher, a live statement of
  how the edit relates to the original, per-item notes, reject and reopen,
  keyboard shortcuts, English and Chinese. Routes (`route-review` /
  `chem_route_review`): the `route-check` audit per step with changed centres
  highlighted, editable structures and conditions with an immediate re-audit,
  a decision and comment per step, an answer per audit prompt, approve or send
  back; the result is the reviewer's version of the route.

- Review pages now outlive the CLI command that opens them (background review
  server, results persisted under `~/.codex-chem/reviews/`); `review-result
  --wait` no longer expires a session when the wait times out. `input-review`
  and `review-batch` open the page without `--wait`; `--static` restores the
  preview-only output. `review-mode auto|ask|off` (or `CODEX_CHEM_REVIEW`)
  decides whether review gates open the page, ask first, or stay silent, and
  the skill now sends image-derived structures and routes to review by default.
  `review-demo` opens two practice sessions with planted mistakes.

- From a first end-to-end run in a fresh session: `input-review --smiles ...
  --source-image ...` reviews the assistant's own reading beside the image;
  `route-review --image` shows the target's source drawing; `chem_parse_image`
  and image-based `chem_input_review` no longer return SVG, Molfile and
  per-atom geometry (88k -> 5k characters for one image; `include_drawings`
  restores them); `chem_route_check` / `chem_route_review` accept
  `route_path`, `chem_compare_structures` accepts `a` / `b`, `chem_parse_image`
  accepts `image_path`; SOCl2 or AcCl in an alcohol is treated as in-situ HCl
  (an esterification that threatens Boc and acetals) instead of as an acid
  chloride preparation that acylates free amines and alcohols. New
  `kit.cis_trans` / `model-kit --task cis_trans`: are the substituents on two
  ring atoms on the same face (both test sessions had written this by hand).

- Literature search returns abstracts, citation counts and open-access links
  from OpenAlex, with Crossref as fallback and supplement.

- Replace the default mechanism painter with the same RDKit skeletal engine
  used by ordinary reaction schemes; keep the old canvas behind an explicit
  `layout.renderer: "legacy"` option. Compact the actual reacting sites,
  preserve unchanged mapped cores, and derive the canvas from its ink bounds.

- Route local cubic electron arrows around atom glyphs, bonds, operators and
  other arrows; block residual collisions and excessive arrow length. Keep
  oxygen lone-pair anchors separate from implicit hydrogen glyphs in `HO−`.

- Keep rounded open V heads for electron pairs. Join the curved body smoothly
  to a short terminal segment along the V's angle bisector, avoiding a shaft
  that crowds one wing on tight arcs. Audit that visible geometry and retain
  single-wing fishhooks for single electrons.

- Remove decorative mechanism-draft arrows; distinguish reaction previews from
  actual mechanism figures. Draw radical dots, chemical subscripts and Kekule
  aromatic bonds; block unsupported stereo and isotope depictions.

- README rewritten as a project introduction: the parts and how they fit,
  features, review-workbench screenshots in English and Chinese, and benchmark
  charts (skill loading, route-check discrimination, OrgReason-Real) with their
  caveats. `README.zh-CN.md` is the Chinese version. The reference sections
  moved unchanged to `docs/` (verification kit, skills and MCP, review
  workbench, figures, external tools); the MCP tool list there now names all 26
  tools. `scripts/make_readme_figures.py` draws the charts.


### Fixed

- Install fixes, each found by installing a clean copy of the repository:
  `mcp` is pinned below 2.0 (a pip install picked mcp 2.x, where `FastMCP` is
  renamed, and the MCP server did not start); `.mcp.json` starts the server with
  `uv run --extra mcp`, so an environment synced without the extra no longer
  breaks it; NumPy and Pillow, imported directly, are now declared.
  `codex-chem doctor` reports the Python side too (interpreter, packages, the
  MCP extra, the packaged review page) with the command that fixes each gap.
  `INSTALL.md` is a step-by-step install with checks and troubleshooting,
  written so that Claude Code or Codex can follow it (`CLAUDE.md` and
  `AGENTS.md` point to it); `scripts/check_install.py` ends it with a real MCP
  handshake. A test checks that every imported package is declared.

- After an E/Z edit on the review page the drawing kept the old
  double-bond geometry, because the edited structure reused every coordinate
  of the original drawing; only the (E)/(Z) label changed. A changed double
  bond is now redrawn (only the smaller side of the bond moves; a ring double
  bond redraws the molecule and turns it onto the old drawing), checked against
  a Molfile round trip on 144 random edits. The queue thumbnail now follows
  edits in progress.

- The N=N of an alkyl azide was reported as an undefined double bond, so
  two azide epimers compared as `same_constitution_stereo_undefined` instead of
  enantiomers or diastereomers.

- Fix same-panel reactant/product mixing in AI traces, preserve explicit mapped
  H, compute trace changes from structures, and embed spec/trace data plus a
  source-spec digest in SVG exports.


### Compatibility

- Legacy 1.0/2.0 specs remain accepted, but chemically inconsistent
  declarations now block publication. The default example uses 2.1. Valid
  electron bookkeeping does not establish kinetics, selectivity or experimental
  mechanism, and editor adapters are not certified lossless round trips.

## [1.0.0] - 2026-09-01

### Added

- Real GFN2-xTB calculation paths for `xtb_opt`, `xtb_reactivity`, and
  `xtb_thermo`, including charge/solvent handling, RDKit atom-index mapping,
  optimized geometries, energies, Fukui indices, frequencies, and free-energy
  corrections.
- Real CREST conformer execution with ensemble energies and populations.
- Full atom and formal-charge balance diagnostics, including implicit
  hydrogens and both catalytic and stoichiometric interpretations of the
  reagent field.
- Canonical multi-engine and multi-image-variant OCSR consensus metadata and
  warnings for uncorroborated single reads.
- A zero-build structure-review editor fallback with live RDKit previews,
  automatic browser opening, and corrected-structure result gating.
- Shared fixed-geometry ACS/RSC/presentation figure styles, SVG/PNG export, and
  publication audits for clipping, overlap, blank output, and bond-length
  drift across molecule, route, and mechanism renderers.
- A Claude Code skill mirror under
  `.claude/skills/organic-chemistry-assistant`, backed by the same `.mcp.json`
  server configuration as Codex and protected by a mirror-consistency test.

### Changed

- Reworked the organic-chemistry skill around structure-grounded reasoning,
  explicit human confirmation, evidence boundaries, and publication checks.
- `input-review --wait` now keeps its localhost editor alive until completion;
  the non-waiting form returns a static preview instead of a dead review URL.

### Fixed

- Hydrogen-only reaction imbalances and formal-charge mismatches are no longer
  reported as balanced.
- Reagent-field atoms and charge are counted when they act as stoichiometric
  partners without forcing catalysts into otherwise balanced equations.
