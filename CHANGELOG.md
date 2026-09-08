# Changelog

All notable changes to this project are documented here.

## [Unreleased]

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
- Add skeletal/ring drawing cases, geometric regressions and a skill workflow
  grounded in inspected published Schemes and ChemDraw drawing controls.
  Export 2× PNG previews from the vector master. Passing tests and geometry
  checks is explicitly separate from visual review or publication readiness.
- Add MechanismSpec 2.1 and `mechanism-validate` / `chem_mechanism_validate`:
  resolved states, atom/isotope/H/charge conservation, simultaneous electron
  replay, source-electron budgets, and graph-edit checks for localized
  main-group pair transfers and paired homolysis.
- Fix same-panel reactant/product mixing in AI traces, preserve explicit mapped
  H, compute trace changes from structures, and embed spec/trace data plus a
  source-spec digest in SVG exports.
- Remove decorative mechanism-draft arrows; distinguish reaction previews from
  actual mechanism figures. Draw radical dots, chemical subscripts and Kekule
  aromatic bonds; block unsupported stereo and isotope depictions.
- Add five curated mechanism classes, negative controls and exact SVG/JSON
  round-trip checks. Expand and synchronize the mechanism skill contract and
  reusable example. These checks are not a model accuracy comparison.
- Compatibility: legacy 1.0/2.0 specs remain accepted, but chemically inconsistent
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
