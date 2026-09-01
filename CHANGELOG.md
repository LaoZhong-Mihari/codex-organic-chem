# Changelog

All notable changes to this project are documented here.

## [Unreleased]

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
