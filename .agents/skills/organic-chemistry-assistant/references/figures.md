# Figures for a human reader

Molecule, route and compact mechanism figures share the RDKit fixed-geometry
painter. Results carry `publication_checks`; blocking issues downgrade `status`
to `blocked_for_publication`.

## Which renderer

- Molecules and single reactions: `chem_draw` (SVG master; PNG derived after the
  audit via `output_file`).
- Routes and schemes: `chem_route_figure` (`style_preset`: `acs` default, `rsc`,
  `presentation`). `route-check --figure-dir` builds the spec from a route file
  and draws every numbered compound; use that while planning. RDKit one-call
  previews are for internal checks only.
- Mechanisms with curved arrows: write a MechanismSpec 2.1
  ([mechanism-spec](mechanism-spec.md)), `chem_mechanism_validate`, then
  `chem_mechanism_render`. Compact rendering is the default;
  `layout.renderer: "legacy"` keeps the old canvas only for compatibility.
  `chem_mechanism_draft` is a rule-template hypothesis with a reaction preview,
  not a mechanism generator and not evidence for a pathway.

Chemistry-aware renderers draw all bonds and arrows; hand-built SVG, PIL or
canvas code may only do page composition. Prose belongs in the caption, not in
the graphic.

## Stale MCP server

A running MCP process can keep old code after local files change. Check the
live server's `chem_mechanism_spec_example`: if it returns a pre-2.1 example,
omits `layout.renderer: "compact"` or lacks the validator, use the checkout's
CLI for **both** validation and rendering. Compact output must report
`figure_audit.renderer: "rdkit_compact_mechanism"`; fall back to the CLI if a
stale server still returns `codex_mechanism_canvas`. A newer validator alone
does not show that the renderer was reloaded.

## Mechanism figures

- Resolve ambiguous or image-derived inputs through
  [structure-input](structure-input.md). Explicit text examples, tool tests and
  assistant-authored teaching drafts can go straight to normalization and
  validation. Manuscript claims need confirmed input.
- Render, inspect the image, and read the emitted trace back. Use
  `mechanism.trace.json` or the SVG's embedded chemical bundle to continue
  reasoning; a screenshot or an OCSR read of the figure loses electron movement
  and state identity. Keep spec, SVG and trace from the same export.
- For a publication-like mechanism, look at actual published Schemes and record
  their DOI and Scheme number. Use skeletal structures, ordinary charge glyphs,
  a shared bond scale and short local electron arrows. Move or rotate the
  reagent next to the reacting atom before trying a larger curve. Keep the
  common scaffold aligned across states. Details:
  [mechanism-graphics](mechanism-graphics.md).

## Look before you ship

Open each exported image at its intended reading size and magnified. Inspect
arrow tails and heads, implicit H, charges, lone pairs, plus signs and every
crossing, and compare with a normal reaction figure of the same structures. A
non-blank PNG, a valid graph or a passing geometric audit does not establish a
good layout. Fix visible defects and render again before showing the result.

Report status honestly: confirmed, draft or blocked, with the audit's warnings.
`blocked_for_publication` output is fixed or reported as blocked, never shipped.
A figure is publication-ready only when the input was confirmed,
`publication_checks` passes and a human has reviewed it.
