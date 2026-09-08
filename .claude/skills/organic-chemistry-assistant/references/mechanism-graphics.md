# Mechanism graphics: references, composition and visual review

Chemical validity and visual legibility are separate requirements. Do not call
a mechanism publication-ready because electron accounting or software tests
passed. A model upgrade alone does not repair a poor depiction/layout engine.

## Look at real Schemes before choosing a composition

When publication quality is requested, inspect complete figures from primary
papers, including an electron-pushing mechanism and an ordinary reaction
scheme. Record the exact DOI and Scheme number. Read the PDF image, not only
its text extraction. Borrow graphical conventions, never the paper's chemical
claims for an unrelated proposed mechanism.

References inspected for this implementation:

- [Beilstein J. Org. Chem. 2017, 13, 1230–1238, DOI 10.3762/bjoc.13.122](https://www.beilstein-journals.org/bjoc/content/pdf/1860-5397-13-122.pdf):
  Scheme 3 (PDF page 3) and Scheme 8 (PDF page 7). Ordinary reactions and the
  mechanism use the same skeletal language. In Scheme 8, local electron arrows
  stay near the bonds that change; much larger progression arrows connect
  distinct intermediate structures. A cyclic arrangement is useful for the
  proposed catalytic sequence, not automatically for a two-step acid/base example.
- [Beilstein J. Org. Chem. 2013, 9, 2635–2640, DOI 10.3762/bjoc.9.299](https://www.beilstein-journals.org/bjoc/content/pdf/1860-5397-9-299.pdf):
  Scheme 3 (PDF page 5). Compact catalyst/intermediate placement and consistent
  skeletons support reading the cycle. Its long reaction-progress arcs must not
  be mistaken for a template for lone-pair electron-pushing arrows.
- [JOC author guidelines, Appendix 2](https://researcher-resources.acs.org/publish/author_guidelines?coden=joceah),
  checked 2026-09-07: use drawing software such as ChemDraw; check lettering and
  line weight at final column size. The guide specifies Arial/Helvetica as
  suitable, lines at least 0.5 pt and black/white raster line art at 1200 dpi.
  The project's 2× PNG is a review preview, not a certified journal raster export.

ChemDraw's official installed help, **Objects → Arrows**, describes dragging
adjustment handles, editing arcs, and making an arbitrary spline with an
arrowhead. This supports independent control of curvature and end tangents;
it does not establish which software the above papers' authors used.

## Author chemical state, then art direction

1. Supply and validate the complete mapped states and electron moves using
   [MechanismSpec 2.1](mechanism-spec.md). Do not edit connectivity to fix layout.
2. Use `journal_style: "acs"`, `presentation_mode: "publication"` and
   `layout.renderer: "compact"`. This is a consistent internal style preset,
   not certification of compliance with a particular journal.
3. Let RDKit draw the skeletal structures, labels, charges and multiple bonds.
   Do not manually rebuild an entire molecule as bold `CH2`/`CH3` text and lines.
   Single-carbon reagents may retain an explicit terminal methyl label.
4. Place the reactive atoms near each other. The compact renderer sizes the
   canvas to the ink; `panel_width` no longer spreads molecules into equal cells.
   For one interacting pair, it can relocate the smaller reagent near the
   larger substrate. A '+' may be omitted inside that encounter cluster;
   spectators and product species remain explicit and separated.
5. Keep the common graph's orientation stable between steps. The renderer
   constrains unchanged mapped cores; geometry-changing bond edits can require
   new 2D coordinates. Such coordinates do not prove a transition-state angle,
   antiperiplanar alignment, or a stereochemical outcome.
6. Route local cubic electron arrows from actual lone-pair/bond objects to atom
   glyph boundaries or bonds. Implicit H and charge glyphs are obstacles, not
   substitute electron endpoints. Use full heads for pairs and fishhooks for
   single electrons. Use reaction arrows only for state progression.
7. Keep paired electron arrowheads as rounded, open V shapes with equal wings.
   The visible incoming shaft must bisect the V. On tight arcs, matching the
   derivative only at the tip can still make the shaft crowd one wing: use a
   short terminal segment on the bisector and join the cubic body tangentially.
   Check both that join and the V in the exported SVG, and audit the actual
   painted path. Do not replace the V with a filled triangle to fix alignment.
   Keep single-electron fishhooks as one open wing; reserve filled heads for
   ordinary reaction progression. The head proportions are an internal style,
   not a journal requirement.

Useful art-direction fields (they affect layout, not chemistry):

```json
{
  "molecules": [{
    "smiles": "...",
    "depiction": {
      "bond_axis": {"atom_maps": [3, 8], "angle_deg": -150},
      "rotate_deg": 0,
      "offset_bonds": [0, 0]
    }
  }],
  "composition": {"species_gap_bonds": 1.35, "auto_compact": true},
  "arrows": [{
    "from_bond": [3, 4], "to_atom_map": 4,
    "routing": {"side": "left", "lift_bonds": 0.8, "lock_side": true}
  }]
}
```

These are fragments, not a complete valid spec. `bond_axis.angle_deg` uses
chemical-coordinate angles (positive y upward); `rotate_deg` is an additional
rotation before the optional absolute bond-axis constraint. `offset_bonds`
uses page coordinates (positive y downward). Specifying an offset on either
partner disables automatic pair relocation. Use the smallest necessary control.
Routing side is relative to the source→target chord in page coordinates;
inspect the result instead of assuming that "left" means above the molecule.
Old pixel `curvature` supplies only a preferred side in compact mode.
Old cell sizes, `flip_x/y`, per-atom textual overrides, custom charge placement,
and broad teaching annotations belong to the legacy canvas; do not depend on
them for the compact artwork. The molecular graph remains authoritative.

## Inspect, correct, then read back

Check all of the following in the rendered PNG/SVG at normal reading size and
at 2× or 3× magnification:

- Molecules and ordinary reaction schemes have consistent bond and font scales.
- The nucleophilic site is near the electrophile; an arrow does not span an
  empty page cell or loop around an unrelated part of a large substrate.
- Every tail visibly touches the intended bond or sits beside the donor atom's
  lone pair, including in `HO−`, where it must not appear beside H.
- Every head ends at the intended atom/bond, clearing glyphs and formal charge.
- At high magnification, paired heads remain unfilled with round joins/caps.
  Their two wings make equal angles with the incoming shaft; the curved-to-straight
  transition has no kink, and the shaft reaches the common V tip without a gap.
- Electron arrows avoid labels, bonds, plus signs, progression arrows and each
  other. Their heads are distinguishable as pair heads or fishhooks.
- Repeated intermediates retain their identifying skeleton/branch orientation.
- No relevant atom, proton, charge or radical was hidden just to make room.

The automatic audit blocks sampled ink collisions and electron arrows longer
than 3.5 bond lengths. That threshold is an internal heuristic, not an ACS or
IUPAC rule. Passing it does not excuse a visibly poor 3.4-bond arrow. Aim for
the short local arcs seen in the references; relocate a reagent before lifting
a curve higher. If automatic layout still fails, use explicit art direction
or a chemical editor and report the remaining defect accurately.

Keep the SVG master, spec and trace together. SVG carries chemical metadata;
PNG is only the visual preview. Native editor adapters carry graphs/anchors,
but their application-specific routing and typography have not been shown to
match the SVG. Do not advertise a lossless ChemDraw/ChemDoodle round trip.
If an editor changes the chemistry or visible layout, update the spec and
regenerate the bundle before further AI reasoning.
