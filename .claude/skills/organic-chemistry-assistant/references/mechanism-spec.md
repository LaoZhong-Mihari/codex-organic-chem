# Mechanisms that can be drawn and read back

This is the project's MechanismSpec 2.1 contract, not a claim of a universal
IUPAC interchange format. The diagram uses conventional electron-pair arrows
and fishhooks; the JSON preserves the chemistry for validation and further use.

## Authoring a step

Start from [the complete SN2 spec](../assets/sn2.mechanism.json). Preserve its
field names; free-text electron descriptions are not arrow anchors.

| Field | Meaning |
| --- | --- |
| `spec_version` | `"2.1"` requires complete supported semantic validation. |
| `panels` | One elementary step per panel, containing its before and after structures. |
| `molecules[].smiles` / `molfile` | One connected species per entry; mapped atoms, actual charges and stereochemistry. |
| `starting_state.molecule_indices` | One-based indices of ONLY the starting species in this panel. |
| `expected_state.molecule_indices` | Disjoint one-based indices of ONLY its resulting species. |
| `reaction_arrow_after` | Layout split; keep consistent with the state selections. |
| `arrows` | Electron movements on the starting state, committed together. |
| `graph_edits` | Complete changed bonds, formal charges and radical populations. |
| `evidence`, `uncertainty` | Concise provenance and unresolved chemical assumptions. |

Every atom must have a positive map number, unique within each state and
conserved between states. Reusing the map on the product side is expected;
renumbering an atom between steps is not. Indices are one-based molecule
positions, not RDKit atom indices. Use explicit `from_molecule_index` and
`to_molecule_index` whenever maps appear on both sides of a panel.

A transferred proton must be a mapped vertex, e.g. `[O:4][H:5]`.
`[OH:4]` has an attached H count but no addressable H atom: it cannot represent
the identity of a moving proton. Preserve non-transferring implicit/bracket H
counts. Include the acid/base molecule rather than inventing a bare proton
unless the model being discussed explicitly uses it. Include counterions and
byproducts consistently; metadata does not repair an unbalanced graph.

Adjacent steps must agree exactly on the previous product and next starting
state. Carry a reagent through earlier panels as a spectator when necessary.
This avoids silently adding/removing chemical matter between steps. Do not
include alternative pathways or resonance contributors as sequential states:
make separate candidate specs and label resonance separately from reaction.

## Electron anchors and edits

Use `kind: "electron_pair", electron_count: 2` for a full arrowhead, or
`kind: "single_electron", electron_count: 1` for a fishhook. Coordinates and
curvature control appearance only, never chemical identity.

| Movement | Required fields and effect |
| --- | --- |
| Lone pair to atom | `from_lone_pair_atom_map`, `to_atom_map`; forms/increases their bond. |
| Bond to endpoint | `from_bond: [a,b]`, `to_atom_map: b`; decreases bond order, pair ends on b. |
| Bond to adjacent bond | `from_bond: [a,b]`, `to_bond: [b,c]`; decreases a-b and increases b-c together. |
| Homolysis | Two fishhooks from the same bond, one to each endpoint; one radical on each atom. |

Add a matching `lone_pairs` entry for every LP source, with `molecule_index`,
`atom_map`, chemically correct `count`, and typically `visible_count: 1`.
The other LPs can remain hidden. Charges come from the molecule; a display
override must agree with the structure. Avoid arbitrary labels that could
change the meaning of an element, hydrogen count or charge.

Graph edits use numeric orders, not bond-type strings:

```json
[
  {"type": "form_bond", "atom_maps": [1, 2], "order": 1, "source_arrow": 1},
  {"type": "break_bond", "atom_maps": [2, 3], "source_arrow": 2},
  {"type": "formal_charge", "atom_map": 1, "from": -1, "to": 0},
  {"type": "formal_charge", "atom_map": 3, "from": 0, "to": -1}
]
```

For a bond-order change use `type: "change_bond", atom_maps: [a,b], from: 2,
to: 1`. For radicals use `type: "radical", atom_map: a, from: 0, to: 1`.
`source_arrow` is an optional one-based provenance link; omit it when multiple
arrows contribute rather than providing a false link. A `proton_transfer`
annotation is supplemental: still list the actual H bond and charge edits.
Do not author redundant `bond_changes` / `charge_changes` summaries: the trace
computes those from the structures independently.

An SN2 step has attack and leaving-group arrows together. Carbonyl addition
has attack and pi-to-oxygen movement together; protonation is a subsequent
step with an explicit donor. E2 has base-to-H, C-H-to-C-C and C-X-to-X arrows
together. Do not split a concerted step into pentavalent-carbon intermediates.

## Validate, view, read back

```bash
uv run codex-chem mechanism-validate --spec mechanism.spec.json
uv run codex-chem mechanism-render --spec mechanism.spec.json --output-dir out
```

Check `chem_mechanism_spec_example` on the live MCP server. If it returns a
pre-2.1 example, lacks `chem_mechanism_validate`, or omits `semantic_checks`, use
the updated checkout's CLI for BOTH validation and rendering. Do not pass a
newly validated spec to an old renderer. A long-running MCP process requires
reconnection to load changed code; a new on-disk skill alone does not reload it.
For the current drawing implementation, also check that the example advertises
`layout.renderer: "compact"` and a rendered result reports
`figure_audit.renderer: "rdkit_compact_mechanism"`. A 2.1 semantic validator can
still be paired with an old custom painter in a process started before this update.

Follow [the graphical reference and review workflow](mechanism-graphics.md).
Compact mode uses the same RDKit painter as ordinary reaction schemes, places
small reagents near reactive sites, and checks electron-arrow ink collisions.

Keep these outputs together:

- `mechanism.spec.json`: authored chemical state and layout instructions.
- `mechanism.svg` and its derived PNG: human-readable chemical figure.
- `mechanism.trace.json`: resolved before/after atom/bond state, observed
  changes, electron moves, validation result, and the source spec SHA-256.

SVG `<metadata id="codex-mechanism">` embeds the same spec and trace as a JSON
bundle. Parse XML, read the metadata text, then parse JSON; do not infer
chemistry from SVG paths. Hash canonical UTF-8 spec JSON with sorted keys,
compact separators and unescaped Unicode. External visual edits can make the
embedded bundle stale: re-export from an updated spec after editing. The hash
links chemical data, not arbitrary edits to the visible SVG paths.

Read `semantic_checks` separately from `publication_checks.figure_audit`.
Inspect electron origins/endpoints, legible H/charge/radical labels, crossings,
and panel correspondence. A valid graph and a nonblank PNG still need visual
inspection; neither is evidence that the pathway occurs experimentally.

## Current limits and useful follow-on reasoning

The semantic checker supports localized main-group pair transfers and paired
homolysis. General radical propagation, aromatic electron replay, metal
coordination/redox cycles, pericyclic orbital selection rules, transition
states, energies and rates need separate analysis. Unsupported cases must stay
explicitly unsupported. The export pipeline conservatively blocks stereo wedges,
E/Z and isotope cases whose complete depiction/adapter path has not been verified.
Use a chemical editor for these;
do not strip the information to obtain an `ok` result. Native editor exports
are review adapters, not verified lossless editor round trips.

For further chemical reasoning, read the resolved product's charge, radical
population and bonds, identify remaining donor/acceptor sites, and propose the
next step conditionally. A useful explanation names the changed atom maps,
why the reagent/conditions could favor that step, a close competing path, and
what evidence would distinguish them. For an alkoxide after carbonyl addition,
protonation requires a donor in the state; for E2, antiperiplanar access and
competition with substitution remain separate questions.

The repository's `scripts/benchmark_mechanisms.py` runs five curated mechanism
classes and deliberate errors. These test the representation/toolchain, not a
model's success rate. To compare model versions, use the same held-out prompts,
conditions, generation budget, model identifiers, and independent rubric;
retain first-attempt outputs and report revisions separately.

Sources for representation and implementation:
[ChemDoodle JSON chemical-object pushers](https://web.chemdoodle.com/docs/chemdoodle-json-format)
and [RDKit molecule and stereochemistry documentation](https://www.rdkit.org/docs/RDKit_Book.html).
