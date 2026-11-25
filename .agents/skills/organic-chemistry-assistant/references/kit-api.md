# The kit: Python API, CLI and formats

Three equivalent surfaces. Prefer Python when you need several checks or any
logic of your own; plain RDKit is always available alongside.

```python
from codex_organic_chem import kit      # run with: uv run python
```

| Python (`kit.`) | CLI (`codex-chem`) | MCP | Answers |
| --- | --- | --- | --- |
| `compare(a, b)` | `compare --a --b` | `chem_compare_structures` | `identical`, `enantiomers`, `diastereomers` (which centres, on both molecules), `same_constitution_stereo_undefined`, `tautomers`, `constitutional_isomers`, `different_formula` (delta) |
| `stereo_audit(s)` | `compare --a` | same, one SMILES | every stereo element, defined or not, CIP label, `meso`, how many isomers fit the drawing |
| `enumerate_stereoisomers(s)` | `enumerate --task stereoisomers` | `chem_enumerate` | distinct stereoisomers, enantiomer pairs, achiral (meso) forms |
| `substitution_products(s, "Cl", n_substitutions=1)` | `enumerate --task substitution --substituent Cl --n 1` | `chem_enumerate` | distinct mono/di-substitution products, with and without stereoisomers |
| `symmetry_environments(s)` | `enumerate --task environments` | `chem_enumerate` | symmetry-distinct C and H (expected signal counts), enantiotopic and diastereotopic H |
| `functional_groups(s)` | `enumerate --task groups` | `chem_enumerate` | functional groups, protecting groups, liabilities |
| `mechanism_steps_check(steps, reactants, products, conditions)` | `mechanism-steps-check --steps ...` | `chem_mechanism_steps_check` | stepwise intermediate audit: valence, radical parity, proton-transfer skeleton, deltas, cascade/product closure |
| `crop_image(path, boxes, grid=False)` | `crop-image path [--grid] [--boxes ...]` | `chem_crop_image` | coordinate ruler overlay and 2x upscaled sub-image crops from dense pages |
| `route_check(route)` | `route-check --route f.json [--figure-dir d] [--issues-only]` | `chem_route_check` | step-by-step audit (below) |
| `face(s, atom)` | `model-kit --task face --atoms i` | `chem_model_kit` | steric shielding of the Re/Si faces of a trigonal atom, per conformer, and the atoms responsible |
| `cis_trans(s, a, b)` | `model-kit --task cis_trans --atoms a,b` | `chem_model_kit` | whether the substituents on two atoms of one ring are on the same face, over all conformers |
| `distance(s, a, b)` | `model-kit --task distance --atoms a,b` | `chem_model_kit` | closest approach, its energy cost, population within a threshold |
| `dihedral(s, [a,b,c,d])` | `model-kit --task dihedral --atoms a,b,c,d` | `chem_model_kit` | anti/syn-periplanar and gauche populations |
| `planarity(s, atoms=None)` | `model-kit --task planarity` | `chem_model_kit` | deviation from a plane for ring systems or chosen atoms (xTB-refined when installed) |
| `embed_check(s)` | `model-kit --task embed` | `chem_model_kit` | whether the drawn configuration can exist in 3D |
| `isomer_energies([s1, s2], method="xtb")` | `model-kit --task isomer_energies --smiles s1 --smiles s2 --method xtb` | `chem_model_kit` | relative energies of isomers at their best conformer |
| `atom_table(s)` | `model-kit --task atoms` | `chem_model_kit` | atom indices, map numbers, neighbours |

Atom indices are those of the SMILES in reading order (implicit H come after
all heavy atoms). Print `atom_table` instead of counting by eye, or put atom-map
numbers on the atoms you care about and pass `by_map=True` / `--by-map`.

Other tools: `normalize`, `draw`, `compute` (`descriptors`, `conformers`,
`charges`, `xtb_opt`, `xtb_reactivity`, `xtb_thermo`, `crest`; `--solvent` is an
ALPB name), `reaction-analyze`, `mechanism-validate`, `mechanism-render`,
`route-figure`, `literature-search` (`--year-from`, abstracts from OpenAlex),
`parse-image`, `parse-scheme`, `doctor`.

Review page for the user ([human-review](human-review.md)): `input-review`,
`review-batch`, `route-review` open it and return a session; `review-result
--session-id ID [--wait]` collects the answer; `review-mode [auto|ask|off]`
shows or saves when it opens; `review-demo` opens two practice sessions. MCP:
`chem_input_review`, `chem_structure_review_batch`, `chem_route_review`,
`chem_structure_review_result`, `chem_review_mode`.

Hosts without a shell can enable the MCP tool `chem_python` by starting the
server with `CODEX_CHEM_ENABLE_PYTHON=1`. It runs arbitrary Python as the local
user, which is why it is off by default.

## Route format

```json
{
  "target": "C[C@H](N)CC(=O)OC",
  "compounds": {
    "1": {"smiles": "C[C@@H](O)CC(=O)OC", "name": "methyl (R)-3-hydroxybutanoate"},
    "2": "C[C@@H](OS(C)(=O)=O)CC(=O)OC"
  },
  "steps": [
    {
      "step": 1,
      "reactants": ["1"],
      "product": "2",
      "conditions": "MsCl, Et3N, CH2Cl2, 0 °C",
      "reagents": ["CS(=O)(=O)Cl"],
      "byproducts": ["Cl"],
      "why": "O-sulfonylation; the C-O bond is untouched, so the centre is retained."
    }
  ]
}
```

- `compounds`: every numbered compound, starting materials included, as SMILES
  with full stereochemistry. `reactants` and `product` refer to these ids (an
  inline SMILES is accepted).
- `conditions`: the text over the arrow. Reagent names are matched against
  compatibility rules and against the formula change, so write real reagent
  names (`TBSCl`, `DIBAL-H`, `Pd/C`, `DIAD`), not "protect" or "reduce".
- `reagents` / `byproducts` (optional SMILES): when both are given the audit
  requires the equation to close exactly and a miss is `blocking`. Use them for
  any step whose formula change is not a single textbook gain or loss.
- Atom-mapped SMILES on both sides of a step replace the heuristic atom mapping
  with yours; do this for rearrangements.

`uv run codex-chem route-check --example` prints a complete valid route. The
MCP tools `chem_route_check` and `chem_route_review` take the route inline
(`route`) or as the path of the file you wrote (`route_path`).

## Reading the audit

`status`: `blocked` (at least one `blocking` issue), `review`, or `consistent`.
Per step: `atom_balance` (formula delta and whether the conditions explain it),
`mapping` (how many product atoms were traced to reactants), `stereo` (retained,
inverted with a substituent replaced, changed without a bond change, destroyed,
created, undefined), `functional_groups` (changed and
`bystanders_that_must_survive`), `compatibility_checks`. Route level:
`final_vs_target`, `longest_linear_sequence`, `starting_materials`,
`protecting_groups_in_final_product`, `stereo_summary`.
