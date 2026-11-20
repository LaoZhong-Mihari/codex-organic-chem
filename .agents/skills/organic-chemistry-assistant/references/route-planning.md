# Planning and auditing a synthetic route

A route is finished when it has survived its own audit, not when it reads well.
Work in this order. Do not ask the user to approve each disconnection unless
they want to steer; commit to one route and state assumptions in one line each.

## 1. Fix the target

- One SMILES with every stereocentre and double bond defined.
  `codex-chem compare --a "<target>"` must report `fully_defined: true` and the
  number of stereocentres you expect.
- Target given as an image: [structure-input](structure-input.md) first. A route
  to the wrong diastereomer is a failed route however good the chemistry.

## 2. Strategy before tactics

Write two or three strategically different plans, five lines each: key bonds
disconnected, order of ring formation, and where every stereocentre comes from
(chiral pool, substrate control, reagent or catalyst control). At least one must
try to be shorter: convergence, hidden symmetry, a cascade, a catalytic
asymmetric step instead of an auxiliary, a less obvious key disconnection.
Pick one, and say in a line each why the others lose. Defaulting to the most
precedented sequence is how routes end up long and uninteresting.

## 3. Write the route as data

Number every compound in order of appearance, starting materials included, and
write the route file ([format](kit-api.md#route-format)): each step has
`reactants`, `product`, `conditions` (the short form written over an arrow) and
`why` (one sentence, including how any new stereocentre is set). List atom-
contributing `reagents` and `byproducts` as SMILES wherever the audit should
close the equation exactly.

## 4. Audit, fix, repeat

```bash
uv run codex-chem route-check --route route.json --figure-dir out
```

- Clear every `blocking` issue: these are wrong structures, not opinions.
- Answer every `review` item in one line or change the route:
  - `formula_change_*`: a structure has the wrong atom count, or the conditions
    omit the reagent that supplies the atoms. Fix or declare `reagents`/`byproducts`.
  - `configuration_changed_without_bond_change`,
    `inversion_where_only_the_heteroatom_was_derivatised`,
    `retention_under_inverting_conditions`: a wedge is wrong, or the mechanism is
    not what the conditions say.
  - `atoms_relocated`: a substituent moved between reactant and product. Either
    name the rearrangement or correct the position.
  - `undefined_stereo_in_product`: draw the isomer that forms, or say why the
    mixture is inconsequential (the centre is destroyed in step N).
  - `selectivity_among_equal_groups`, `bystander_group_vs_reagent`: say what
    protects the group that must survive, reorder steps, or change the reagent.
  - `starting_material_availability`: check it is purchasable and affordable, or
    add the steps that make it.
- Open `out/route.png`. Every drawn structure must be the one you meant.

## 5. The chemistry the audit cannot judge

For each step, with the audit's bystander list in front of you:

- Do the other functional and protecting groups survive these exact conditions?
- Is every stereocentre of the target set somewhere, and is none scrambled
  later (centres alpha to carbonyls under base or activation, centres between
  two acceptors, allylic centres under Pd/H2, retro-aldol of beta-hydroxy
  carbonyls, beta-elimination)?
- Does a 3D claim decide the step? Then look: `model-kit --task face` for the
  open face of a carbonyl or alkene, `distance` for whether a ring closure can
  be reached, `embed`/`isomer_energies` for a strained or equilibrating ring
  fusion. State that this is a steric ground-state argument.
- Does site selectivity decide the step? `compute --task xtb_reactivity`.

## 6. Precedent and supply

When retrieval is allowed, search each key or unusual step
(`literature-search`, then web search): exact transformation first, then the
closest analogue. Attach DOI and grade (exact precedent, close analogue, same
class, review). A step that failed in a closely related system is the most
valuable thing a search can find. Without retrieval, say so and mark recalled
precedent as unverified; never write a citation from memory as if checked.

## 7. Present

1. Retrosynthetic overview, five lines or fewer.
2. Numbered steps: Reactants / Product / Conditions / Why.
3. `Changed after checking:` one line listing what the audit made you fix.
4. Total steps, longest linear sequence, protecting-group operations.
5. The three steps most likely to fail, each with a fallback that changes the
   **tactic** (a different bond construction or order of events), not merely the
   catalyst or solvent.
6. What was verified (audit status, calculations, literature grade) and that
   nothing has been run in a laboratory.

Then open the route for the user, `uv run codex-chem route-review --route
route.json --intro "<what you most want checked>"`
([human-review](human-review.md)), and say so in one line. They can correct a
structure or the conditions of a step, mark steps and send the route back.
This does not hold up your answer. When they say they are done, or a
background `review-result --session-id ID --wait` returns, continue from their
version of the route (`route_result.route`), answer every flagged step and
comment, and re-run the audit.
