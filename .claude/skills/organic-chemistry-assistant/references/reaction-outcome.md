# Predicting a product, a selectivity or the better conditions

The common failure is not an illegal structure. It is committing to the textbook
pathway and then "verifying" it with checks that any well-formed product would
pass. Enumerate first, discriminate second.

## 1. Inventory

Normalize each substrate and list what can react:
`codex-chem enumerate --smiles "<substrate>" --task groups` gives functional
groups, protecting groups and liabilities (epimerisable centres, beta-leaving
groups). Add acidic and basic sites, strained bonds and symmetry yourself.

## 2. Enumerate candidates before choosing

Write every plausible product as a stereo-defined SMILES:

- the textbook outcome for this functional-group pair;
- what **this specific reagent or catalyst family** is known to do differently:
  chain walking and alkene/alkyne isomerisation before bond formation,
  migratory insertion, 1,2- versus 1,4-addition, over-reaction, ring closure onto
  a nearby group, rearrangement, a change of regiochemistry with ligand;
- the alternative regio- and stereoisomers of the same reaction.

If only one candidate comes to mind for a metal- or photocatalysed reaction, the
enumeration is not finished.

## 3. Discriminate with the right kind of evidence

| What decides it | Evidence |
| --- | --- |
| catalyst/ligand/reagent-controlled pathway | precedent: `literature-search` and web search for the exact catalyst system with this substrate class. No structure tool can settle this. |
| which of several similar sites reacts | `compute --task xtb_reactivity` (Fukui indices, charges per atom) |
| which face or how close two atoms get | `model-kit --task face\|distance\|dihedral` |
| equilibrating conditions (thermodynamic product) | `model-kit --task isomer_energies` (use `--method xtb`; gaps under ~1.5 kcal/mol decide nothing) |
| stereospecific mechanism (SN2, E2, syn/anti addition, pericyclic rules) | reason it, then confirm the drawn product with `compare` against the alternative |

Without retrieval, say that the catalyst-specific behaviour is being inferred
from memory of the reaction class, and keep the runner-up alive.

## 4. Check the candidate you chose

- `reaction-analyze` on the full equation: `mass_charge_balance` must close with
  the stated stoichiometry.
- `compare` with one SMILES: every stereo element defined, the count you expect.
- If the question asks for a number, `enumerate` or code it.

Passing these shows the answer is well formed. It does not show it is observed.

## 5. Answer

- The predicted product(s) as SMILES, with the mechanism in two or three lines.
- The runner-up and the observation that would tell them apart (a diagnostic
  NMR signal or coupling, HRMS, an isotope-labelling outcome).
- What the prediction rests on: precedent with DOI and grade, a calculation with
  its number, or inference; and the condition under which it would flip.
