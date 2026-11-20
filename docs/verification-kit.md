# Verification kit

[Back to README](../README.md)

## Verification Kit

The skill asks the model to reason first and then try to break its own answer
with checks that can fail. Parsing a SMILES the model just wrote cannot fail in
an informative way; these can:

```bash
uv run codex-chem compare --a "C[C@H](O)CC" --b "C[C@@H](O)CC"          # enantiomers
uv run codex-chem compare --a "CC1(C)C2CCC1(C)C(=O)C2"                    # stereo audit: 2 isomers fit
uv run codex-chem enumerate --smiles "CC1CCCCC1" --task substitution     # 5 constitutional, 12 with stereo
uv run codex-chem enumerate --smiles "OCC(O)CO" --task environments      # 2 C and 5 H environments
uv run codex-chem route-check --example > route.json
uv run codex-chem route-check --route route.json --figure-dir out --issues-only
uv run codex-chem model-kit --task face --smiles "O=C1C[C@@H]2CC[C@H]1C2" --atoms 1
uv run codex-chem model-kit --task isomer_energies --method xtb \
  --smiles "C1CC[C@H]2CCCC[C@@H]2C1" --smiles "C1CC[C@H]2CCCC[C@H]2C1"
```

`route-check` reads a route as data (numbered compounds as SMILES, steps with
reactants, product and conditions) and reports, per step: the formula change
and whether a reagent named in the conditions accounts for it; each
stereocentre as retained, inverted, destroyed, created or undefined, judged by
spatial arrangement under an atom mapping rather than by CIP letters; the
groups that are drawn unchanged although the reagent class attacks them; and
atoms that cannot be traced from reactant to product. At route level it
compares the last product with the target centre by centre and lists step
count, longest linear sequence, starting materials and leftover protecting
groups. `blocking` issues are facts; `review` items are heuristic prompts.
The format and every field are described in the
[kit reference](../.agents/skills/organic-chemistry-assistant/references/kit-api.md).

`model-kit` covers what a physical model is used for: which face of a carbonyl
or alkene is more shielded and by which atoms, whether two atoms can approach,
anti-periplanar populations, planarity, whether a drawn ring fusion can exist,
and which isomer is lower in energy. It is a steric ground-state model and says
so in every result.

```bash
uv run python scripts/benchmark_route_check.py --output-dir outputs/route_check_benchmark
```

corrupts one compound at a time in three curated routes (inverted centre,
dropped wedge, E/Z flip, lost carbon, extra CH2, O/N swap, shifted substituent,
omitted last step) and counts how often a tool's output changes. Current
result: 280 of 282 corruptions change the audit at blocking or review level; the
earlier balance and normalization checks change for 67 and report every correct
step as unbalanced. This shows that the feedback discriminates. It does not
show that a model uses it; that needs a model-in-the-loop comparison.

A `chem-route-critic` subagent (`.claude/agents/`) gives a route an independent
second reading in Claude Code. Hosts without a shell can enable the opt-in MCP
tool `chem_python` (see [DEPENDENCIES.md](../DEPENDENCIES.md)).

## Stepwise Synthesis Suggestions

The default route workflow is now: commit to one full route, write it as data
and run `route-check` (see the skill's route-planning reference).
`synthesis-suggest` remains for interactive brainstorming, one decision layer at
a time, when the user wants to steer each disconnection.

```bash
codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O"
codex-chem synthesis-suggest --target-smiles "CC(=O)Oc1ccccc1C(=O)O" --confirmed
```

Each option includes local tool data and, when enabled, literature records
(OpenAlex abstracts, Crossref metadata). The next step is to choose or reject one disconnection before any
route expansion.
