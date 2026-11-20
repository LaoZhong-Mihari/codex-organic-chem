# Structures from images and ambiguous input

Misreading a drawing is the most expensive error in this whole workflow: every
later check then confirms the wrong molecule. Hashed bonds read as wedges, a
missed stereocentre, a wrong ring size and an abbreviation expanded wrongly are
the usual cases.

## Reading a structure yourself

1. Transcribe to SMILES with every wedge, hash and double-bond geometry.
2. Re-read the image against your transcription, element by element: ring sizes,
   each substituent and its position, and at **every** stereocentre whether the
   bond is a wedge or a hash and which end is the narrow one.
3. `codex-chem compare --a "<smiles>"` (stereo audit): the number of
   stereocentres, the undefined ones and the formula must match the drawing. A
   printed formula, mass or name in the source is an independent check; use it.
   On a ring, `model-kit --task cis_trans --atoms a,b` says whether two
   substituents ended up on the same face, which is what the wedges encode.
4. Draw it back (`codex-chem draw --smiles ... --output png --output-file x.png`)
   and look at both images side by side. Orientation may differ; wedges at
   corresponding atoms must not.
5. State the reading in one line, including any assumption ("C7 hydroxyl read as
   hashed; relative configuration only"). An explicit assumption is what lets a
   reviewer catch a wrong diastereomer at a glance.
6. Send it to the user's review page with the image attached
   ([below](#human-confirmation)) before anything is built on it.

Never infer stereochemistry an image does not show. Resolve `R`/`X` groups and
abbreviations by graph replacement, never by string substitution.

## OCSR tools

- `chem_parse_image` (crops) and `chem_parse_scheme` (multi-compound schemes).
  Preprocessing variants and all configured engines run automatically; ranking
  rewards independent engines agreeing on one canonical structure.
- Read the top candidate's `metadata.consensus`: `multi_engine` is strong;
  `single_read` means compare the rendered candidate against the source image
  before trusting it.
- Use OCSR as a second reader: `compare` its candidate with your own
  transcription. Agreement is evidence; a `diastereomers` result names the
  centre to look at again.

## Human confirmation

Your reading and an OCSR reading looked at the same pixels and share the same
blind spots (a faint hash, an ambiguous wedge end). Their agreement lowers the
risk; it does not confirm the structure. The person who drew it can, in a few
seconds. So unless nobody is there to answer, reading an image ends here:

```bash
uv run codex-chem input-review --smiles "<your reading>" --source-image target.png \
  --note "C4 azide read as a wedge, cis to the ester" --wait
```

Give the command ten minutes or run it in the background; it returns when the
user has finished. Several structures go in one `review-batch`; MCP:
`chem_input_review(smiles, source_image, note)`,
`chem_structure_review_batch`. The page shows your structure beside the image
and the user flips a wrong wedge with one click. Continue only from the
returned `isomeric_smiles`, never from your pre-review SMILES. If the answer is
`ask_user_first`, put that one question to the user. Details:
[human-review](human-review.md).

Do not open it when nobody can answer (batch or autonomous runs, the user said
not to interrupt) or the command returns `review_disabled`: state the reading
and its assumptions and label the structure as not confirmed by a person.
Explicit machine-readable text, tool tests and assistant-authored teaching
examples need no review.
