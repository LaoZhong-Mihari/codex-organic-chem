# Putting a person in the loop

Your checks catch bookkeeping. They cannot tell you that the wedge in the user's
drawing was a hash, or that step 6 failed in their hands last month. The review
page is how the user tells you, in seconds and without typing SMILES: it opens
in their browser, shows what you are about to rely on next to its source, and
lets them correct it in place. Their corrections come back as data.

## When to open it

Call the command at these points. Do not decide for the user whether they want
it: the command reads their setting (`codex-chem review-mode`: `auto`, `ask`,
`off`) and tells you what to do.

1. **Structures you read from an image, a scheme or an ambiguous name**, before
   anything is built on them. One session for all of them.
2. **A route the user will act on**, after your own `route-check` has no
   `blocking` issue. Present the route in chat as well; the page is where they
   correct it.
3. **You cannot decide between readings** (two diastereomers, an unclear ring
   size): one item with `candidates` and a `question`.
4. **The user asks** to see, check or edit structures: pass `--yes`
   (`confirmed=true`).

Skip it when no person is there to answer (batch or headless runs, the user
said not to interrupt), and for structures the user typed as SMILES or that
are textbook examples of your own. Then state your reading and assumptions and
label the result "not confirmed by a person".

What the command returns:

- `pending` with `session_id`: the page is open (`browser_opened: false` means
  give the user `review_url`). Say in one line what you want checked.
- `ask_user_first`: put `user_prompt` to the user as a single question. On yes,
  call again with `--yes`; on no, continue unconfirmed.
- `review_disabled`: do not open or offer it; continue unconfirmed.

## Opening a review

```bash
# your reading of one image, beside the image
uv run codex-chem input-review --smiles "<your reading>" --source-image crop.png --note "..." --wait
# several structures
uv run codex-chem review-batch --input items.json --lang zh
# let OCSR read the image; its other candidates become alternatives
uv run codex-chem input-review --image-path crop.png --note "..."
# a route (same file as route-check); --image shows the target's source drawing beside it
uv run codex-chem route-review --route route.json --intro "..." --image target.png
```

MCP: `chem_input_review`, `chem_structure_review_batch`, `chem_route_review`
(takes `route` or `route_path`). `--wait` blocks up to nine minutes for the
answer, so give the command a ten-minute timeout; without it the command
returns at once and `review-result` collects the answer.

`items.json`:

```json
{"title": "Scheme 2 structures",
 "intro": "Read from your scheme. The wedges on 7 and 9 are the ones I am least sure of.",
 "images": ["scheme2.png"],
 "items": [
  {"id": "7", "label": "7", "smiles": "...", "image_path": "crops/7.png",
   "note": "C4-OH read as hashed; ring read as six-membered",
   "question": "Is the C4 hydroxyl a wedge or a hash?",
   "candidates": [{"label": "C4 epimer", "smiles": "..."}]}]}
```

- `image_path` is the crop the structure came from; it is shown beside your
  structure, which is what makes a wrong wedge visible at a glance. Always
  attach it when there is one.
- `note` is your reading and its assumptions; `question` is the one thing you
  most want answered; `candidates` are the readings you could not rule out.
- Use the user's language: `--lang zh` (`language: "zh"`) when they write
  Chinese.
- Send only what needs a person. Do not ask them to confirm atom counts or
  anything else a check already settles.

On the page the user can click a stereocentre to flip it, set an open centre or
double bond with R/S and E/Z buttons, pick an alternative, edit SMILES, paste a
Molfile or draw in Ketcher, and leave a note per item. On a route they can also
change the conditions of a step, mark each step right or wrong, answer each
audit prompt and approve the route or send it back; every edit re-runs the
audit in front of them.

## Collecting the answer

```bash
uv run codex-chem review-result --session-id <id> --wait    # waits up to 9 min
```

MCP: `chem_structure_review_result(session_id, wait=true)`.

- The page stays open after the command that created it returns, and a wait
  that times out (`wait_timed_out`) leaves the session open. Where the host can
  run a command in the background, start the wait there and continue with work
  that does not depend on the answer.
- Structures from an image are a gate: do not build a route, mechanism or
  expensive calculation on them while the session is `pending`. If the wait
  times out, tell the user the page is waiting for them and stop there.
- A route review does not block your answer: give the route, say the page is
  open, and revise when the result arrives.

Using the result:

- Structures: continue only from each item's returned `isomeric_smiles`.
  `modified` items say what changed in `changed_from_input` (for example
  `C7: R -> S`): carry that correction into everything derived from the old
  structure. `rejected`, `skipped` and `error` items are not confirmed; read the
  `comment`, re-derive and ask again.
- Routes: continue from `route_result.route`, the reviewer's version. Address
  every step in `step_decisions` marked `flagged`, every entry in
  `issue_answers` marked `problem`, and both comment fields; then re-run
  `route-check`. `verdict: approved` approves the bookkeeping they looked at,
  not the chemistry they did not comment on.
- Say in one line what the person changed, the same way you report what your
  own checks made you change.
