import json
import os
import subprocess
import sys
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from codex_organic_chem import review_chem as rc
from codex_organic_chem import review_daemon
from codex_organic_chem import service
from codex_organic_chem import structure_review as sr
from codex_organic_chem.review_demo import MISREAD_PROLINE, TRUE_PROLINE, demo_items, demo_route
from codex_organic_chem.review_route import RouteReview

FIXTURES = Path(__file__).parent / "fixtures" / "routes"
PROLINE = "COC(=O)[C@@H]1C[C@@H](O)CN1C(=O)OC(C)(C)C"


def _post(url: str, payload: dict, token: str | None = None) -> dict:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["X-Review-Token"] = token
    with urlopen(Request(url, data=json.dumps(payload).encode(), headers=headers), timeout=20) as response:
        return json.loads(response.read())


def _api(session: dict, path: str) -> str:
    return f"{session['api_url']}/api/review-sessions/{session['session_id']}{path}"


# ---------------------------------------------------------------- chemistry


def test_describe_lists_stereo_elements_in_the_atom_order_of_the_input():
    described = rc.describe(PROLINE)
    assert described["status"] == "ok" and described["smiles"] == PROLINE
    assert described["formula"] == "C11H19NO5" and described["n_undefined"] == 0
    assert [(e["id"], e["name"], e["cip"]) for e in described["stereo"]] == [("a4", "C5", "S"), ("a6", "C7", "R")]
    assert "atom-selector" in described["svg"]  # atoms are clickable in the page


def test_one_click_edits_change_exactly_what_was_clicked_and_keep_atom_names():
    flipped = rc.apply_edit(PROLINE, "invert", "a6")
    after = rc.describe(flipped, reference=PROLINE)
    assert after["versus_reference"] == {"relation": "diastereomers", "changed_atoms": [6], "details": ["C7: R -> S"]}
    assert [e["name"] for e in after["stereo"]] == ["C5", "C7"]  # numbering survives the edit
    assert rc.describe(rc.apply_edit(flipped, "invert", "a6"), reference=PROLINE)["versus_reference"]["relation"] == "identical"
    assert rc.describe(rc.apply_edit(PROLINE, "mirror"), reference=PROLINE)["versus_reference"]["relation"] == "enantiomers"

    cleared = rc.describe(rc.apply_edit(PROLINE, "clear", "a4"), reference=PROLINE)
    assert cleared["n_undefined"] == 1 and cleared["versus_reference"]["details"] == ["C5: no longer defined"]
    for wanted in ("R", "S"):
        assigned = rc.describe(rc.apply_edit(cleared["smiles"], "set", "a4", wanted))
        assert assigned["stereo"][0]["cip"] == wanted and assigned["n_undefined"] == 0


def test_double_bond_geometry_can_be_set_flipped_and_cleared():
    open_diene = rc.describe("CC=CC(O)C/C=C\\C")
    assert {e["id"]: e["specified"] for e in open_diene["stereo"]} == {"a3": False, "b1-2": False, "b6-7": True}
    for wanted in ("E", "Z"):
        assigned = rc.describe(rc.apply_edit(open_diene["smiles"], "set", "b1-2", wanted))
        assert {e["id"]: e["cip"] for e in assigned["stereo"]}["b1-2"] == wanted
    e_isomer = rc.apply_edit(open_diene["smiles"], "set", "b1-2", "E")
    flipped = rc.describe(rc.apply_edit(e_isomer, "invert", "b6-7"), reference=e_isomer)
    assert {e["id"]: e["cip"] for e in flipped["stereo"]}["b6-7"] == "E"
    assert flipped["versus_reference"]["details"] == ["C7=C8: E/Z flipped"]
    assert rc.describe(rc.apply_edit(e_isomer, "clear", "b1-2"))["n_undefined"] == 2


def test_describe_reports_bad_input_and_ring_fixed_centres_as_data():
    assert rc.describe("not a smiles")["status"] == "invalid"
    assert rc.describe("CCO", reference="CCN")["versus_reference"]["formula_delta"] == "-1H -1N +1O"
    norbornane = rc.describe("C1CC2CCC1C2")
    assert norbornane["n_undefined"] == 0 and all(e["constrained"] for e in norbornane["stereo"])
    with pytest.raises(rc.StructureError):
        rc.apply_edit("CCO", "set", "a1", "R")
    # an azide is linear: its N=N is not an open double bond
    assert rc.describe("C[C@H](N=[N+]=[N-])CC")["n_stereo_elements"] == 1


# ------------------------------------------------------- structure sessions


def test_structure_session_carries_context_for_the_reviewer(tmp_path):
    image = Path(demo_items()[1]["image_path"])
    session = service.chem_structure_review_batch(
        items=[
            {
                "id": "2",
                "label": "compound 2",
                "smiles": MISREAD_PROLINE,
                "image_path": str(image),
                "note": "C4-OH read as hashed",
                "question": "wedge or hash?",
                "candidates": [TRUE_PROLINE, {"label": "broken", "smiles": "C1CC"}],
            }
        ],
        title="Scheme 1",
        intro="check the wedges",
        language="zh",
    )
    assert session["kind"] == "structures" and session["title"] == "Scheme 1" and "lang=zh" in session["review_url"]
    item = session["items"][0]
    assert item["note"] == "C4-OH read as hashed" and item["question"] == "wedge or hash?"
    assert [c["label"] for c in item["candidates"]] == ["Alternative 1"]
    assert any("alternative 2" in warning for warning in session["warnings"])
    assert "svg" not in json.dumps(session)  # drawings stay out of the assistant's context

    with urlopen(f"{_api(session, '/assets/' + item['image'])}?token={session['review_token']}") as response:
        assert response.headers["Content-Type"] == "image/png" and response.read()[:4] == b"\x89PNG"
    with pytest.raises(HTTPError):
        urlopen(_api(session, "/assets/" + item["image"]) + "?token=wrong")

    page = urlopen(session["review_url"]).read().decode()
    assert "Structure review" in page and "click a stereocentre to flip it" in page


def test_page_api_edits_and_returns_what_the_reviewer_changed():
    session = service.chem_structure_review_batch(items=[{"id": "2", "smiles": MISREAD_PROLINE}])
    token = session["review_token"]
    with pytest.raises(HTTPError):
        _post(_api(session, "/describe"), {"smiles": "CCO"})  # no token
    edited = _post(_api(session, "/edit"), {"smiles": MISREAD_PROLINE, "op": "invert", "target": "a6", "reference": MISREAD_PROLINE}, token)
    assert edited["versus_reference"]["relation"] == "diastereomers"
    _post(_api(session, "/comment"), {"comment": "scan is faint"}, token)
    final = _post(
        _api(session, "/items/2"),
        {"status": "confirmed", "reviewed_smiles": edited["smiles"], "comment": "OH is a wedge in my drawing"},
        token,
    )
    assert final["status"] == "completed"

    result = service.chem_structure_review_result(session["session_id"])
    item = result["items"][0]
    assert item["status"] == "modified" and item["comment"] == "OH is a wedge in my drawing"
    assert rc.describe(item["isomeric_smiles"], reference=TRUE_PROLINE)["versus_reference"]["relation"] == "identical"
    assert item["changed_from_input"]["relation"] == "diastereomers"
    assert result["comment"] == "scan is faint"
    assert result["summary"]["changed_items"] == ["2"] and result["summary"]["all_usable"] is True
    assert "never your pre-review SMILES" in result["next_action"]


def test_rejected_and_reopened_items():
    session = service.chem_structure_review_batch(items=[{"id": "a", "smiles": "CCO"}, {"id": "b", "smiles": "CCN"}])
    token = session["review_token"]
    rejected = sr.submit_review_item(session["session_id"], "a", {"status": "rejected", "comment": "ring is 7-membered"}, token)
    assert rejected["items"][0]["status"] == "rejected" and rejected["items"][0]["canonical_smiles"] is None
    reopened = sr.submit_review_item(session["session_id"], "a", {"status": "pending"}, token)
    assert reopened["items"][0]["status"] == "pending" and reopened["status"] == "pending"
    sr.submit_review_item(session["session_id"], "a", {"status": "rejected", "comment": "ring is 7-membered"}, token)
    final = sr.submit_review_item(session["session_id"], "b", {"status": "confirmed", "reviewed_smiles": "CCN"}, token)
    assert final["status"] == "completed"
    assert final["summary"]["unusable_items"] == ["a"] and final["summary"]["all_usable"] is False


def test_a_short_wait_leaves_the_session_open_and_results_are_on_disk():
    session = service.chem_structure_review_batch(items=[{"id": "a", "smiles": "CCO"}], timeout_s=60)
    waited = service.chem_structure_review_result(session["session_id"], wait=True, timeout_s=0.05)
    assert waited["status"] == "pending" and waited["wait_timed_out"] is True
    sr.submit_review_item(session["session_id"], "a", {"status": "confirmed", "reviewed_smiles": "CCO"}, session["review_token"])

    stored = json.loads((sr.reviews_dir() / f"{session['session_id']}.json").read_text())
    assert stored["status"] == "completed" and stored["items"][0]["canonical_smiles"] == "CCO"
    # another process (or this one after a restart) reads the same result
    with sr._SESSIONS_LOCK:
        sr._SESSIONS.pop(session["session_id"])
    assert service.chem_structure_review_result(session["session_id"], wait=True, timeout_s=1)["status"] == "completed"
    with pytest.raises(KeyError):
        service.chem_structure_review_result("../../etc/passwd")


# ----------------------------------------------------------- route sessions


def test_route_review_reaudits_after_a_correction_and_returns_the_reviewers_route():
    session = service.chem_route_review(route=demo_route(), intro="check step 5")
    assert session["kind"] == "route" and session["status"] == "pending"
    token = session["review_token"]

    page = json.loads(urlopen(_api(session, "") + f"?token={token}&view=page").read())
    route = page["route"]
    codes = {step["key"]: [issue["code"] for issue in step["issues"]] for step in route["steps"]}
    assert "retention_under_inverting_conditions" in codes["5"]
    assert "configuration_changed_without_bond_change" in codes["6"]
    assert route["compounds"]["7"]["svg"].startswith("<?xml") or "<svg" in route["compounds"]["7"]["svg"]

    wrong = route["compounds"]["7"]["smiles"]
    fixed = _post(_api(session, "/edit"), {"smiles": wrong, "op": "invert", "target": "a6"}, token)["smiles"]
    page = _post(_api(session, "/route/compound"), {"id": "7", "smiles": fixed}, token)
    codes = {step["key"]: [issue["code"] for issue in step["issues"]] for step in page["route"]["steps"]}
    assert "retention_under_inverting_conditions" not in codes["5"] and codes["6"] == []
    assert page["route"]["compounds"]["7"]["changed"] and page["route"]["n_edited_compounds"] == 1
    step5 = next(step for step in page["route"]["steps"] if step["key"] == "5")
    assert step5["product_svg"] and any(chip["kind"] == "inverted" for chip in step5["stereo"])

    _post(_api(session, "/route/step"), {"key": "3", "decision": "flagged", "comment": "use LiOH-free conditions"}, token)
    _post(_api(session, "/route/step"), {"key": "5", "decision": "approved"}, token)
    page = _post(_api(session, "/route/step"), {"key": "4", "conditions": "Ms2O, pyridine, 0 °C"}, token)
    issue_key = next(step for step in page["route"]["steps"] if step["key"] == "3")["issues"][0]["key"]
    _post(_api(session, "/route/issue"), {"key": issue_key, "answer": "problem"}, token)
    done = _post(_api(session, "/complete"), {"verdict": "needs_changes", "comment": "rework step 3"}, token)
    assert done["status"] == "completed"

    result = service.chem_structure_review_result(session["session_id"])
    outcome = result["route_result"]
    assert outcome["verdict"] == "needs_changes" and outcome["needs_revision"] is True
    assert outcome["route"]["compounds"]["7"]["smiles"] == fixed
    assert outcome["changed_compounds"]["7"]["relation"] == "enantiomers" or outcome["changed_compounds"]["7"]["relation"] == "diastereomers"
    assert outcome["changed_conditions"]["4"]["to"] == "Ms2O, pyridine, 0 °C"
    assert {"step": "3", "decision": "flagged", "comment": "use LiOH-free conditions"} in outcome["step_decisions"]
    assert outcome["issue_answers"] == [{"issue": issue_key, "answer": "problem", "comment": ""}]
    assert result["summary"]["flagged_steps"] == ["3"] and result["comment"] == "rework step 3"
    assert "route_result.route" in result["next_action"]
    assert "svg" not in json.dumps(result)


def test_route_review_survives_an_unparseable_edit_and_rejects_bad_routes():
    review = RouteReview(json.loads((FIXTURES / "roche_polyketide_fragment.json").read_text()))
    assert review.view()["n_blocking"] == 0
    with pytest.raises(ValueError):
        review.set_compound("2", "C1CC")
    review.reset_compound()
    assert review.result()["changed_compounds"] == {}
    assert service.chem_route_review(route={"steps": []})["status"] == "error"


# -------------------------------------------------- review mode and the gate


def test_review_mode_decides_whether_the_page_opens(monkeypatch):
    assert service.chem_review_mode()["review_mode"] == "auto"
    assert service.chem_review_mode("ask")["review_mode"] == "ask"
    asked = service.chem_structure_review_batch(items=[{"smiles": "CCO"}])
    assert asked["status"] == "ask_user_first" and "Open it?" in asked["user_prompt"]
    assert service.chem_structure_review_batch(items=[{"smiles": "CCO"}], confirmed=True)["status"] == "pending"
    static = service.chem_input_review(smiles="CCO")
    assert "session_id" not in static and static["review_mode"] == "ask" and "Open it?" in static["user_prompt"]

    assert service.chem_review_mode("off")["review_mode"] == "off"
    off = service.chem_route_review(route=demo_route())
    assert off["status"] == "review_disabled" and "not confirmed by a person" in off["next_action"]
    assert service.chem_route_review(route=demo_route(), confirmed=True)["status"] == "pending"

    monkeypatch.setenv("CODEX_CHEM_REVIEW", "auto")  # the environment wins over the saved setting
    status = service.chem_review_mode()
    assert status["review_mode"] == "auto" and status["saved_setting"] == "off"
    assert service.chem_review_mode("sometimes")["status"] == "error"


# ---------------------------------------------------- detached review server


def _cli(*args: str) -> subprocess.CompletedProcess:
    environment = os.environ.copy()
    environment["CODEX_CHEM_REVIEW_OPEN_BROWSER"] = "0"
    return subprocess.run(
        [sys.executable, "-m", "codex_organic_chem", *args], capture_output=True, text=True, env=environment, timeout=120
    )


def test_cli_review_page_outlives_the_command_that_opened_it(tmp_path):
    route_file = tmp_path / "route.json"
    route_file.write_text(json.dumps(demo_route()))
    try:
        opened = _cli("route-review", "--route", str(route_file), "--lang", "zh")
        assert opened.returncode == 0, opened.stderr
        session = json.loads(opened.stdout)
        assert session["status"] == "pending" and session["server"] == "detached"
        assert opened.stderr.startswith("Review URL: http://127.0.0.1:")
        assert session["collect_with"].endswith(f"--session-id {session['session_id']} --wait")

        # the command has exited; the page and its API are still served
        assert "Structure review" in urlopen(session["review_url"]).read().decode()
        pending = json.loads(_cli("review-result", "--session-id", session["session_id"]).stdout)
        assert pending["status"] == "pending"
        _post(_api(session, "/complete"), {"verdict": "approved", "comment": "fine"}, session["review_token"])

        collected = json.loads(_cli("review-result", "--session-id", session["session_id"], "--wait", "--timeout-s", "5").stdout)
        assert collected["status"] == "completed" and collected["route_result"]["verdict"] == "approved"
        assert json.loads(_cli("review-result", "--session-id", "nope").stdout)["status"] == "error"

        assert json.loads(_cli("review-mode", "ask").stdout)["review_mode"] == "ask"
        asked = json.loads(_cli("review-batch", "--input", str(route_file)).stdout)
        assert asked["status"] == "ask_user_first"
    finally:
        review_daemon.stop_daemon()


def test_demo_sessions_open_with_their_planted_mistakes():
    demo = service.chem_review_demo(language="zh", detached=False)
    assert set(demo["sessions"]) == {"structures", "route"}
    structures = service.chem_structure_review_result(demo["sessions"]["structures"]["session_id"])
    assert structures["items"][1]["input_smiles"] == MISREAD_PROLINE and structures["language"] == "zh"
    assert rc.describe(MISREAD_PROLINE, reference=TRUE_PROLINE)["versus_reference"]["relation"] == "diastereomers"


# ------------------------------------------ what a first real session exposed


def test_your_own_reading_is_reviewed_beside_its_source_image():
    image = demo_items()[1]["image_path"]
    opened = service.chem_input_review(smiles=TRUE_PROLINE, source_image=image, note="OH read as a wedge", language="zh")
    item = opened["items"][0]
    assert item["image"] and item["note"] == "OH read as a wedge" and item["label"] == "Structure read from the image"
    assert "lang=zh" in opened["review_url"]


def test_tool_results_leave_drawings_out_and_accept_the_obvious_argument_names(tmp_path):
    compact = service._without_drawings(
        {"candidates": [{"isomeric_smiles": "CCO", "svg": "<svg/>", "molblock": "M  END",
                         "metadata": {"rank": 1, "atoms": [1], "bonds": [2]}}], "atoms": "kept outside metadata"}
    )
    assert compact == {"candidates": [{"isomeric_smiles": "CCO", "metadata": {"rank": 1}}], "atoms": "kept outside metadata"}
    missing = service.chem_parse_image(image_path=str(tmp_path / "missing.png"))
    assert missing["path"].endswith("missing.png") and not missing.get("candidates")

    assert service.chem_compare_structures(a="CCO", b="COC")["relation"] == "constitutional_isomers"
    route_file = FIXTURES / "roche_polyketide_fragment.json"
    assert service.chem_route_check(route_path=str(route_file))["status"] == "consistent"
    assert service.chem_route_check()["status"] == "error"
    assert service.chem_route_review(route_path=str(route_file))["kind"] == "route"
    assert service.chem_route_review(route_path=str(tmp_path / "nope.json"))["status"] == "error"


def test_esterification_in_an_alcohol_is_read_as_in_situ_acid():
    from codex_organic_chem import kit

    def audit(conditions, reactant, product):
        route = {"target": product, "compounds": {"1": reactant, "2": product},
                 "steps": [{"step": 1, "reactants": ["1"], "product": "2", "conditions": conditions}]}
        return [issue["code"] for issue in kit.route_check(route)["issues"]]

    amino_acid, methyl_ester = "OC(=O)[C@@H]1C[C@@H](O)CN1", "COC(=O)[C@@H]1C[C@@H](O)CN1"
    assert audit("SOCl2 (1.5 equiv), MeOH, 0 °C to reflux", amino_acid, methyl_ester) == []
    assert audit("NaBH4, MeOH", amino_acid, methyl_ester) == ["formula_change_not_matched_by_conditions"]
    boc_acid, boc_ester = "OC(=O)[C@@H]1C[C@@H](O)CN1C(=O)OC(C)(C)C", "COC(=O)[C@@H]1C[C@@H](O)CN1C(=O)OC(C)(C)C"
    assert audit("MeOH, SOCl2", boc_acid, boc_ester) == ["bystander_group_vs_reagent"]  # HCl would take the Boc off
    assert "bystander_group_vs_reagent" in audit("SOCl2, DMF (cat.), CH2Cl2", "OC(=O)CCCO", "ClC(=O)CCCO")


def test_cli_wait_gives_up_without_closing_the_session():
    try:
        waited = _cli("input-review", "--smiles", "CCO", "--wait", "--wait-s", "1", "--timeout-s", "60")
        result = json.loads(waited.stdout)
        assert result["status"] == "pending" and result["wait_timed_out"] is True
        assert "review-result --session-id" in result["collect_with"]
    finally:
        review_daemon.stop_daemon()


@pytest.mark.parametrize(
    "reference, edits, atoms_that_may_move",
    [
        ("CC(O)C(C)C=CC", [("set", "b5-6", "E")], 2),
        ("CC(O)C(C)C=CC", [("set", "b5-6", "Z")], 2),
        ("C/C=C/C(=O)OCc1ccccc1", [("invert", "b1-2", None)], 2),
        ("C/C=C/C=C/C(=O)O", [("invert", "b1-2", None), ("invert", "b3-4", None)], 4),
        ("CC/C(C)=C/C(=O)OC", [("set", "b2-4", "Z")], 4),
        ("C1CCCC/C=C/CCCC1", [("invert", "b5-6", None)], 11),  # ring bond: redrawn whole
        ("C/C=C/[C@H](O)CC", [("invert", "a3", None)], 0),  # a centre needs no new coordinates
    ],
)
def test_the_drawing_follows_an_e_z_edit_and_moves_as_little_as_it_can(reference, edits, atoms_that_may_move):
    from rdkit import Chem
    from codex_organic_chem import structure_compare as sc

    start = rc.describe(reference)["smiles"]
    template = rc._layout(sc.parse(start))
    edited = start
    for op, target, value in edits:
        edited = rc.apply_edit(edited, op, target, value)
    laid = rc.lay_out(sc.parse(edited), template)
    assert rc._double_bonds_drawn_right(laid)  # the picture shows the isomer the label names
    # independent of that check: a Molfile reader takes E/Z from the coordinates
    shown = Chem.MolToSmiles(Chem.MolFromMolBlock(Chem.MolToMolBlock(laid)))
    assert shown == Chem.MolToSmiles(sc.parse(edited))

    query = Chem.Mol(template)
    Chem.RemoveStereochemistry(query)
    match = laid.GetSubstructMatch(query, useChirality=False)
    old, new = template.GetConformer(), laid.GetConformer()
    moved = sum(
        1
        for i, j in enumerate(match)
        if abs(old.GetAtomPosition(i).x - new.GetAtomPosition(j).x) > 0.05
        or abs(old.GetAtomPosition(i).y - new.GetAtomPosition(j).y) > 0.05
    )
    assert moved <= atoms_that_may_move


def test_e_and_z_are_no_longer_drawn_alike():
    import re

    reference = "CC(O)C(C)C=CC"
    drawings = {}
    for wanted in ("E", "Z"):
        described = rc.describe(rc.apply_edit(reference, "set", "b5-6", wanted), reference=reference)
        drawings[wanted] = sorted(re.findall(r"class='bond-\d+ atom-\d+ atom-\d+' d='([^']+)'", described["svg"]))
    assert drawings["E"] and drawings["E"] != drawings["Z"]
