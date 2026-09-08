"""Regressions for failures found by comparing the actual rendered figures."""

import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from codex_organic_chem.figure_style import resolve_preset
from codex_organic_chem.mechanism_scene import (
    _head,
    _panel,
    _sample,
    _shaft_points,
    add,
    inside,
    mul,
    norm,
    sub,
    unit,
)
from codex_organic_chem.service import chem_mechanism_render

FIXTURES = Path(__file__).parent / "fixtures/mechanisms"


def fixture(name="sn2"):
    return json.loads((FIXTURES / f"{name}.json").read_text())


@pytest.mark.parametrize(
    "name",
    [
        "sn2",
        "proton_transfer",
        "e2",
        "carbonyl_addition",
        "homolysis",
        "cyclohexanone",
        "sn2_skeletal",
    ],
)
def test_compact_examples_have_local_arrows_and_shared_bond_scale(name):
    r = chem_mechanism_render(fixture(name))
    audit = r["publication_checks"]["figure_audit"]
    assert not audit["blocking_issues"]
    assert audit["ready_for_visual_review"]
    assert not audit["publication_ready"]
    assert audit["bond_geometry"]["target_bond_px"] == 30
    assert abs(audit["bond_geometry"]["min_median_bond_px"] - 30) < 0.2
    assert abs(audit["bond_geometry"]["max_median_bond_px"] - 30) < 0.2
    assert all(
        a["length_bonds"] < 2.4 and a["ink_contacts"] == 0
        for a in audit["electron_arrows"]
    )


def test_canvas_width_no_longer_stretches_molecules_and_electron_arrows():
    s = fixture()
    original = chem_mechanism_render(s)["publication_checks"]["figure_audit"]
    s["layout"]["panel_width"] = 3000
    wider = chem_mechanism_render(s)["publication_checks"]["figure_audit"]
    assert wider["canvas"] == original["canvas"]
    assert wider["electron_arrows"] == original["electron_arrows"]


def test_long_arrows_are_blocking_instead_of_merely_warning():
    s = fixture()
    s["panels"][0]["composition"] = {"species_gap_bonds": 8, "auto_compact": False}
    r = chem_mechanism_render(s)
    assert r["status"] == "blocked_for_publication"
    audit = r["publication_checks"]["figure_audit"]
    assert not audit["ready_for_visual_review"]
    assert any("bond lengths long" in e for e in audit["blocking_issues"])


def test_arrow_through_plus_sign_is_blocking():
    s = fixture()
    s["panels"][0]["composition"] = {"auto_compact": False}
    s["panels"][0]["arrows"][0]["routing"] = {"lift_bonds": 0, "lock_side": True}
    r = chem_mechanism_render(s)
    assert r["status"] == "blocked_for_publication"
    assert any("operator ink" in e for e in r["publication_checks"]["blocking_issues"])


def test_hydroxide_lone_pair_stays_at_oxygen_not_the_implicit_hydrogen():
    panel = _panel(fixture("sn2_skeletal")["panels"][0], resolve_preset("acs"))
    attack = panel.curves[0]
    oxygen = attack.source.point
    assert attack.lp is not None
    assert norm(sub(attack.lp, oxygen)) < 12
    assert all(
        not inside(attack.lp, b, 1) for b in attack.source.molecule.label_boxes()
    )


def test_invalid_unmapped_input_can_still_return_a_blocked_preview():
    s = {"spec_version": "2.1", "panels": [{"molecules": [{"smiles": "CCO"}]}]}
    r = chem_mechanism_render(s)
    assert r["status"] == "blocked_for_publication"
    assert "data-painter='rdkit-shared'" in r["svg"]


@pytest.mark.parametrize("angle", [0, 30, 90, 137, 180, 265])
def test_open_v_is_bisected_by_terminal_tangent_at_every_orientation(angle):
    tip = (37.5, -21.0)
    direction = (math.cos(math.radians(angle)), math.sin(math.radians(angle)))
    control = sub(tip, mul(direction, 13))
    left, actual_tip, right = _head(tip, control, 2)
    assert actual_tip == tip
    midpoint = mul(add(left, right), 0.5)
    assert unit(sub(tip, midpoint)) == pytest.approx(direction)
    assert norm(sub(left, tip)) == pytest.approx(norm(sub(right, tip)))
    a, b = unit(sub(left, tip)), unit(sub(right, tip))
    assert unit(add(a, b)) == pytest.approx(mul(direction, -1))
    assert norm(sub(midpoint, tip)) == pytest.approx(4.8)
    assert norm(sub(left, right)) == pytest.approx(5.0)


@pytest.mark.parametrize("name", [p.stem for p in sorted(FIXTURES.glob("*.json"))])
def test_visible_v_join_is_centered_smooth_and_used_by_collision_audit(name):
    for item in fixture(name)["panels"]:
        panel = _panel(item, resolve_preset("acs"))
        for curve in panel.curves:
            if curve.electrons != 2:
                continue
            shaft = _shaft_points(curve)
            left, tip, right = curve.head
            backward = unit(add(unit(sub(left, tip)), unit(sub(right, tip))))
            assert shaft[0] == curve.points[0]
            assert tip == curve.points[-1]
            assert unit(sub(shaft[1], shaft[0])) == pytest.approx(
                unit(sub(curve.points[1], curve.points[0]))
            )
            assert unit(sub(shaft[-1], tip)) == pytest.approx(backward)
            # The cubic's final tangent, straight finish and V bisector coincide.
            assert unit(sub(shaft[-1], shaft[-2])) == pytest.approx(mul(backward, -1))
            body = _sample(shaft)
            assert curve.samples[: len(body)] == body
            assert curve.samples[-1] == tip
            assert curve.length == pytest.approx(
                sum(norm(sub(a, b)) for a, b in zip(curve.samples, curve.samples[1:]))
            )


@pytest.mark.parametrize("name", [p.stem for p in sorted(FIXTURES.glob("*.json"))])
def test_exported_svg_shaft_actually_reaches_the_v_tip_along_its_bisector(name):
    root = ET.fromstring(chem_mechanism_render(fixture(name))["svg"])
    shafts = [n for n in root.iter() if n.get("class") == "mech-arrow"]
    heads = [n for n in root.iter() if "electron-head" in n.get("class", "").split()]
    assert len(shafts) == len(heads)
    for shaft, head in zip(shafts, heads):
        if shaft.get("data-electrons") != "2":
            continue
        numbers = list(map(float, re.findall(r"-?\d+(?:\.\d+)?", shaft.get("d"))))
        assert len(numbers) == 10  # Cubic body plus the visible bisecting finish.
        left, tip, right = [
            tuple(map(float, p.split(","))) for p in head.get("points").split()
        ]
        control, join, end = (
            tuple(numbers[4:6]),
            tuple(numbers[6:8]),
            tuple(numbers[8:10]),
        )
        assert end == tip
        bisector = unit(add(unit(sub(left, tip)), unit(sub(right, tip))))
        assert unit(sub(join, tip)) == pytest.approx(bisector, abs=0.005)
        assert unit(sub(join, control)) == pytest.approx(mul(bisector, -1), abs=0.015)


@pytest.mark.parametrize(
    "name,head_class,vertices",
    [("sn2", "arrow-pair", 3), ("homolysis", "arrow-radical", 2)],
)
def test_pair_heads_are_rounded_open_vs_and_fishhooks_remain_single_wings(
    name, head_class, vertices
):
    svg = ET.fromstring(chem_mechanism_render(fixture(name))["svg"])
    heads = [node for node in svg.iter() if head_class in node.get("class", "").split()]
    assert heads
    assert all(node.tag == "{http://www.w3.org/2000/svg}polyline" for node in heads)
    assert all(len(node.get("points").split()) == vertices for node in heads)
    css = svg.find("{http://www.w3.org/2000/svg}style").text
    assert ".mech-arrow,.electron-head{fill:none" in css
    assert "stroke-linecap:round;stroke-linejoin:round" in css
