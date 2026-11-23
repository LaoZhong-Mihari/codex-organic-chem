from codex_organic_chem.service import chem_route_figure, chem_route_figure_spec_example


def test_route_figure_renders_composed_svg(tmp_path):
    spec = {
        "title": "Two-step test route",
        "subtitle": "Renderer should compose molecules, arrows, and labels.",
        "rows": [
            {
                "items": [
                    {"id": "s0", "label": "Starting material", "code": "S0", "smiles": "CCO"},
                    {"type": "arrow", "label": "oxidation", "sublabel": "monitor"},
                    {"id": "p1", "label": "Product", "code": "P1", "smiles": "CC=O"},
                ]
            }
        ],
        "footnotes": ["Review stereochemistry before final use."],
    }

    result = chem_route_figure(spec=spec, output_dir=str(tmp_path), basename="test_route")

    assert result["status"] == "ok"
    assert result["kind"] == "route_figure"
    assert "codex_route_figure_composed_svg" == result["publication_checks"]["renderer"]
    assert "<svg" in result["svg"]
    assert "Two-step test route" in result["svg"]
    assert "class=\"molLabel\"" in result["svg"]
    assert "class=\"arrowLabel\"" in result["svg"]
    assert "oxidation" in result["svg"]
    assert "Starting material" in result["svg"]
    assert (tmp_path / "test_route.svg").exists()
    assert any(path.endswith("test_route.svg") for path in result["output_files"])
    assert [mol["id"] for mol in result["molecules"]] == ["s0", "p1"]


def test_route_figure_example_spec_is_renderable():
    result = chem_route_figure(chem_route_figure_spec_example())

    assert result["status"] == "ok"
    assert "Example protected route" in result["svg"]
    assert result["molecules"]


def test_route_figure_rejects_invalid_molecule():
    spec = {"rows": [{"items": [{"id": "bad", "label": "Bad", "smiles": "not-a-smiles"}]}]}

    result = chem_route_figure(spec)

    assert result["status"] == "error"
    assert result["warnings"]


# --- scheme layout: alignment, subscripts, captions, row wrapping -------------

import pytest  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdFMCS  # noqa: E402

from codex_organic_chem.figure_style import prepare_mol  # noqa: E402
from codex_organic_chem.route_figure import _align_to_previous, render_route_figure, subscript_formulas  # noqa: E402

PROLINES = [  # three steps on one Boc-protected proline core
    ("1", "COC(=O)[C@@H]1C[C@@H](O)CN1"),
    ("2", "COC(=O)[C@@H]1C[C@@H](O)CN1C(=O)OC(C)(C)C"),
    ("3", "COC(=O)[C@@H]1C[C@H](OS(C)(=O)=O)CN1C(=O)OC(C)(C)C"),
]


def proline_route(**extra):
    items = []
    for index, (label, smiles) in enumerate(PROLINES):
        if index:
            items.append({"type": "arrow", "label": "step", "sublabel": ""})
        items.append({"label": label, "smiles": smiles})
    return {"rows": [{"items": items}], **extra}


def test_formula_digits_become_subscripts_and_nothing_else():
    assert subscript_formulas("Boc2O, Et3N") == "Boc₂O, Et₃N"
    assert subscript_formulas("H2 (1 atm), Pd/C") == "H₂ (1 atm), Pd/C"
    assert subscript_formulas("THF/H2O, 0 °C") == "THF/H₂O, 0 °C"
    assert subscript_formulas("Pd(PPh3)4, 10 mol% L1") == "Pd(PPh₃)₄, 10 mol% L1"
    for unchanged in ["2,2-dimethoxypropane", "18-crown-6", "1M HCl", "SN2", "(SN2)", "E1cb", "-78 °C, 12 h"]:
        assert subscript_formulas(unchanged) == unchanged


def test_arrow_text_is_subscripted_unless_switched_off():
    spec = {"rows": [{"items": [{"label": "1", "smiles": "CCO"}, {"type": "arrow", "label": "K2CO3", "sublabel": "CH2Cl2"},
                                {"label": "2", "smiles": "CC=O"}]}]}
    assert "K₂CO₃" in render_route_figure(spec)["svg"]
    assert "K2CO3" in render_route_figure({**spec, "subscript_formulas": False})["svg"]


def test_shared_core_keeps_its_coordinates_and_stereo():
    template = prepare_mol(smiles=PROLINES[1][1])
    mol = Chem.MolFromSmiles(PROLINES[2][1])
    aligned, info = _align_to_previous(mol, [{"mol": template, "name": "2"}])
    assert info["aligned_to"] == "2" and info["aligned_atoms"] >= 15
    query = Chem.MolFromSmarts(
        rdFMCS.FindMCS([template, aligned], ringMatchesRingOnly=True, completeRingsOnly=True).smartsString
    )
    t_conf, a_conf = template.GetConformer(), aligned.GetConformer()
    for t_idx, a_idx in zip(template.GetSubstructMatch(query), aligned.GetSubstructMatch(query)):
        if template.GetAtomWithIdx(t_idx).IsInRing():
            assert (t_conf.GetAtomPosition(t_idx) - a_conf.GetAtomPosition(a_idx)).Length() < 1e-3
    lengths = [(a_conf.GetAtomPosition(b.GetBeginAtomIdx()) - a_conf.GetAtomPosition(b.GetEndAtomIdx())).Length()
               for b in aligned.GetBonds()]
    assert max(lengths) / min(lengths) < 1.05  # one bond length throughout
    round_trip = Chem.MolFromMolBlock(Chem.MolToMolBlock(aligned))
    assert Chem.MolToSmiles(round_trip) == Chem.MolToSmiles(mol)  # the drawn wedges match the input


def test_route_figure_reports_alignment_and_can_switch_it_off():
    aligned = render_route_figure(proline_route())
    assert aligned["status"] == "ok"
    assert [m.get("aligned_to") for m in aligned["molecules"]] == [None, "1", "2"]
    plain = render_route_figure(proline_route(align_structures=False))
    assert all("aligned_to" not in m for m in plain["molecules"])


def test_unrelated_structures_are_not_forced_into_alignment():
    spec = {"rows": [{"items": [{"label": "1", "smiles": "c1ccc2ccccc2c1"}, {"type": "arrow", "label": "x"},
                                {"label": "2", "smiles": "OCC(O)CO"}]}]}
    assert all("aligned_to" not in m for m in render_route_figure(spec)["molecules"])


def test_molfile_items_keep_their_own_layout_but_guide_the_next_structure():
    spec = proline_route()
    spec["rows"][0]["items"][2] = {"label": "2", "molfile": Chem.MolToMolBlock(prepare_mol(smiles=PROLINES[1][1]))}
    molecules = render_route_figure(spec)["molecules"]
    assert "aligned_to" not in molecules[1]
    assert molecules[2]["aligned_to"] == "2"


def test_long_names_wrap_instead_of_running_off_the_canvas():
    spec = {"rows": [{"items": [{"label": "1", "code": "trans-4-hydroxy-L-proline methyl ester hydrochloride",
                                 "smiles": "CCO"}]}]}
    result = render_route_figure(spec)
    assert result["status"] == "ok"
    assert result["svg"].count('class="code"') >= 2


def test_a_long_title_widens_the_canvas():
    result = render_route_figure({"title": "A deliberately long title " * 4, "rows": [{"items": [{"label": "1", "smiles": "CCO"}]}]})
    assert result["status"] == "ok"
    assert result["layout"]["canvas"][0] > max(result["layout"]["row_widths"])


@pytest.mark.parametrize("max_width", [700, 900])
def test_long_rows_wrap_at_arrows_into_even_rows(max_width):
    items = []
    for index in range(7):
        if index:
            items.append({"type": "arrow", "label": f"step {index}"})
        items.append({"label": str(index + 1), "smiles": "c1ccccc1C(=O)N" + "C" * index})
    single = render_route_figure({"rows": [{"items": items}]})
    wrapped = render_route_figure({"rows": [{"items": items}], "max_width": max_width})
    assert single["layout"]["rows"] == 1
    widths = wrapped["layout"]["row_widths"]
    assert len(widths) >= 2 and max(widths) <= max_width
    assert min(widths) > 0.6 * max(widths)  # no row left with a single step while the others are full
    assert wrapped["status"] == "ok"
    rows_start = wrapped["svg"].count('class="rxnArrow"')
    assert rows_start == 6  # wrapping moves arrows, it never drops one
