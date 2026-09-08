import statistics

from codex_organic_chem.figure_audit import FigureElement, audit_figure
from codex_organic_chem.figure_style import render_molecule, resolve_preset
from codex_organic_chem.service import chem_draw, chem_mechanism_render, chem_mechanism_spec_example, chem_route_figure

MOLS = {
    "small": "CCO",
    "aromatic": "CC(=O)Oc1ccccc1C(=O)O",
    "tall": "CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(cc1)C(=O)NC1CCOC1",
    "stereo": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
}


def test_fixed_bond_length_is_uniform_across_molecule_sizes():
    preset = resolve_preset("acs")
    target = preset["bond_length_px"]
    for smiles in MOLS.values():
        depiction = render_molecule(smiles=smiles, preset=preset)
        med = statistics.median(depiction.bond_lengths_px)
        assert abs(med - target) < target * 0.05, f"{smiles}: {med} vs {target}"
        assert not depiction.warnings


def test_depiction_reports_geometry_for_audit():
    depiction = render_molecule(smiles="CC(=O)Oc1ccccc1C(=O)O", preset="acs")
    assert depiction.atom_points
    assert depiction.label_boxes  # O/OH labels exist
    assert depiction.ink_box is not None
    moved = depiction.translated_geometry(100.0, 50.0)
    assert moved["ink_box"][0] == depiction.ink_box[0] + 100.0


def test_audit_blocks_overlapping_structures():
    elements = [
        FigureElement(kind="molecule", label="a", box=(0, 0, 100, 100), bond_lengths_px=[30.0]),
        FigureElement(kind="molecule", label="b", box=(50, 0, 150, 100), bond_lengths_px=[30.0]),
    ]
    checks = audit_figure(renderer="test", width=200, height=120, elements=elements)
    assert checks["blocking_issues"]
    assert not checks["ready_for_visual_review"]


def test_audit_blocks_clipped_structures_and_flags_bond_drift():
    elements = [
        FigureElement(kind="molecule", label="clipped", box=(150, 10, 260, 90), bond_lengths_px=[30.0]),
        FigureElement(kind="molecule", label="tiny-bonds", box=(10, 10, 90, 90), bond_lengths_px=[14.0]),
    ]
    checks = audit_figure(renderer="test", width=200, height=100, elements=elements, target_bond_px=30.0)
    assert any("clipped" in issue for issue in checks["blocking_issues"])
    assert any("inconsistent" in warning for warning in checks["warnings"])


def test_audit_passes_clean_layout():
    elements = [
        FigureElement(kind="molecule", label="a", box=(10, 10, 90, 90), bond_lengths_px=[30.0, 30.0]),
        FigureElement(kind="molecule", label="b", box=(120, 10, 190, 90), bond_lengths_px=[30.0]),
    ]
    checks = audit_figure(renderer="test", width=220, height=110, elements=elements, target_bond_px=30.0)
    assert checks["blocking_issues"] == []
    assert checks["ready_for_visual_review"] is True
    assert checks["publication_ready"] is False  # human review still required


def test_route_figure_audit_runs_automatically(tmp_path):
    spec = {
        "title": "Audit integration",
        "rows": [
            {
                "items": [
                    {"id": "a", "label": "A", "smiles": "CCO"},
                    {"type": "arrow", "label": "step"},
                    {"id": "b", "label": "B", "smiles": "CC=O"},
                ]
            }
        ],
    }
    result = chem_route_figure(spec=spec, output_dir=str(tmp_path))
    checks = result["publication_checks"]
    assert checks["audit"] == "automatic_geometry_and_raster"
    assert checks["blocking_issues"] == []
    geometry = checks["bond_geometry"]
    assert geometry["min_median_bond_px"] == geometry["max_median_bond_px"]


def test_route_figure_audit_blocks_forced_undersized_canvas():
    spec = {
        "width": 260,
        "height": 160,
        "rows": [
            {
                "items": [
                    {"id": "a", "label": "A", "smiles": "CC(C)(C)OC(=O)N1CCC(CC1)c1ccccc1"},
                    {"type": "arrow", "label": "step"},
                    {"id": "b", "label": "B", "smiles": "O=C(O)c1ccccc1O"},
                ]
            }
        ],
    }
    result = chem_route_figure(spec=spec)
    assert result["status"] == "blocked_for_publication"
    assert result["publication_checks"]["blocking_issues"]


def test_mechanism_render_includes_figure_audit():
    result = chem_mechanism_render(chem_mechanism_spec_example())
    audit = result["publication_checks"]["figure_audit"]
    assert audit["renderer"] == "rdkit_compact_mechanism"
    per_element = audit["bond_geometry"]["per_element"]
    assert per_element
    medians = [entry["median_bond_px"] for entry in per_element]
    assert max(medians) - min(medians) < 2.0


def test_chem_draw_png_derives_from_svg_master(tmp_path):
    png = tmp_path / "aspirin.png"
    result = chem_draw(smiles="CC(=O)Oc1ccccc1C(=O)O", output="png", output_file=str(png))
    if result["status"] == "unavailable":  # no rasterizer on this machine
        assert result["warnings"]
        return
    assert result["status"] == "ok"
    assert png.exists()
    assert "<svg" in result["svg_master"]
    assert result["raster_check"]["status"] == "ok"


def test_chem_draw_png_without_output_file_is_an_error():
    result = chem_draw(smiles="CCO", output="png")
    assert result["status"] == "error"
