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
