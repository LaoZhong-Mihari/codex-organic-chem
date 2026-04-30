from codex_organic_chem.service import chem_figure_tool_status, chem_tool_doctor


def test_figure_tool_status_includes_chemdraw_and_open_source_editor():
    result = chem_figure_tool_status()
    assert "chemdraw" in result["tools"]
    assert "ketcher" in result["tools"]
    assert result["tools"]["chemdraw"]["category"] == "chemically_aware_editor"
    assert result["tools"]["ketcher"]["category"] == "open_source_structure_editor"
    assert "chemdoodle_web" in result["tools"]
    assert "marvin_js" in result["tools"]
    assert "object_bound_renderer_available" in result["summary"]
    assert result["boundaries"]


def test_doctor_includes_figure_tool_statuses():
    result = chem_tool_doctor()
    assert "publication_figure_guide" in result
    assert "figure_tool_statuses" in result
    assert "chemdraw" in result["figure_tool_statuses"]["tools"]
