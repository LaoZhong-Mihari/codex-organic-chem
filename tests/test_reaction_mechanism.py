import re

from codex_organic_chem.service import (
    chem_mechanism_draft,
    chem_mechanism_render,
    chem_mechanism_spec_example,
    chem_reaction_analyze,
    chem_synthesis_suggest,
)


def test_reaction_sanity_check_sn2():
    result = chem_reaction_analyze("CBr.[OH-]>>CO.[Br-]", mode="sanity_check")
    assert result["reaction_smiles"] == "CBr.[OH-]>>CO.[Br-]"
    assert "possible_SN2_substitution" in result["analysis"]["reaction_classes"]
    assert result["warnings"]


def test_retro_hints_for_ester():
    result = chem_reaction_analyze("CC(=O)Oc1ccccc1C(=O)O>>CC(=O)Oc1ccccc1C(=O)O", mode="retro")
    hints = result["analysis"]["retrosynthesis_hints"]
    assert any("ester" in hint["disconnection"] for hint in hints)


def test_mechanism_draft_has_steps_and_boundaries():
    result = chem_mechanism_draft("CBr.[OH-]>>CO.[Br-]", style="stepwise")
    assert result["steps"]
    assert result["evidence_boundary"]["tool_facts"]
    assert "rule-based hypotheses" in result["warnings"][0]
    svg = result["steps"][0]["rendered_svg"]
    assert svg.index("<svg") < svg.index("<defs")


def test_publication_mechanism_requires_confirmation():
    result = chem_mechanism_draft("CBr.[OH-]>>CO.[Br-]", quality="publication")
    assert result["status"] == "awaiting_user_confirmation"
    assert result["confirmation_required"] is True


def test_publication_mechanism_confirmed_has_package():
    result = chem_mechanism_draft("CBr.[OH-]>>CO.[Br-]", quality="publication", structure_confirmed=True)
    assert result["status"] == "ok"
    assert result["publication_package"]["panels"]
    assert result["publication_package"]["publication_ready"] is False
    assert result["publication_package"]["recommended_next_tool"] == "chem_mechanism_render"
    assert result["step_explanations"][0]["why_it_is_plausible"]


def test_synthesis_suggest_is_stepwise_and_requires_confirmation():
    blocked = chem_synthesis_suggest("CC(=O)Oc1ccccc1C(=O)O", literature=False)
    assert blocked["status"] == "awaiting_user_confirmation"
    result = chem_synthesis_suggest("CC(=O)Oc1ccccc1C(=O)O", confirmed=True, literature=False)
    assert result["status"] == "awaiting_route_decision"
    assert result["options"]
    assert result["options"][0]["single_step_only"] is True


def test_mechanism_render_draws_atom_mapped_canvas(tmp_path):
    spec = chem_mechanism_spec_example()
    result = chem_mechanism_render(spec, output_dir=str(tmp_path))
    assert result["status"] == "ok"
    assert result["spec_version"] == "2.0"
    assert "HO" in result["svg"]
    assert ">H3C<" in result["svg"]
    assert ">CH3<" in result["svg"]
    assert "lp → C" not in result["svg"]
    assert "δ+" not in result["svg"]
    assert "nucleophile" not in result["svg"]
    assert "alkyl bromide" not in result["svg"]
    assert "[OH-:1].[CH3:2][Br:3]" not in str(spec["panels"][0]["molecules"])
    assert "mechanism.svg" in result["output_files"][0]
    assert (tmp_path / "mechanism.svg").exists()
    assert (tmp_path / "mechanism.trace.json").exists()
    assert (tmp_path / "mechanism.chemdoodle.json").exists()
    assert (tmp_path / "mechanism.chemdoodle.html").exists()
    assert (tmp_path / "mechanism.ketcher.json").exists()
    assert (tmp_path / "mechanism.marvin.json").exists()
    assert result["mechanism_accuracy_mode"] == "explicit_intermediate_atom_mapped_canvas"
    assert result["publication_checks"]["ready_for_human_review"] is True
    assert "editor_followup" in result
    assert result["editor_followup"]["recommended_sequence"]
    assert "<fragment" in result["cdxml"]
    assert "<n " in result["cdxml"]
    assert "<b " in result["cdxml"]
    assert "<arrow" in result["cdxml"]
    assert "class='lone-pair-object'" in result["svg"]
    assert result["svg"].count("class='lone-pair-object'") == 1
    assert "class='formal-charge-marker'" in result["svg"]
    assert "class='charge-ring'" in result["svg"]
    assert "data-routed-curvature" in result["svg"]
    assert "data-overlap-score" in result["svg"]
    arrow_paths = re.findall(r"<path class='mech-arrow' d='M ([^']+)'", result["svg"])
    assert arrow_paths
    first_arrow_numbers = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", arrow_paths[0])]
    assert first_arrow_numbers[-2] < 212.6
    assert abs(first_arrow_numbers[-2] - 212.6) > 5.0
    assert "data-source-object='p1:m1:lp:1:0'" in result["svg"]
    assert "data-target-object='p1:m2:atom:2'" in result["svg"]
    assert "data-source-object='p1:m2:bond:2-3'" in result["svg"]
    assert "Pusher" in result["chemdoodle_json_string"]
    assert result["chemdoodle_json"]["s"][0]["o1"]
    assert result["chemdoodle_json"]["s"][0]["o2"]
    assert result["chemdoodle_json"]["s"][0]["e"] == 2
    assert any(atom.get("p") == 1 for molecule in result["chemdoodle_json"]["m"] for atom in molecule["a"])
    assert sum(1 for molecule in result["chemdoodle_json"]["m"] for atom in molecule["a"] if atom.get("p")) == 1
    assert result["chemdoodle_json"]["metadata"]["atom_label_overrides"]
    assert any(atom.get("display_label") == "HO" for molecule in result["chemdoodle_json"]["m"] for atom in molecule["a"])
    assert "ChemDoodle.ViewerCanvas" in result["chemdoodle_html"]
    assert "ChemDoodle.readJSON" in result["chemdoodle_html"]
    assert result["mechanism_trace"]["steps"][0]["electron_moves"][0]["source"]["kind"] == "lone_pair"
    assert result["mechanism_trace"]["steps"][0]["electron_moves"][0]["source"]["molecule_index"] == 1
    assert result["mechanism_trace"]["steps"][0]["bond_changes"]["formed"]
    assert result["ketcher_adapter"]["adapter"] == "ketcher_input_correction"
    assert result["marvin_adapter"]["adapter"] == "marvin_review_target"


def test_mechanism_render_blocks_hidden_lone_pair_source():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["lone_pairs"] = []
    result = chem_mechanism_render(spec)
    assert result["status"] == "blocked_for_publication"
    assert any("source lone pair is not visible" in issue for issue in result["publication_checks"]["blocking_issues"])


def test_mechanism_render_warns_on_atom_center_arrow_source():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["arrows"][0].pop("from_lone_pair_atom_map")
    spec["panels"][0]["arrows"][0]["from_atom_map"] = 1
    result = chem_mechanism_render(spec)
    assert result["status"] == "ok_with_warnings"
    assert any("source is an atom center" in warning for warning in result["warnings"])


def test_mechanism_render_warns_on_disconnected_fragments_with_arrows():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["molecules"] = [{"smiles": "[OH-:1].[CH3:2][Br:3]", "label": "bad combined reactants"}]
    result = chem_mechanism_render(spec)
    assert result["status"] == "ok_with_warnings"
    assert any("contains disconnected fragments" in warning for warning in result["warnings"])


def test_mechanism_render_supports_legacy_v1_spec():
    spec = chem_mechanism_spec_example()
    spec.pop("spec_version")
    spec["panels"][0].pop("graph_edits")
    result = chem_mechanism_render(spec)
    assert result["status"] == "ok"
    assert result["spec_version"] == "1.0"
    assert result["mechanism_trace"]["steps"][0]["electron_moves"]


def test_mechanism_render_blocks_graph_edit_missing_atom_map():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["graph_edits"].append({"type": "form_bond", "atom_maps": [1, 99]})
    result = chem_mechanism_render(spec)
    assert result["status"] == "blocked_for_publication"
    assert any("graph_edit" in issue and "99" in issue for issue in result["publication_checks"]["blocking_issues"])


def test_mechanism_render_handles_single_electron_pusher():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["arrows"] = [
        {
            "from_bond": [2, 3],
            "to_atom_map": 3,
            "kind": "single_electron",
            "electron_count": 1,
            "label": "fishhook",
        }
    ]
    spec["panels"][0]["graph_edits"] = [{"type": "radical", "atom_map": 3}]
    result = chem_mechanism_render(spec)
    assert result["status"] == "ok"
    assert "arrow-radical" in result["svg"]
    assert "data-source-object='p1:m2:bond:2-3'" in result["svg"]
    assert "data-target-object='p1:m2:atom:3'" in result["svg"]
    assert "SingleElectron" in result["cdxml"]
    assert result["chemdoodle_json"]["s"][0]["e"] == 1
    assert result["chemdoodle_json"]["s"][0]["o1"] == "p1m2b0"
    assert result["chemdoodle_json"]["s"][0]["o2"] == "p1m2a1"
    assert result["mechanism_trace"]["steps"][0]["electron_moves"][0]["electron_count"] == 1


def test_mechanism_trace_records_proton_transfer_and_resonance_metadata():
    spec = chem_mechanism_spec_example()
    spec["panels"][0]["graph_edits"] = [
        {"type": "proton_transfer", "from_atom_map": 1, "to_atom_map": 2},
        {"type": "resonance", "atom_maps": [1, 2, 3]},
    ]
    spec["panels"][0]["proton_transfers"] = [{"from_atom_map": 1, "to_atom_map": 2}]
    result = chem_mechanism_render(spec)
    edits = result["mechanism_trace"]["steps"][0]["graph_edits"]
    assert {edit["type"] for edit in edits} == {"proton_transfer", "resonance"}
    assert result["mechanism_trace"]["steps"][0]["proton_counterion_handling"]["proton_transfers"]


def test_mechanism_render_preserves_atom_label_style_and_routes_arrows_to_boundary():
    spec = {
        "spec_version": "2.0",
        "journal_style": "acs",
        "presentation_mode": "publication",
        "layout": {
            "columns": 1,
            "panel_width": 360,
            "panel_height": 150,
            "show_footer": False,
            "show_notes": False,
            "show_molecule_labels": False,
        },
        "panels": [
            {
                "title": "Hydrogen label placement",
                "molecules": [
                    {
                        "smiles": "[O:1][CH2:2][Br:3]",
                        "show_carbons": True,
                    }
                ],
                "charges": [{"atom_map": 2, "label": "+", "color": "#111"}],
            }
        ],
    }
    result = chem_mechanism_render(spec)
    svg = result["svg"]
    assert result["status"] == "ok"
    assert "class='atom-hydrogen'" not in svg
    assert ">O<" in svg
    assert ">CH2<" in svg
    assert ">Br<" in svg
    assert "class='formal-charge-marker'" in svg
