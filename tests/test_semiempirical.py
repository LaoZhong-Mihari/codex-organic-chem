"""Tests for the GFN2-xTB / CREST records.

The previous CREST test only asserted that ``status`` was one of a set of
strings, which a stub that never ran CREST passed happily for months. So the
integration tests here skip when a binary is missing, but when it is present
they assert on the numbers: a record that runs nothing can no longer pass.
"""

from __future__ import annotations

import textwrap

import pytest

from codex_organic_chem.semiempirical import (
    _boltzmann_populations,
    _parse_ensemble,
    parse_fukui_block,
)
from codex_organic_chem.service import chem_compute
from codex_organic_chem.utils import executable_status


def _record(result, method):
    return next(item for item in result["records"] if item["method"] == method)


def _have(binary: str, flag: str = "--version") -> bool:
    return executable_status(binary, (flag,)).status == "available"


xtb_required = pytest.mark.skipif(not _have("xtb"), reason="xTB is not installed")
crest_required = pytest.mark.skipif(not _have("crest"), reason="CREST is not installed")


ATOMS = [
    {"atom_index": 0, "symbol": "C", "is_heavy": True},
    {"atom_index": 1, "symbol": "O", "is_heavy": True},
    {"atom_index": 2, "symbol": "H", "is_heavy": False},
]

FUKUI_STDOUT = textwrap.dedent(
    """
     #        f(+)     f(-)     f(0)
    Fukui functions:
     1C      0.100    0.020    0.060
     2O      0.300    0.400    0.350
     3H      0.010    0.005    0.007

    other output
    """
)


def test_fukui_rows_map_to_rdkit_atom_indices():
    rows, warnings = parse_fukui_block(FUKUI_STDOUT, ATOMS)
    assert warnings == []
    assert [row["atom_index"] for row in rows] == [0, 1, 2]
    oxygen = rows[1]
    assert oxygen["symbol"] == "O"
    assert oxygen["f_plus_electrophilic_site"] == 0.300
    assert oxygen["f_minus_nucleophilic_site"] == 0.400
    assert oxygen["f_zero_radical_site"] == 0.350


def test_fukui_rows_are_dropped_when_labels_disagree_with_the_molecule():
    """A silent index offset would attribute reactivity to the wrong atom."""
    shifted = FUKUI_STDOUT.replace(" 2O ", " 2N ")
    rows, warnings = parse_fukui_block(shifted, ATOMS)
    assert rows == []
    assert any("labelled N but atom 1 is O" in message for message in warnings)


def test_missing_fukui_table_is_reported_not_silently_empty():
    rows, warnings = parse_fukui_block("no table here", ATOMS)
    assert rows == []
    assert warnings == ["xTB produced no Fukui table."]


def test_ensemble_energies_come_from_each_frame_comment_line(tmp_path):
    path = tmp_path / "crest_conformers.xyz"
    frames = []
    for energy in (-10.000000, -9.998400, -9.990000):
        frames.append(f"2\n{energy:>18.8f}\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0")
    path.write_text("\n".join(frames) + "\n", encoding="utf-8")
    assert _parse_ensemble(path) == [-10.0, -9.9984, -9.99]


def test_missing_ensemble_file_yields_no_energies(tmp_path):
    assert _parse_ensemble(tmp_path / "absent.xyz") == []


def test_boltzmann_populations_rank_and_sum_to_one():
    populations = _boltzmann_populations([-10.0, -9.9984, -9.99], 298.15)
    assert populations[0]["relative_energy_kcal_mol"] == 0.0
    assert populations[1]["relative_energy_kcal_mol"] == pytest.approx(1.004, abs=0.01)
    assert populations[0]["population"] > populations[1]["population"] > populations[2]["population"]
    assert sum(item["population"] for item in populations) == pytest.approx(1.0, abs=0.01)


def test_compute_rejects_unknown_task():
    rec = _record(chem_compute("CCO", tasks=["quantum_astrology"]), "quantum_astrology")
    assert rec["status"] == "unavailable"


@xtb_required
def test_xtb_opt_returns_energy_from_json_not_stdout_scraping():
    rec = _record(chem_compute("CCO", tasks=["xtb_opt"]), "xtb_opt")
    assert rec["status"] == "ok"
    results = rec["results"]
    assert results["total_energy_hartree"] < 0
    assert results["homo_lumo_gap_ev"] > 0
    assert results["optimized_xyz"].splitlines()[0].strip() == "9"
    assert {entry["atom_index"] for entry in results["partial_charges"]} == set(range(9))


@xtb_required
def test_xtb_passes_formal_charge_so_ions_are_not_computed_as_neutral():
    """Hydroxide is -4.68 Eh; without --chrg xTB returns -4.43 Eh for neutral OH."""
    rec = _record(chem_compute("[OH-]", tasks=["xtb_opt"]), "xtb_opt")
    assert rec["status"] == "ok"
    assert rec["parameters"]["net_charge"] == -1
    assert rec["results"]["total_energy_hartree"] == pytest.approx(-4.68, abs=0.05)


@xtb_required
def test_reactivity_record_names_the_carbonyl_carbon_as_the_electrophile():
    rec = _record(chem_compute("CC(=O)C", tasks=["xtb_reactivity"]), "xtb_reactivity")
    assert rec["status"] == "ok"
    results = rec["results"]
    carbons = [
        site for site in results["most_electrophilic_sites"] if site["symbol"] == "C"
    ]
    # Acetone SMILES CC(=O)C: index 1 is the carbonyl carbon, 0 and 3 the methyls.
    assert carbons[0]["atom_index"] == 1
    nucleophiles = results["most_nucleophilic_sites"]
    assert nucleophiles[0]["symbol"] == "O"
    assert results["global_reactivity"]["ionization_potential_ev"] > 0
    carbonyl = next(entry for entry in results["per_atom"] if entry["atom_index"] == 1)
    assert carbonyl["partial_charge"] > 0
    assert "f_plus_electrophilic_site" in carbonyl


@xtb_required
def test_thermo_record_reports_a_true_minimum_for_ethane():
    rec = _record(chem_compute("CC", tasks=["xtb_thermo"]), "xtb_thermo")
    assert rec["status"] == "ok"
    results = rec["results"]
    assert results["imaginary_frequency_count"] == 0
    assert results["zero_point_energy_hartree"] > 0
    # Free energy includes G(RRHO), so it sits above the electronic energy.
    assert results["total_free_energy_hartree"] > results["total_energy_hartree"]


@xtb_required
def test_solvent_reaches_xtb_and_changes_the_energy():
    gas = _record(chem_compute("CC(=O)C", tasks=["xtb_opt"]), "xtb_opt")
    water = _record(chem_compute("CC(=O)C", tasks=["xtb_opt"], solvent="water"), "xtb_opt")
    assert water["parameters"]["solvent"] == "water"
    assert water["results"]["total_energy_hartree"] != gas["results"]["total_energy_hartree"]


@crest_required
def test_crest_actually_runs_and_returns_an_ensemble():
    rec = _record(chem_compute("CCO", tasks=["crest"]), "crest_conformer_search")
    assert rec["status"] == "ok", rec["warnings"]
    results = rec["results"]
    assert results["conformer_count"] >= 1
    assert results["lowest_energy_hartree"] < 0
    assert results["conformers"][0]["relative_energy_kcal_mol"] == 0.0
    assert 0 < results["dominant_population"] <= 1.0
    assert results["best_conformer_xyz"]
    # The old stub reported "installed, run it yourself" and never executed.
    assert "execution_policy" not in results
