"""GFN2-xTB and CREST records that return per-atom, structure-mapped numbers.

The point of this module is not "more calculations": it is to hand the caller
quantitative values keyed to RDKit atom indices, so a claim about which site
reacts can be checked against a number instead of a remembered pattern.

Three facts drive the implementation:

* xTB writes ``xtbout.json`` with ``--json``; scraping the human-readable
  stdout for energies is unnecessary and fragile.
* ``Chem.AddHs`` preserves heavy-atom indices and ``MolToXYZBlock`` writes atoms
  in mol order, so xTB's per-atom output maps back to RDKit indices directly.
* Formal charge must be passed explicitly. Without ``--chrg`` xTB treats
  hydroxide as neutral and returns an energy for the wrong species
  (-4.4284 Eh vs -4.6814 Eh).
"""

from __future__ import annotations

import json
import math
import re
import tempfile
from pathlib import Path
from typing import Any

from .models import CalculationRecord
from .rdkit_tools import RDKIT_AVAILABLE, parse_molecule
from .utils import executable_status, run_command

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through service tests
    from rdkit import Chem
    from rdkit.Chem import AllChem
else:  # pragma: no cover
    Chem = None
    AllChem = None


HARTREE_TO_KCAL = 627.5094740631
GAS_CONSTANT_KCAL = 1.987204258640832e-3  # kcal/(mol*K)
ROOM_TEMPERATURE_K = 298.15


class GeometryError(RuntimeError):
    """Raised when no 3D geometry can be produced for a semiempirical run."""


def _embed_geometry(smiles: str) -> tuple[Any, str, int, int, list[str]]:
    """Build a 3D geometry plus the electronic bookkeeping xTB needs.

    Returns (mol_with_hydrogens, xyz_block, net_charge, unpaired_electrons,
    warnings). Atom order in the xyz block matches ``mol_with_hydrogens``.
    """
    if not RDKIT_AVAILABLE:
        raise GeometryError("RDKit is required to build semiempirical input coordinates.")
    parsed = parse_molecule(smiles=smiles)
    if parsed.mol is None:
        raise GeometryError(parsed.warnings[0] if parsed.warnings else f"Invalid SMILES: {smiles}")
    warnings = list(parsed.warnings)
    mol = Chem.AddHs(parsed.mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 48879
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise GeometryError("RDKit failed to generate a 3D conformer.")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            warnings.append("Force-field pre-optimization failed; xTB started from the raw embedding.")
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    unpaired = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    return mol, Chem.MolToXYZBlock(mol), charge, unpaired, warnings


def _atom_table(mol: Any) -> list[dict[str, Any]]:
    return [
        {
            "atom_index": atom.GetIdx(),
            "symbol": atom.GetSymbol(),
            "is_heavy": atom.GetAtomicNum() > 1,
        }
        for atom in mol.GetAtoms()
    ]


def _base_args(charge: int, unpaired: int, solvent: str | None) -> list[str]:
    args = ["--chrg", str(charge)]
    if unpaired:
        args += ["--uhf", str(unpaired)]
    if solvent:
        args += ["--alpb", solvent]
    return args


def _unavailable(method: str, params: dict, message: str) -> CalculationRecord:
    return CalculationRecord(method=method, parameters=params, status="unavailable", warnings=[message])


def _read_xtb_json(work_dir: Path) -> dict[str, Any]:
    path = work_dir / "xtbout.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        return {}


FUKUI_ROW = re.compile(r"^\s*(\d+)([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")


def parse_fukui_block(stdout: str, atoms: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str]]:
    """Map xTB's Fukui table onto RDKit atom indices.

    xTB prints ``1C`` style row labels in input order, which is mol order, so
    row N belongs to atom index N-1. Row symbols are checked against the mol
    rather than trusted, because a silent offset would misattribute reactivity
    to the wrong atom.
    """
    lines = stdout.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip().startswith("Fukui functions:"))
    except StopIteration:
        return [], ["xTB produced no Fukui table."]
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    for line in lines[start + 1 :]:
        match = FUKUI_ROW.match(line)
        if not match:
            if rows:
                break
            continue
        position = int(match.group(1))
        symbol = match.group(2)
        index = position - 1
        if index >= len(atoms):
            warnings.append(f"xTB reported Fukui row {position} beyond the molecule's atom count.")
            break
        if atoms[index]["symbol"].upper() != symbol.upper():
            warnings.append(
                f"Fukui row {position} is labelled {symbol} but atom {index} is "
                f"{atoms[index]['symbol']}; per-atom values were dropped rather than misattributed."
            )
            return [], warnings
        rows.append(
            {
                "atom_index": index,
                "symbol": atoms[index]["symbol"],
                # f(+) tracks added electron density: susceptibility to
                # nucleophilic attack, i.e. the electrophilic centre.
                "f_plus_electrophilic_site": float(match.group(3)),
                # f(-) tracks removed electron density: susceptibility to
                # electrophilic attack, i.e. the nucleophilic centre.
                "f_minus_nucleophilic_site": float(match.group(4)),
                "f_zero_radical_site": float(match.group(5)),
            }
        )
    return rows, warnings


def _rank_sites(rows: list[dict[str, Any]], atoms: list[dict[str, Any]], key: str, limit: int = 5) -> list[dict]:
    heavy = [row for row in rows if atoms[row["atom_index"]]["is_heavy"]]
    ranked = sorted(heavy, key=lambda row: row[key], reverse=True)[:limit]
    return [{"atom_index": row["atom_index"], "symbol": row["symbol"], "value": round(row[key], 4)} for row in ranked]


def _parse_labelled_float(text: str, label: str) -> float | None:
    for line in text.splitlines():
        if label in line:
            for token in reversed(line.replace(":", " ").split()):
                try:
                    return float(token)
                except ValueError:
                    continue
    return None


def xtb_opt_record(
    smiles: str,
    solvent: str | None = None,
    timeout_s: int = 300,
) -> CalculationRecord:
    """Optimize the geometry and report energies read from ``xtbout.json``."""
    method = "xtb_opt"
    params = {"smiles": smiles, "task": "opt", "method": "GFN2-xTB", "solvent": solvent}
    status = executable_status("xtb", ("--version",))
    if status.status != "available":
        return _unavailable(method, params, "xTB binary is not installed or not on PATH.")
    try:
        mol, xyz, charge, unpaired, warnings = _embed_geometry(smiles)
    except GeometryError as exc:
        return CalculationRecord(method=method, parameters=params, status="error", warnings=[str(exc)])
    params.update({"net_charge": charge, "unpaired_electrons": unpaired})

    with tempfile.TemporaryDirectory(prefix="codex-chem-xtb-opt-") as tmp:
        work = Path(tmp)
        (work / "input.xyz").write_text(xyz, encoding="utf-8")
        try:
            code, stdout, stderr = run_command(
                ["xtb", "input.xyz", "--opt", "--gfn", "2", "--json", *_base_args(charge, unpaired, solvent)],
                timeout_s=timeout_s,
                cwd=str(work),
            )
        except Exception as exc:
            return CalculationRecord(
                method=method,
                parameters=params,
                status="error",
                warnings=[*warnings, f"xTB execution failed: {exc}"],
            )
        payload = _read_xtb_json(work)
        out_xyz = work / "xtbopt.xyz"
        optimized = out_xyz.read_text(encoding="utf-8", errors="replace") if out_xyz.exists() else None

    charges = payload.get("partial charges") or []
    atoms = _atom_table(mol)
    results = {
        "exit_code": code,
        "total_energy_hartree": payload.get("total energy"),
        "electronic_energy_hartree": payload.get("electronic energy"),
        "homo_lumo_gap_ev": payload.get("HOMO-LUMO gap / eV"),
        "dipole_au": payload.get("dipole / a.u."),
        "optimized_xyz": optimized,
        "index_convention": "partial_charges are indexed by RDKit atom index of the hydrogen-added molecule.",
    }
    if charges and len(charges) == len(atoms):
        results["partial_charges"] = [
            {"atom_index": atom["atom_index"], "symbol": atom["symbol"], "charge": round(float(value), 4)}
            for atom, value in zip(atoms, charges)
        ]
    if optimized is None:
        warnings.append("xTB produced no xtbopt.xyz; the optimization did not converge.")
    if stderr.strip() and code != 0:
        warnings.append(f"xTB stderr: {stderr.strip().splitlines()[-1]}")
    return CalculationRecord(
        method=method,
        tool_version=status.version,
        parameters=params,
        output_files=["xtbopt.xyz"] if optimized else [],
        results=results,
        status="ok" if code == 0 else "error",
        warnings=warnings,
    )


def xtb_reactivity_record(
    smiles: str,
    solvent: str | None = None,
    timeout_s: int = 300,
) -> CalculationRecord:
    """Per-atom electronic reactivity: Fukui indices, partial charges, IP/EA.

    This is the record that makes "which site reacts" checkable: every value
    carries the RDKit atom index it belongs to.
    """
    method = "xtb_reactivity"
    params = {"smiles": smiles, "method": "GFN2-xTB", "solvent": solvent}
    status = executable_status("xtb", ("--version",))
    if status.status != "available":
        return _unavailable(method, params, "xTB binary is not installed or not on PATH.")
    try:
        mol, xyz, charge, unpaired, warnings = _embed_geometry(smiles)
    except GeometryError as exc:
        return CalculationRecord(method=method, parameters=params, status="error", warnings=[str(exc)])
    atoms = _atom_table(mol)
    params.update({"net_charge": charge, "unpaired_electrons": unpaired})

    with tempfile.TemporaryDirectory(prefix="codex-chem-xtb-fukui-") as tmp:
        work = Path(tmp)
        xyz_path = work / "input.xyz"
        xyz_path.write_text(xyz, encoding="utf-8")
        shared = _base_args(charge, unpaired, solvent)
        try:
            code, stdout, stderr = run_command(
                ["xtb", "input.xyz", "--vfukui", "--gfn", "2", "--json", *shared],
                timeout_s=timeout_s,
                cwd=str(work),
            )
        except Exception as exc:
            return CalculationRecord(
                method=method,
                parameters=params,
                status="error",
                warnings=[*warnings, f"xTB execution failed: {exc}"],
            )
        payload = _read_xtb_json(work)
        fukui_rows, fukui_warnings = parse_fukui_block(stdout, atoms)
        warnings.extend(fukui_warnings)

        # --vomega runs GFN1 regardless of --gfn 2, so it is a separate record
        # field with its own method label rather than mixed into the GFN2 block.
        global_indices: dict[str, Any] = {}
        try:
            _, omega_out, _ = run_command(
                ["xtb", "input.xyz", "--vomega", *shared],
                timeout_s=timeout_s,
                cwd=str(work),
            )
            global_indices = {
                "method": "GFN1-xTB (xtb --vomega ignores --gfn 2)",
                "ionization_potential_ev": _parse_labelled_float(omega_out, "delta SCC IP"),
                "electron_affinity_ev": _parse_labelled_float(omega_out, "delta SCC EA"),
                "global_electrophilicity_index_ev": _parse_labelled_float(
                    omega_out, "Global electrophilicity index"
                ),
            }
        except Exception as exc:
            warnings.append(f"Global reactivity indices unavailable: {exc}")

    charges = payload.get("partial charges") or []
    if charges and len(charges) != len(atoms):
        warnings.append("xTB partial-charge count does not match the molecule; charges were dropped.")
        charges = []
    per_atom = []
    fukui_by_index = {row["atom_index"]: row for row in fukui_rows}
    for atom in atoms:
        entry = {"atom_index": atom["atom_index"], "symbol": atom["symbol"]}
        if charges:
            entry["partial_charge"] = round(float(charges[atom["atom_index"]]), 4)
        row = fukui_by_index.get(atom["atom_index"])
        if row:
            entry.update(
                {
                    "f_plus_electrophilic_site": row["f_plus_electrophilic_site"],
                    "f_minus_nucleophilic_site": row["f_minus_nucleophilic_site"],
                    "f_zero_radical_site": row["f_zero_radical_site"],
                }
            )
        per_atom.append(entry)

    results = {
        "exit_code": code,
        "total_energy_hartree": payload.get("total energy"),
        "homo_lumo_gap_ev": payload.get("HOMO-LUMO gap / eV"),
        "dipole_au": payload.get("dipole / a.u."),
        "global_reactivity": global_indices,
        "per_atom": per_atom,
        "most_electrophilic_sites": _rank_sites(fukui_rows, atoms, "f_plus_electrophilic_site"),
        "most_nucleophilic_sites": _rank_sites(fukui_rows, atoms, "f_minus_nucleophilic_site"),
        "most_radical_sites": _rank_sites(fukui_rows, atoms, "f_zero_radical_site"),
        "index_convention": (
            "atom_index refers to RDKit atom indices of the hydrogen-added molecule; "
            "heavy-atom indices match the input SMILES ordering. Site rankings list heavy atoms only."
        ),
    }
    warnings.append(
        "Fukui indices and partial charges rank sites on a single force-field geometry; "
        "they are screening evidence, not a substitute for transition-state or experimental data."
    )
    if stderr.strip() and code != 0:
        warnings.append(f"xTB stderr: {stderr.strip().splitlines()[-1]}")
    return CalculationRecord(
        method=method,
        tool_version=status.version,
        parameters=params,
        results=results,
        status="ok" if code == 0 and per_atom else "error",
        warnings=warnings,
    )


def xtb_thermo_record(
    smiles: str,
    solvent: str | None = None,
    timeout_s: int = 600,
) -> CalculationRecord:
    """Optimize, then run a Hessian for free energy and imaginary-mode count.

    ``--ohess`` would do both in one call but crashes with a Fortran format
    error in xtb 6.7.1, so the two steps are run separately.
    """
    method = "xtb_thermo"
    params = {"smiles": smiles, "method": "GFN2-xTB", "solvent": solvent}
    status = executable_status("xtb", ("--version",))
    if status.status != "available":
        return _unavailable(method, params, "xTB binary is not installed or not on PATH.")
    try:
        _, xyz, charge, unpaired, warnings = _embed_geometry(smiles)
    except GeometryError as exc:
        return CalculationRecord(method=method, parameters=params, status="error", warnings=[str(exc)])
    params.update({"net_charge": charge, "unpaired_electrons": unpaired})

    with tempfile.TemporaryDirectory(prefix="codex-chem-xtb-thermo-") as tmp:
        work = Path(tmp)
        (work / "input.xyz").write_text(xyz, encoding="utf-8")
        shared = _base_args(charge, unpaired, solvent)
        try:
            opt_code, _, _ = run_command(
                ["xtb", "input.xyz", "--opt", "--gfn", "2", *shared],
                timeout_s=timeout_s,
                cwd=str(work),
            )
            optimized = work / "xtbopt.xyz"
            geometry_arg = "xtbopt.xyz" if optimized.exists() else "input.xyz"
            if not optimized.exists():
                warnings.append("Geometry optimization produced no xtbopt.xyz; the Hessian used the input geometry.")
            hess_code, hess_out, hess_err = run_command(
                ["xtb", geometry_arg, "--hess", "--gfn", "2", "--json", *shared],
                timeout_s=timeout_s,
                cwd=str(work),
            )
        except Exception as exc:
            return CalculationRecord(
                method=method,
                parameters=params,
                status="error",
                warnings=[*warnings, f"xTB execution failed: {exc}"],
            )
        payload = _read_xtb_json(work)
        optimized_xyz = (work / "xtbopt.xyz").read_text(encoding="utf-8", errors="replace") if (work / "xtbopt.xyz").exists() else None

    imaginary = _parse_labelled_float(hess_out, "# imaginary freq.")
    imaginary_count = int(imaginary) if imaginary is not None else None
    results = {
        "opt_exit_code": opt_code,
        "hessian_exit_code": hess_code,
        "total_energy_hartree": payload.get("total energy"),
        "total_free_energy_hartree": _parse_labelled_float(hess_out, "total free energy"),
        "zero_point_energy_hartree": _parse_labelled_float(hess_out, "zero point energy"),
        "g_rrho_contribution_hartree": _parse_labelled_float(hess_out, "G(RRHO) contrib."),
        "homo_lumo_gap_ev": payload.get("HOMO-LUMO gap / eV"),
        "imaginary_frequency_count": imaginary_count,
        "optimized_xyz": optimized_xyz,
    }
    if imaginary_count:
        warnings.append(
            f"{imaginary_count} imaginary frequency/frequencies remain: this geometry is not a minimum, "
            "so the reported free energy does not describe a stable species."
        )
    warnings.append("GFN2-xTB thermochemistry is semiempirical; use it for ranking, not for reported energetics.")
    if hess_err.strip() and hess_code != 0:
        warnings.append(f"xTB stderr: {hess_err.strip().splitlines()[-1]}")
    return CalculationRecord(
        method=method,
        tool_version=status.version,
        parameters=params,
        results=results,
        status="ok" if hess_code == 0 else "error",
        warnings=warnings,
    )


def _boltzmann_populations(energies_hartree: list[float], temperature_k: float) -> list[dict[str, float]]:
    if not energies_hartree:
        return []
    lowest = min(energies_hartree)
    relative = [(value - lowest) * HARTREE_TO_KCAL for value in energies_hartree]
    weights = [math.exp(-delta / (GAS_CONSTANT_KCAL * temperature_k)) for delta in relative]
    total = sum(weights) or 1.0
    return [
        {
            "conformer": index + 1,
            "energy_hartree": energies_hartree[index],
            "relative_energy_kcal_mol": round(relative[index], 3),
            "population": round(weights[index] / total, 4),
        }
        for index in range(len(energies_hartree))
    ]


def _parse_ensemble(path: Path) -> list[float]:
    """Read absolute energies from the comment line of each xyz frame."""
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    energies: list[float] = []
    cursor = 0
    while cursor < len(lines):
        header = lines[cursor].strip()
        if not header:
            cursor += 1
            continue
        try:
            count = int(header)
        except ValueError:
            break
        if cursor + 1 >= len(lines):
            break
        try:
            energies.append(float(lines[cursor + 1].split()[0]))
        except (ValueError, IndexError):
            pass
        cursor += count + 2
    return energies


def crest_conformer_record(
    smiles: str,
    solvent: str | None = None,
    timeout_s: int = 600,
    threads: int = 4,
    quick: bool = True,
) -> CalculationRecord:
    """Run a real CREST ensemble search and report populations.

    Previously this reported "installed, run it yourself"; a capability the
    caller cannot invoke is indistinguishable from a missing one.
    """
    method = "crest_conformer_search"
    params = {
        "smiles": smiles,
        "method": "CREST/GFN2-xTB",
        "solvent": solvent,
        "timeout_s": timeout_s,
        "quick": quick,
    }
    status = executable_status("crest", ("--version",))
    if status.status != "available":
        return _unavailable(method, params, "CREST binary is not installed or not on PATH.")
    try:
        mol, xyz, charge, unpaired, warnings = _embed_geometry(smiles)
    except GeometryError as exc:
        return CalculationRecord(method=method, parameters=params, status="error", warnings=[str(exc)])
    heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    params.update({"net_charge": charge, "unpaired_electrons": unpaired, "heavy_atom_count": heavy})
    if heavy > 30:
        warnings.append(
            f"{heavy} heavy atoms: CREST cost grows steeply with size, so this run may hit the timeout. "
            "Raise timeout_s deliberately or screen with RDKit conformers first."
        )

    with tempfile.TemporaryDirectory(prefix="codex-chem-crest-") as tmp:
        work = Path(tmp)
        (work / "input.xyz").write_text(xyz, encoding="utf-8")
        args = ["crest", "input.xyz", "--gfn2", "-T", str(int(threads)), "--chrg", str(charge)]
        if unpaired:
            args += ["--uhf", str(unpaired)]
        if solvent:
            args += ["--alpb", solvent]
        if quick:
            args.append("--quick")
        try:
            code, stdout, stderr = run_command(args, timeout_s=timeout_s, cwd=str(work))
        except Exception as exc:
            return CalculationRecord(
                method=method,
                parameters=params,
                status="timeout" if "timed out" in str(exc).lower() or "timeout" in str(exc).lower() else "error",
                warnings=[
                    *warnings,
                    f"CREST did not finish within {timeout_s}s: {exc}. "
                    "Increase timeout_s, reduce the molecule, or use RDKit conformers for a quick estimate.",
                ],
            )
        energies = _parse_ensemble(work / "crest_conformers.xyz")
        best = work / "crest_best.xyz"
        best_xyz = best.read_text(encoding="utf-8", errors="replace") if best.exists() else None

    populations = _boltzmann_populations(energies, ROOM_TEMPERATURE_K)
    results = {
        "exit_code": code,
        "conformer_count": len(energies),
        "lowest_energy_hartree": min(energies) if energies else None,
        "energy_window_kcal_mol": round(populations[-1]["relative_energy_kcal_mol"], 3) if populations else None,
        "conformers": populations,
        "dominant_population": populations[0]["population"] if populations else None,
        "best_conformer_xyz": best_xyz,
        "temperature_k": ROOM_TEMPERATURE_K,
    }
    if not energies:
        warnings.append("CREST finished but produced no conformer ensemble; treat the run as failed.")
    elif len(energies) == 1:
        warnings.append("CREST found a single conformer; the molecule is either rigid or the search was too coarse.")
    if quick:
        warnings.append("--quick trades ensemble completeness for speed; rerun with quick=false before quantitative use.")
    if stderr.strip() and code != 0:
        warnings.append(f"CREST stderr: {stderr.strip().splitlines()[-1]}")
    return CalculationRecord(
        method=method,
        tool_version=status.version,
        parameters=params,
        results=results,
        status="ok" if code == 0 and energies else "error",
        warnings=warnings,
    )
