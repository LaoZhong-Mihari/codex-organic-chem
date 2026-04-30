from __future__ import annotations

import os
import tempfile
from pathlib import Path

from .figure_tools import ACS_PUBLICATION_STANDARD, FIGURE_TOOL_GUIDE, figure_tool_statuses
from .models import CalculationRecord
from .rdkit_tools import RDKIT_AVAILABLE, parse_molecule, rdkit_version
from .utils import bundled_tool_prefix, executable_status, run_command

if RDKIT_AVAILABLE:  # pragma: no cover - covered through higher-level tests
    from rdkit import Chem
    from rdkit.Chem import AllChem
else:  # pragma: no cover
    Chem = None
    AllChem = None


INSTALL_GUIDE = {
    "rdkit": {
        "purpose": "Core molecule parsing, validation, drawing, descriptors, conformers, and charges.",
        "required_for": ["normalize", "draw", "compute descriptors/conformers/charges", "reaction heuristics"],
        "macos": ["uv sync --extra dev --extra mcp"],
        "notes": "Installed as the Python rdkit package in this project/tool environment.",
    },
    "open-babel": {
        "purpose": "Format conversion and interoperability with legacy chemistry files.",
        "required_for": ["future format conversion workflows"],
        "binary": "obabel",
        "macos": ["brew install open-babel"],
        "notes": "Homebrew formula is open-babel; executable is obabel.",
    },
    "xtb": {
        "purpose": "Semiempirical GFN2-xTB optimization and energy estimates.",
        "required_for": ["chem_compute task=xtb_opt"],
        "binary": "xtb",
        "macos": ["brew tap grimme-lab/qc", "brew install xtb"],
        "fallback": [
            "micromamba create -y -p ~/.local/share/codex-organic-chem/conda-tools -c conda-forge xtb"
        ],
        "notes": "For larger systems, set OMP_NUM_THREADS and OMP_STACKSIZE as appropriate.",
    },
    "crest": {
        "purpose": "CREST conformer/rotamer ensemble searches built around xTB.",
        "required_for": ["chem_compute task=crest"],
        "binary": "crest",
        "macos": ["brew tap grimme-lab/qc", "brew install crest"],
        "fallback": [
            "brew install micromamba",
            "micromamba create -y -p ~/.local/share/codex-organic-chem/conda-tools -c conda-forge crest xtb",
        ],
        "notes": "On macOS, conda-forge is often more reliable than compiling from Homebrew source.",
    },
    "osra": {
        "purpose": "Traditional optical chemical structure recognition fallback.",
        "required_for": ["chem_parse_image fallback when MolScribe/RxnScribe adapters are absent"],
        "binary": "osra",
        "macos": [
            "Install Docker Desktop, then use an OSRA Docker image/wrapper, or build OSRA from source.",
            "Set CODEX_CHEM_OSRA_PATH=/path/to/osra if installed outside PATH.",
        ],
        "notes": "No stable Homebrew formula was available during setup on this machine.",
    },
    "molscribe": {
        "purpose": "Modern molecule image-to-structure recognition.",
        "required_for": ["higher-quality chem_parse_image molecule OCSR"],
        "macos": [
            "Install MolScribe in a separate ML environment following its upstream instructions.",
            "Set CODEX_CHEM_MOLSCRIBE_CMD='your-command --input {input}'",
        ],
        "notes": "Kept as an adapter command because model weights and PyTorch versions vary by machine.",
    },
    "rxnscribe": {
        "purpose": "Modern reaction scheme image recognition.",
        "required_for": ["higher-quality chem_parse_image reaction OCSR"],
        "macos": [
            "Install RxnScribe/OpenChemIE-style reaction OCR in a separate ML environment.",
            "Set CODEX_CHEM_RXNSCRIBE_CMD='your-command --input {input}'",
        ],
        "notes": "Kept as an adapter command because reaction OCR stacks are heavy and change quickly.",
    },
}


def tool_statuses() -> dict[str, dict]:
    statuses = {
        "rdkit": {
            "name": "rdkit",
            "status": "available" if RDKIT_AVAILABLE else "unavailable",
            "version": rdkit_version(),
            "message": None if RDKIT_AVAILABLE else "Install the rdkit Python package.",
        },
        "osra": executable_status("osra", ("--version",)).to_dict(),
        "obabel": executable_status("obabel", ("-V",)).to_dict(),
        "xtb": executable_status("xtb", ("--version",)).to_dict(),
        "crest": executable_status("crest", ("--version",)).to_dict(),
        "molscribe_command": {
            "name": "molscribe_command",
            "status": "available" if os.environ.get("CODEX_CHEM_MOLSCRIBE_CMD") else "unavailable",
            "message": "Set CODEX_CHEM_MOLSCRIBE_CMD to enable a custom MolScribe adapter.",
        },
        "rxnscribe_command": {
            "name": "rxnscribe_command",
            "status": "available" if os.environ.get("CODEX_CHEM_RXNSCRIBE_CMD") else "unavailable",
            "message": "Set CODEX_CHEM_RXNSCRIBE_CMD to enable a custom RxnScribe adapter.",
        },
    }
    for key, guide in INSTALL_GUIDE.items():
        status_key = "obabel" if key == "open-babel" else key
        if status_key in statuses:
            statuses[status_key]["purpose"] = guide["purpose"]
            statuses[status_key]["install"] = {
                "macos": guide.get("macos", []),
                "fallback": guide.get("fallback", []),
                "notes": guide.get("notes"),
            }
    statuses["tool_env"] = {
        "name": "codex-organic-chem conda-tools",
        "status": "available" if bundled_tool_prefix().exists() else "unavailable",
        "path": str(bundled_tool_prefix()),
        "message": "Used automatically for crest/xtb when binaries are not on PATH.",
    }
    return statuses


def doctor_report() -> dict:
    statuses = tool_statuses()
    missing = [
        key
        for key in ["obabel", "xtb", "crest", "osra", "molscribe_command", "rxnscribe_command"]
        if statuses.get(key, {}).get("status") != "available"
    ]
    ready = [
        key
        for key in ["rdkit", "obabel", "xtb", "crest", "osra", "molscribe_command", "rxnscribe_command"]
        if statuses.get(key, {}).get("status") == "available"
    ]
    return {
        "summary": {
            "ready": ready,
            "missing_or_optional": missing,
            "core_rdkit_ready": statuses["rdkit"]["status"] == "available",
            "calculation_stack_ready": statuses["xtb"]["status"] == "available"
            and statuses["crest"]["status"] == "available",
            "image_ocr_stack_ready": statuses["osra"]["status"] == "available"
            or statuses["molscribe_command"]["status"] == "available"
            or statuses["rxnscribe_command"]["status"] == "available",
        },
        "statuses": statuses,
        "install_guide": INSTALL_GUIDE,
        "publication_figure_guide": FIGURE_TOOL_GUIDE,
        "acs_publication_standard": ACS_PUBLICATION_STANDARD,
        "figure_tool_statuses": figure_tool_statuses(),
        "macos_one_shot": "scripts/install_external_tools_macos.sh",
        "notes": [
            "RDKit is required for the core assistant and is installed with uv.",
            "Open Babel, xTB, and CREST add conversion and computation capability.",
            "OSRA/MolScribe/RxnScribe are optional OCSR adapters for molecule/reaction images.",
            "Publication mechanism figures require explicit intermediates, atom-map anchored arrows, mechanism-relevant lone pairs, charges, and partial charges; use codex-chem mechanism-render before final ChemDraw/Illustrator/Inkscape polish.",
            "Unavailable tools are reported explicitly; the assistant should not invent results.",
        ],
    }


def _smiles_to_xyz(smiles: str, path: Path) -> list[str]:
    parsed = parse_molecule(smiles=smiles)
    if parsed.mol is None:
        return parsed.warnings
    mol = Chem.AddHs(parsed.mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 48879
    conf_id = AllChem.EmbedMolecule(mol, params)
    if conf_id < 0:
        return ["RDKit failed to generate a 3D conformer for xTB input."]
    try:
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
        except Exception:
            pass
    xyz = Chem.MolToXYZBlock(mol, confId=conf_id)
    path.write_text(xyz, encoding="utf-8")
    return parsed.warnings


def xtb_opt_record(smiles: str, timeout_s: int = 180) -> CalculationRecord:
    status = executable_status("xtb", ("--version",))
    params = {"smiles": smiles, "task": "opt", "method": "GFN2-xTB"}
    if status.status != "available":
        return CalculationRecord(
            method="xtb_opt",
            tool_version=None,
            parameters=params,
            status="unavailable",
            warnings=["xTB binary is not installed or not on PATH."],
        )
    if not RDKIT_AVAILABLE:
        return CalculationRecord(
            method="xtb_opt",
            tool_version=status.version,
            parameters=params,
            status="unavailable",
            warnings=["RDKit is required to build xTB input coordinates."],
        )
    with tempfile.TemporaryDirectory(prefix="codex-chem-xtb-") as tmp:
        tmp_path = Path(tmp)
        xyz_path = tmp_path / "input.xyz"
        warnings = _smiles_to_xyz(smiles, xyz_path)
        if not xyz_path.exists():
            return CalculationRecord(
                method="xtb_opt",
                tool_version=status.version,
                parameters=params,
                input_files=[str(xyz_path)],
                status="error",
                warnings=warnings,
            )
        try:
            code, stdout, stderr = run_command(
                ["xtb", str(xyz_path), "--opt", "--gfn", "2"],
                timeout_s=timeout_s,
                cwd=str(tmp_path),
            )
        except Exception as exc:  # pragma: no cover - external only
            return CalculationRecord(
                method="xtb_opt",
                tool_version=status.version,
                parameters=params,
                input_files=[str(xyz_path)],
                status="error",
                warnings=warnings + [f"xTB execution failed: {exc}"],
            )
        out_xyz = tmp_path / "xtbopt.xyz"
        energy = None
        for line in stdout.splitlines():
            if "TOTAL ENERGY" in line.upper():
                parts = line.split()
                for part in reversed(parts):
                    try:
                        energy = float(part)
                        break
                    except ValueError:
                        continue
        results = {
            "exit_code": code,
            "total_energy_hartree": energy,
            "stdout_tail": stdout.splitlines()[-20:],
            "stderr_tail": stderr.splitlines()[-20:],
        }
        if out_xyz.exists():
            results["optimized_xyz"] = out_xyz.read_text(encoding="utf-8", errors="replace")
        return CalculationRecord(
            method="xtb_opt",
            tool_version=status.version,
            parameters=params,
            input_files=[str(xyz_path)],
            output_files=[str(out_xyz)] if out_xyz.exists() else [],
            results=results,
            status="ok" if code == 0 else "error",
            warnings=warnings,
        )


def crest_record(smiles: str) -> CalculationRecord:
    status = executable_status("crest", ("--version",))
    if status.status != "available":
        return CalculationRecord(
            method="crest_conformer_search",
            tool_version=None,
            parameters={"smiles": smiles},
            status="unavailable",
            warnings=["CREST binary is not installed or not on PATH."],
        )
    return CalculationRecord(
        method="crest_conformer_search",
        tool_version=status.version,
        parameters={"smiles": smiles},
        results={
            "binary_available": True,
            "execution_policy": "manual_or_future_explicit_run",
            "suggested_manual_command": "crest input.xyz --gfn2",
        },
        status="available",
        warnings=[
            "CREST is installed, but automated CREST execution is disabled in MVP to avoid long-running jobs. "
            "Use RDKit conformers for quick estimates or run CREST manually with an explicit job policy."
        ],
    )
