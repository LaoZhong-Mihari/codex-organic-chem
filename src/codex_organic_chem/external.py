from __future__ import annotations

import os

from .figure_tools import ACS_PUBLICATION_STANDARD, FIGURE_TOOL_GUIDE, figure_tool_statuses
from .ocsr_adapters import OCSR_ADAPTERS, default_adapter_command
from .rdkit_tools import RDKIT_AVAILABLE, rdkit_version
from .utils import bundled_tool_prefix, executable_status


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
        "purpose": "Semiempirical GFN2-xTB geometries, per-atom reactivity, and thermochemistry.",
        "required_for": [
            "chem_compute task=xtb_opt",
            "chem_compute task=xtb_reactivity (Fukui indices, partial charges, IP/EA)",
            "chem_compute task=xtb_thermo (free energy, imaginary-mode count)",
        ],
        "binary": "xtb",
        "macos": ["brew tap grimme-lab/qc", "brew install xtb"],
        "fallback": [
            "micromamba create -y -p ~/.local/share/codex-organic-chem/conda-tools -c conda-forge xtb"
        ],
        "notes": "For larger systems, set OMP_NUM_THREADS and OMP_STACKSIZE as appropriate.",
    },
    "crest": {
        "purpose": "CREST conformer/rotamer ensemble searches built around xTB.",
        "required_for": ["chem_compute task=crest (ensemble + Boltzmann populations)"],
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
            "micromamba create -y -p ~/.local/share/codex-organic-chem/ocsr-tools/osra-osx64 -c edbeard -c mcs07 -c conda-forge --platform osx-64 osra=2.1.0",
            "micromamba install -y -p ~/.local/share/codex-organic-chem/ocsr-tools/osra-osx64 -c edbeard -c mcs07 -c conda-forge --platform osx-64 openbabel=2.4.1",
            "Set CODEX_CHEM_OSRA_PATH=/path/to/osra if installed outside PATH.",
        ],
        "notes": "On Apple Silicon macOS this uses the osx-64 conda package through Rosetta; the bundled path is auto-detected.",
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
    "decimer": {
        "purpose": "DECIMER Image Transformer SMILES recognition for chemical depictions.",
        "required_for": ["optional multi-tool OCSR ensemble"],
        "macos": [
            "Install DECIMER in a separate ML environment following upstream instructions.",
            "Set CODEX_CHEM_DECIMER_CMD='your-command --input {input}'",
        ],
        "notes": "Kept as an adapter command because TensorFlow/PyTorch stacks vary by machine.",
    },
    "molgrapher": {
        "purpose": "Graph-based molecular recognition with abbreviation/OCR metadata.",
        "required_for": ["optional multi-tool OCSR ensemble"],
        "macos": [
            "Install MolGrapher in a separate ML environment following upstream instructions.",
            "Set CODEX_CHEM_MOLGRAPHER_CMD='your-command --input {input}'",
        ],
        "notes": "Useful when atom/bond/bbox metadata is needed for scheme-level resolution.",
    },
    "openchemie": {
        "purpose": "Chemistry literature figure extraction and reaction diagram parsing.",
        "required_for": ["optional full-scheme and reaction OCSR"],
        "macos": [
            "Install OpenChemIE in a separate ML environment following upstream instructions.",
            "Set CODEX_CHEM_OPENCHEMIE_CMD='your-command --input {input}'",
        ],
        "notes": "Use as a configured adapter; unavailable status should not block other tools.",
    },
    "chemschematicresolver": {
        "purpose": "Chemical schematic diagram and label resolution.",
        "required_for": ["optional schematic label/structure resolution"],
        "macos": [
            "Install the bundled csr-osx64 environment and pyosra/osra_rgroup extension under ~/.local/share/codex-organic-chem/ocsr-tools.",
            "Set CODEX_CHEM_CSR_CMD='scripts/csr_adapter.py --input {input}' or CODEX_CHEM_ENABLE_CSR_DEFAULT=1 when you explicitly want CSR in the adapter pass.",
        ],
        "notes": "The official conda package is linux-64 only; this macOS setup uses an x86_64 Python 3.6 environment through Rosetta. Q8 crops produced no CSR candidates, so CSR is installed but off by default.",
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
    for spec in OCSR_ADAPTERS:
        key = f"{spec.name}_command"
        command = os.environ.get(spec.env_var) or default_adapter_command(spec)
        statuses[key] = {
            "name": key,
            "status": "available" if command else "unavailable",
            "env_var": spec.env_var,
            "command": command,
            "message": (
                f"Using configured/default {spec.name} adapter."
                if command
                else f"Set {spec.env_var} to enable a custom {spec.name} adapter."
            ),
            "purpose": spec.purpose,
        }
    for key, guide in INSTALL_GUIDE.items():
        status_key = "obabel" if key == "open-babel" else key
        if key in {"molscribe", "decimer", "molgrapher", "openchemie", "chemschematicresolver", "rxnscribe"}:
            status_key = f"{key}_command"
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
    ocsr_keys = [f"{spec.name}_command" for spec in OCSR_ADAPTERS]
    missing = [
        key
        for key in ["obabel", "xtb", "crest", "osra", *ocsr_keys]
        if statuses.get(key, {}).get("status") != "available"
    ]
    ready = [
        key
        for key in ["rdkit", "obabel", "xtb", "crest", "osra", *ocsr_keys]
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
            or any(statuses[key]["status"] == "available" for key in ocsr_keys),
        },
        "statuses": statuses,
        "install_guide": INSTALL_GUIDE,
        "publication_figure_guide": FIGURE_TOOL_GUIDE,
        "acs_publication_standard": ACS_PUBLICATION_STANDARD,
        "figure_tool_statuses": figure_tool_statuses(),
        "macos_one_shot": "scripts/install_external_tools_macos.sh",
        "notes": [
            "RDKit is required for the core assistant and is installed with uv.",
            "xTB and CREST add per-atom reactivity (Fukui indices), thermochemistry, and conformer ensembles; Open Babel adds format conversion.",
            "OSRA and the configured MolScribe/DECIMER/MolGrapher/OpenChemIE/ChemSchematicResolver/RxnScribe commands are optional OCSR adapters.",
            "Publication mechanism figures require explicit intermediates, atom-map anchored arrows, mechanism-relevant lone pairs, charges, and partial charges; use codex-chem mechanism-render before final ChemDraw/Illustrator/Inkscape polish.",
            "Unavailable tools are reported explicitly; the assistant should not invent results.",
        ],
    }


# xTB and CREST execution lives in semiempirical.py; this module reports what is
# installed and how to install it.
