from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shlex


@dataclass(frozen=True)
class OcsrAdapterSpec:
    name: str
    env_var: str
    kinds: tuple[str, ...]
    purpose: str


OCSR_ADAPTERS: tuple[OcsrAdapterSpec, ...] = (
    OcsrAdapterSpec(
        name="molscribe",
        env_var="CODEX_CHEM_MOLSCRIBE_CMD",
        kinds=("molecule",),
        purpose="Modern molecule image-to-graph recognition.",
    ),
    OcsrAdapterSpec(
        name="decimer",
        env_var="CODEX_CHEM_DECIMER_CMD",
        kinds=("molecule",),
        purpose="DECIMER Image Transformer SMILES recognition for chemical depictions.",
    ),
    OcsrAdapterSpec(
        name="molgrapher",
        env_var="CODEX_CHEM_MOLGRAPHER_CMD",
        kinds=("molecule",),
        purpose="Graph-based molecular recognition with abbreviation OCR metadata.",
    ),
    OcsrAdapterSpec(
        name="openchemie",
        env_var="CODEX_CHEM_OPENCHEMIE_CMD",
        kinds=("molecule", "reaction"),
        purpose="Chemistry literature figure extraction and reaction diagram parsing.",
    ),
    OcsrAdapterSpec(
        name="chemschematicresolver",
        env_var="CODEX_CHEM_CSR_CMD",
        kinds=("molecule",),
        purpose="Chemical schematic diagram and label resolution.",
    ),
    OcsrAdapterSpec(
        name="rxnscribe",
        env_var="CODEX_CHEM_RXNSCRIBE_CMD",
        kinds=("reaction",),
        purpose="Reaction diagram parsing.",
    ),
)


def adapter_specs_for_kind(kind: str) -> tuple[OcsrAdapterSpec, ...]:
    if kind == "auto":
        return OCSR_ADAPTERS
    return tuple(spec for spec in OCSR_ADAPTERS if kind in spec.kinds)


def default_adapter_command(spec: OcsrAdapterSpec) -> str | None:
    if spec.name != "chemschematicresolver":
        return None
    if os.environ.get("CODEX_CHEM_ENABLE_CSR_DEFAULT") != "1":
        return None
    csr_python = (
        Path.home()
        / ".local"
        / "share"
        / "codex-organic-chem"
        / "ocsr-tools"
        / "csr-osx64"
        / "bin"
        / "python"
    )
    repo_root = Path(__file__).resolve().parents[2]
    adapter = repo_root / "scripts" / "csr_adapter.py"
    if csr_python.exists() and adapter.exists():
        return f"{shlex.quote(str(csr_python))} {shlex.quote(str(adapter))} --input {{input}}"
    return None
