from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
from typing import Any


JsonDict = dict[str, Any]


def _clean(value: Any) -> Any:
    if is_dataclass(value):
        return _clean(asdict(value))
    if isinstance(value, dict):
        return {k: _clean(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_clean(v) for v in value]
    return value


@dataclass
class ToolStatus:
    name: str
    status: str
    version: str | None = None
    path: str | None = None
    message: str | None = None

    def to_dict(self) -> JsonDict:
        return _clean(self)


@dataclass
class MoleculeRecord:
    source: str
    canonical_smiles: str | None = None
    isomeric_smiles: str | None = None
    inchi_key: str | None = None
    molblock: str | None = None
    svg: str | None = None
    confidence: float = 0.0
    warnings: list[str] = field(default_factory=list)
    metadata: JsonDict = field(default_factory=dict)

    def to_dict(self) -> JsonDict:
        return _clean(self)


@dataclass
class ReactionRecord:
    source: str
    reaction_smiles: str | None = None
    reactants: list[JsonDict] = field(default_factory=list)
    reagents: list[JsonDict] = field(default_factory=list)
    products: list[JsonDict] = field(default_factory=list)
    atom_mapping: JsonDict = field(default_factory=dict)
    conditions: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    confidence: float = 0.0
    warnings: list[str] = field(default_factory=list)
    analysis: JsonDict = field(default_factory=dict)

    def to_dict(self) -> JsonDict:
        return _clean(self)


@dataclass
class CalculationRecord:
    method: str
    tool_version: str | None = None
    parameters: JsonDict = field(default_factory=dict)
    input_files: list[str] = field(default_factory=list)
    output_files: list[str] = field(default_factory=list)
    results: JsonDict = field(default_factory=dict)
    status: str = "ok"
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        return _clean(self)


@dataclass
class MechanismStep:
    elementary_step: str
    electron_flow_arrows: list[JsonDict] = field(default_factory=list)
    bond_changes: list[str] = field(default_factory=list)
    charge_changes: list[str] = field(default_factory=list)
    rationale: str = ""
    rendered_svg: str | None = None
    confidence: float = 0.0
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        return _clean(self)

