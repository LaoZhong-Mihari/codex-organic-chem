"""How much discriminating feedback does the route audit give?

Model accuracy cannot improve from a check that returns the same thing for a
right and a wrong answer. This benchmark measures exactly that property for the
tools, with no model in the loop: take curated correct routes, corrupt one
compound at a time the way a tired chemist (or a language model writing SMILES)
does, and ask whether the tool output changes.

Two arms are compared on identical corrupted routes:

* ``legacy``  : what the pre-existing tools emit for the same steps, namely
  ``chem_reaction_analyze`` mass/charge-balance status for ``reactants>>product``
  and ``chem_normalize_structure`` warnings for each compound.
* ``audit``   : ``route_check`` issues (blocking / review / info).

A corruption counts as detected when the arm's output for the corrupted route
differs from its output for the clean route. For ``audit`` the stricter
``blocking or review`` rate is reported separately.

This is a tool regression benchmark. It says nothing about whether a model
uses the feedback, and curated corruptions are not a sample of real model
errors.

    uv run python scripts/benchmark_route_check.py --output-dir outputs/route_check_benchmark
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Iterator

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from codex_organic_chem.route_check import route_check  # noqa: E402
from codex_organic_chem.service import chem_normalize_structure, chem_reaction_analyze  # noqa: E402
from codex_organic_chem.structure_compare import canonical, parse, stereo_elements  # noqa: E402

FIXTURES = ROOT / "tests" / "fixtures" / "routes"
MAX_SITES = 4


# --------------------------------------------------------------------------
# corruption operators: each yields corrupted SMILES for one compound


def _sanitized(rw: Chem.RWMol) -> str | None:
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def invert_centre(smiles: str) -> Iterator[str]:
    mol = parse(smiles)
    for element in stereo_elements(mol):
        if element["kind"] == "tetrahedral" and element["specified"]:
            copy_mol = Chem.RWMol(mol)
            copy_mol.GetAtomWithIdx(element["atom_index"]).InvertChirality()
            out = _sanitized(copy_mol)
            if out:
                yield out


def drop_stereo(smiles: str) -> Iterator[str]:
    mol = parse(smiles)
    for element in stereo_elements(mol):
        if element["kind"] == "tetrahedral" and element["specified"]:
            copy_mol = Chem.RWMol(mol)
            copy_mol.GetAtomWithIdx(element["atom_index"]).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            out = _sanitized(copy_mol)
            if out:
                yield out


def flip_double_bond(smiles: str) -> Iterator[str]:
    mol = parse(smiles)
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            copy_mol = Chem.RWMol(mol)
            target = copy_mol.GetBondWithIdx(bond.GetIdx())
            target.SetStereo(Chem.BondStereo.STEREOZ if stereo == Chem.BondStereo.STEREOE else Chem.BondStereo.STEREOE)
            try:
                out = Chem.MolToSmiles(copy_mol.GetMol())
            except Exception:
                continue
            if out != canonical(mol):
                yield out


def delete_methyl(smiles: str) -> Iterator[str]:
    """Drop one terminal carbon (a miscounted chain or a forgotten methyl)."""
    mol = parse(smiles)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            copy_mol = Chem.RWMol(mol)
            copy_mol.RemoveAtom(atom.GetIdx())
            out = _sanitized(copy_mol)
            if out and "." not in out:
                yield out


def insert_methylene(smiles: str) -> Iterator[str]:
    """Insert CH2 into a C-C single bond (homologation / wrong ring size)."""
    mol = parse(smiles)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() != Chem.BondType.SINGLE or a.GetSymbol() != "C" or b.GetSymbol() != "C":
            continue
        if a.GetIsAromatic() or b.GetIsAromatic():
            continue
        copy_mol = Chem.RWMol(mol)
        new = copy_mol.AddAtom(Chem.Atom(6))
        copy_mol.RemoveBond(a.GetIdx(), b.GetIdx())
        copy_mol.AddBond(a.GetIdx(), new, Chem.BondType.SINGLE)
        copy_mol.AddBond(new, b.GetIdx(), Chem.BondType.SINGLE)
        for idx in (a.GetIdx(), b.GetIdx()):
            copy_mol.GetAtomWithIdx(idx).SetChiralTag(
                mol.GetAtomWithIdx(idx).GetChiralTag()
            )
        out = _sanitized(copy_mol)
        if out:
            yield out


def swap_heteroatom(smiles: str) -> Iterator[str]:
    """O <-> N at a divalent position (an OH written as NH2, an ether as an amine)."""
    mol = parse(smiles)
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() or atom.GetFormalCharge():
            continue
        if atom.GetSymbol() == "O" and atom.GetDegree() <= 2 and all(
            b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()
        ):
            new_number = 7
        elif atom.GetSymbol() == "N" and atom.GetDegree() <= 2 and all(
            b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()
        ):
            new_number = 8
        else:
            continue
        copy_mol = Chem.RWMol(mol)
        copy_mol.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(new_number)
        out = _sanitized(copy_mol)
        if out:
            yield out


def move_substituent(smiles: str) -> Iterator[str]:
    """Shift a heteroatom substituent to the neighbouring carbon (same formula)."""
    mol = parse(smiles)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {"O", "N", "Br", "Cl"} or atom.IsInRing() or atom.GetIsAromatic():
            continue
        carbons = [n for n in atom.GetNeighbors() if n.GetSymbol() == "C" and not n.GetIsAromatic()]
        if len(carbons) != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), carbons[0].GetIdx()).GetBondType() != Chem.BondType.SINGLE:
            continue
        origin = carbons[0]
        for target in origin.GetNeighbors():
            if target.GetSymbol() != "C" or target.GetIsAromatic() or target.GetTotalNumHs() == 0:
                continue
            if target.GetIdx() == atom.GetIdx():
                continue
            copy_mol = Chem.RWMol(mol)
            copy_mol.RemoveBond(atom.GetIdx(), origin.GetIdx())
            copy_mol.AddBond(atom.GetIdx(), target.GetIdx(), Chem.BondType.SINGLE)
            for idx in (origin.GetIdx(), target.GetIdx()):
                copy_mol.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                copy_mol.GetAtomWithIdx(idx).SetNoImplicit(False)
                copy_mol.GetAtomWithIdx(idx).SetNumExplicitHs(0)
            out = _sanitized(copy_mol)
            if out and Chem.MolToSmiles(Chem.MolFromSmiles(out), isomericSmiles=False) != canonical(mol, stereo=False):
                yield out


OPERATORS: dict[str, Callable[[str], Iterator[str]]] = {
    "invert_one_stereocentre": invert_centre,
    "drop_one_stereo_label": drop_stereo,
    "flip_E_Z": flip_double_bond,
    "delete_one_carbon": delete_methyl,
    "insert_one_CH2": insert_methylene,
    "swap_O_and_N": swap_heteroatom,
    "move_substituent_to_adjacent_carbon": move_substituent,
}


# --------------------------------------------------------------------------
# arms


def audit_signal(route: dict[str, Any]) -> dict[str, set]:
    report = route_check(route)
    strict, loose = set(), set()
    for issue in report.get("issues", []):
        key = (issue.get("step"), issue["code"], issue["message"])
        loose.add(key)
        if issue["severity"] in {"blocking", "review"}:
            strict.add(key)
    blocking = {
        (i.get("step"), i["code"]) for i in report.get("issues", []) if i["severity"] == "blocking"
    }
    return {"strict": strict, "loose": loose, "blocking": blocking}


def legacy_signal(route: dict[str, Any]) -> set:
    compounds = {
        str(k): (v["smiles"] if isinstance(v, dict) else v) for k, v in route["compounds"].items()
    }
    signal: set = set()
    for key, smiles in compounds.items():
        record = chem_normalize_structure(smiles=smiles)
        for warning in record.get("warnings", []):
            signal.add(("normalize", key, warning))
    for step in route["steps"]:
        reactants = ".".join(compounds[str(r)] for r in step["reactants"])
        reaction = f"{reactants}>>{compounds[str(step['product'])]}"
        analysis = chem_reaction_analyze(reaction, mode="sanity_check").get("analysis", {})
        status = (analysis.get("mass_charge_balance") or {}).get("status")
        signal.add(("balance_status", step["step"], status))
    return signal


def corrupted_routes(route: dict[str, Any]) -> Iterator[tuple[str, str, dict[str, Any]]]:
    produced = [str(step["product"]) for step in route["steps"]]
    final = produced[-1]
    for key in produced:
        entry = route["compounds"][key]
        smiles = entry["smiles"] if isinstance(entry, dict) else entry
        clean = canonical(parse(smiles))
        for name, operator in OPERATORS.items():
            seen: set[str] = set()
            for bad in operator(smiles):
                if bad in seen or canonical(parse(bad)) == clean:
                    continue
                seen.add(bad)
                mutated = copy.deepcopy(route)
                if isinstance(mutated["compounds"][key], dict):
                    mutated["compounds"][key]["smiles"] = bad
                else:
                    mutated["compounds"][key] = bad
                yield name, ("final_product" if key == final else "intermediate"), mutated
                if len(seen) >= MAX_SITES:
                    break
    truncated = copy.deepcopy(route)
    truncated["steps"] = truncated["steps"][:-1]
    yield "omit_last_step", "final_product", truncated


def run(fixtures: list[Path]) -> dict[str, Any]:
    table: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    clean_summary = []
    for path in fixtures:
        route = json.loads(path.read_text(encoding="utf-8"))
        clean_audit = audit_signal(route)
        clean_legacy = legacy_signal(route)
        clean_summary.append(
            {
                "route": path.stem,
                "steps": len(route["steps"]),
                "audit_blocking_on_clean_route": len(clean_audit["blocking"]),
                "audit_review_prompts_on_clean_route": len(clean_audit["strict"]) - len(clean_audit["blocking"]),
                "legacy_steps_flagged_unbalanced_on_clean_route": sum(
                    1 for item in clean_legacy if item[0] == "balance_status" and item[2] == "unbalanced"
                ),
            }
        )
        for name, where, mutated in corrupted_routes(route):
            row = table[(name, where)]
            row["n"] += 1
            audit = audit_signal(mutated)
            if audit["blocking"] - clean_audit["blocking"]:
                row["audit_blocking"] += 1
            if audit["strict"] - clean_audit["strict"]:
                row["audit_blocking_or_review"] += 1
            if audit["loose"] != clean_audit["loose"]:
                row["audit_any_change"] += 1
            if legacy_signal(mutated) != clean_legacy:
                row["legacy_any_change"] += 1
    rows = []
    for (name, where), counts in sorted(table.items()):
        rows.append(
            {
                "corruption": name,
                "where": where,
                "n": counts["n"],
                "legacy_any_change": counts["legacy_any_change"],
                "audit_blocking": counts["audit_blocking"],
                "audit_blocking_or_review": counts["audit_blocking_or_review"],
                "audit_any_change": counts["audit_any_change"],
            }
        )
    total = {
        key: sum(row[key] for row in rows)
        for key in ("n", "legacy_any_change", "audit_blocking", "audit_blocking_or_review", "audit_any_change")
    }
    return {"clean_routes": clean_summary, "rows": rows, "total": total}


def to_markdown(result: dict[str, Any]) -> str:
    def pct(a: int, b: int) -> str:
        return f"{a}/{b} ({100 * a / b:.0f}%)" if b else "-"

    lines = [
        "# Route-audit discrimination benchmark",
        "",
        "Corruptions applied one compound at a time to curated correct routes. "
        "Detected = the arm's output differs from its output on the clean route.",
        "",
        "| corruption | where | n | legacy tools | audit: blocking | audit: blocking or review | audit: any change |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in result["rows"]:
        lines.append(
            f"| {row['corruption']} | {row['where']} | {row['n']} | {pct(row['legacy_any_change'], row['n'])} | "
            f"{pct(row['audit_blocking'], row['n'])} | {pct(row['audit_blocking_or_review'], row['n'])} | "
            f"{pct(row['audit_any_change'], row['n'])} |"
        )
    total = result["total"]
    lines.append(
        f"| **all** | | {total['n']} | {pct(total['legacy_any_change'], total['n'])} | "
        f"{pct(total['audit_blocking'], total['n'])} | {pct(total['audit_blocking_or_review'], total['n'])} | "
        f"{pct(total['audit_any_change'], total['n'])} |"
    )
    lines += ["", "## Clean routes (false-alarm side)", "", "| route | steps | audit blocking | audit review prompts | legacy steps flagged unbalanced |", "|---|---:|---:|---:|---:|"]
    for row in result["clean_routes"]:
        lines.append(
            f"| {row['route']} | {row['steps']} | {row['audit_blocking_on_clean_route']} | "
            f"{row['audit_review_prompts_on_clean_route']} | {row['legacy_steps_flagged_unbalanced_on_clean_route']} |"
        )
    lines += [
        "",
        "Limits: three curated routes; synthetic corruptions, not sampled model errors; a detected "
        "corruption means the tool says something different, not that a model acts on it. "
        "Review prompts on clean routes are legitimate chemistry questions but still cost attention.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--fixtures", type=Path, default=FIXTURES)
    args = parser.parse_args()
    fixtures = sorted(args.fixtures.glob("*.json"))
    result = run(fixtures)
    markdown = to_markdown(result)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "results.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        (args.output_dir / "results.md").write_text(markdown, encoding="utf-8")
    sys.stdout.write(markdown)


if __name__ == "__main__":
    main()
