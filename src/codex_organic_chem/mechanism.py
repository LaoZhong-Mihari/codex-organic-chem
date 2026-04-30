from __future__ import annotations

from html import escape

from .models import MechanismStep
from .reaction import analyze_reaction
from .rdkit_tools import reaction_to_svg


def _step(
    name: str,
    arrows: list[dict],
    bond_changes: list[str],
    charge_changes: list[str],
    rationale: str,
    confidence: float,
    warnings: list[str] | None = None,
) -> MechanismStep:
    return MechanismStep(
        elementary_step=name,
        electron_flow_arrows=arrows,
        bond_changes=bond_changes,
        charge_changes=charge_changes,
        rationale=rationale,
        confidence=confidence,
        warnings=warnings or [],
    )


def _classes_to_steps(classes: list[str], style: str) -> list[MechanismStep]:
    if "possible_SN2_substitution" in classes:
        return [
            _step(
                "concerted backside nucleophilic substitution",
                [{"from": "nucleophile lone pair", "to": "sigma-star C-leaving-group bond"}],
                ["form C-nucleophile bond", "break C-leaving-group bond"],
                ["nucleophile charge decreases or is neutralized", "leaving group departs as an anion"],
                "A strong nucleophile attacks an unhindered saturated carbon anti to the leaving group; bond formation and departure are concerted.",
                0.74,
                ["Verify substrate sterics and competing E2 before trusting SN2 assignment."],
            )
        ]
    if "possible_fischer_esterification" in classes:
        return [
            _step(
                "carbonyl activation by protonation",
                [{"from": "carbonyl oxygen lone pair", "to": "acidic proton"}],
                ["increase electrophilicity of acyl carbon"],
                ["carbonyl oxygen becomes cationic after protonation"],
                "Acid catalysis makes the carboxylic acid carbonyl more susceptible to alcohol attack.",
                0.66,
            ),
            _step(
                "alcohol addition to acyl carbon",
                [{"from": "alcohol oxygen lone pair", "to": "carbonyl carbon"}, {"from": "pi C=O", "to": "oxygen"}],
                ["form C-O bond", "convert carbonyl to tetrahedral intermediate"],
                ["oxygen charge is redistributed through proton transfers"],
                "The alcohol acts as nucleophile and forms the tetrahedral intermediate.",
                0.66,
            ),
            _step(
                "proton transfers and water departure",
                [{"from": "oxygen lone pair", "to": "C-O leaving bond"}],
                ["re-form C=O", "break C-OH2 bond"],
                ["regenerate acid catalyst"],
                "Proton transfer converts hydroxide into water; collapse of the intermediate gives ester.",
                0.62,
                ["Equilibrium and water removal dominate feasibility."],
            ),
        ]
    if "possible_suzuki_coupling" in classes:
        return [
            _step(
                "oxidative addition",
                [{"from": "Pd(0)", "to": "aryl-halide sigma-star"}],
                ["form Pd-aryl and Pd-halide bonds", "break aryl-halide bond"],
                ["Pd oxidation state increases"],
                "The aryl halide enters the palladium catalytic cycle through oxidative addition.",
                0.58,
            ),
            _step(
                "base-assisted transmetalation",
                [{"from": "aryl-boron species", "to": "palladium center"}],
                ["form Pd-aryl bond", "transfer aryl group from boron to palladium"],
                ["boron byproducts form under basic conditions"],
                "Base activates the boronic acid and enables aryl transfer to palladium.",
                0.55,
            ),
            _step(
                "reductive elimination",
                [{"from": "Pd-C bonds", "to": "new C-C bond"}],
                ["form biaryl C-C bond", "release Pd(0)"],
                ["Pd oxidation state decreases"],
                "Two organic ligands on palladium couple and regenerate the catalyst.",
                0.55,
                ["Actual ligand, base, and oxidative-addition scope require literature evidence."],
            ),
        ]
    if "possible_aldol_family" in classes:
        return [
            _step(
                "enolate or enol equivalent formation",
                [{"from": "base", "to": "alpha proton"}, {"from": "C-H bond", "to": "C-C/O pi system"}],
                ["form C=C/O enolate resonance pair"],
                ["alpha carbon or oxygen bears increased electron density"],
                "A carbonyl alpha position is activated as the nucleophilic partner.",
                0.55,
            ),
            _step(
                "carbonyl addition",
                [{"from": "enolate carbon", "to": "electrophilic carbonyl carbon"}, {"from": "pi C=O", "to": "oxygen"}],
                ["form C-C bond", "form alkoxide"],
                ["carbonyl oxygen becomes alkoxide until protonated"],
                "The enolate attacks a second carbonyl compound to give the aldol skeleton.",
                0.55,
            ),
        ]
    if "possible_electrophilic_aromatic_substitution" in classes:
        return [
            _step(
                "electrophile generation",
                [{"from": "activator", "to": "precursor leaving group"}],
                ["generate electrophilic species"],
                ["electrophile becomes strongly electron deficient"],
                "The reaction first forms a sufficiently reactive electrophile.",
                0.52,
            ),
            _step(
                "aromatic attack and sigma-complex formation",
                [{"from": "aromatic pi bond", "to": "electrophile"}],
                ["form aryl-electrophile bond", "temporarily lose aromaticity"],
                ["positive charge delocalizes in sigma complex"],
                "The arene pi system attacks according to directing effects and sterics.",
                0.52,
            ),
            _step(
                "deprotonation restores aromaticity",
                [{"from": "base", "to": "aryl proton"}, {"from": "C-H bond", "to": "aromatic pi system"}],
                ["restore aromatic pi system"],
                ["neutral arene product forms"],
                "Loss of proton restores aromatic stabilization.",
                0.52,
            ),
        ]
    return [
        _step(
            "reactive-center hypothesis",
            [{"from": "electron-rich site", "to": "electron-poor site"}],
            ["candidate bond changes require atom mapping or expert review"],
            ["formal charge changes are unknown"],
            "No high-confidence rule matched. Use atom mapping, literature precedents, or quantum calculations before drawing a detailed mechanism.",
            0.32,
            ["Unclassified mechanism draft is intentionally conservative."],
        )
    ]


def _mechanism_svg(reaction_smiles: str, steps: list[MechanismStep], step_index: int | None = None) -> str | None:
    base_svg, warnings = reaction_to_svg(reaction_smiles, width=1040, height=320)
    if not base_svg:
        return None
    selected_steps = [steps[step_index]] if step_index is not None else steps[:3]
    label = "; ".join(step.elementary_step for step in selected_steps)
    arrow_y = 52 + (24 * (step_index or 0))
    overlay = (
        f'<path d="M 130 {arrow_y} C 290 10, 470 10, 650 {arrow_y}" '
        'stroke="#d62728" stroke-width="4" fill="none" marker-end="url(#arrowhead)"/>'
        f'<text x="28" y="292" font-size="16" font-family="Arial, Helvetica, sans-serif" fill="#222">{escape(label)}</text>'
    )
    defs = (
        '<defs><marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9" refY="3.5" '
        'orient="auto"><polygon points="0 0, 10 3.5, 0 7" fill="#d62728"/></marker></defs>'
    )
    svg_start = base_svg.find("<svg")
    if svg_start >= 0:
        insert_at = base_svg.find(">", svg_start)
        if insert_at >= 0:
            base_svg = base_svg[: insert_at + 1] + defs + base_svg[insert_at + 1 :]
        base_svg = base_svg.replace("</svg>", overlay + "</svg>")
    return base_svg


def _explanation_for_step(step: MechanismStep, index: int) -> dict:
    return {
        "step_number": index,
        "elementary_step": step.elementary_step,
        "what_moves": step.electron_flow_arrows,
        "what_bonds_change": step.bond_changes,
        "what_charges_change": step.charge_changes,
        "why_it_is_plausible": step.rationale,
        "what_to_verify": step.warnings,
        "confidence": step.confidence,
    }


def _publication_package(reaction_smiles: str, steps: list[MechanismStep], confirmed: bool, classes: list[str]) -> dict:
    panels = []
    for index, step in enumerate(steps):
        panels.append(
            {
                "panel_id": f"mechanism-step-{index + 1}",
                "title": step.elementary_step,
                "vector_svg": _mechanism_svg(reaction_smiles, steps, step_index=index),
                "caption": step.rationale,
                "editable_notes": [
                    "SVG is vector output suitable for Illustrator/Inkscape cleanup.",
                    "Arrow endpoints are schematic unless atom-mapped coordinates are supplied.",
                ],
            }
        )
    blocking = []
    if not confirmed:
        blocking.append("Input structure/reaction has not been explicitly confirmed by the user.")
    blocking.append(
        "Template mechanism output is not a final publication mechanism; use chem_mechanism_render with explicit "
        "intermediates, lone pairs, charges, partial charges, and atom-map anchored arrows."
    )
    if not any(cls != "unclassified_rule_based" for cls in classes):
        blocking.append("Reaction class is unclassified; publication-level mechanism requires literature or expert mechanism evidence.")
    if not any(step.confidence >= 0.6 for step in steps):
        blocking.append("Mechanistic confidence is low; do not use as a final figure without additional evidence.")
    return {
        "publication_ready": False,
        "recommended_next_tool": "chem_mechanism_render",
        "quality_standard": "manuscript_vector_package",
        "panels": panels,
        "figure_caption": (
            "Proposed stepwise mechanism. Electron-flow arrows indicate the hypothesized movement of electron pairs; "
            "bond and charge changes are listed per step."
        ),
        "author_checklist": [
            "Confirm all substrates/products and stereochemical labels against the original source.",
            "Replace schematic arrow endpoints with atom-exact endpoints if preparing a final manuscript figure.",
            "Cite literature precedents for the mechanism class and any unusual selectivity or reagent role.",
            "Check mass/charge balance, proton transfers, counterions, solvent, and omitted catalysts.",
        ],
        "blocking_issues": blocking,
    }


def draft_mechanism(
    reaction_smiles: str,
    style: str = "stepwise",
    quality: str = "draft",
    structure_confirmed: bool = False,
) -> dict:
    analysis = analyze_reaction(reaction_smiles, mode="sanity_check")
    classes = analysis.analysis.get("reaction_classes", [])
    steps = _classes_to_steps(classes, style)
    svg = _mechanism_svg(reaction_smiles, steps)
    if svg:
        for step in steps:
            step.rendered_svg = svg
    if quality == "publication" and not structure_confirmed:
        preview, draw_warnings = reaction_to_svg(reaction_smiles)
        return {
            "status": "awaiting_user_confirmation",
            "reaction_preview": {"reaction_smiles": reaction_smiles, "svg": preview},
            "confirmation_required": True,
            "next_action": (
                "Show the rendered reaction to the user. Continue to publication-quality mechanism output only after "
                "the user confirms the reaction is correct or supplies a corrected reaction SMILES/Rxnfile."
            ),
            "warnings": [
                "Publication-quality mechanism generation is blocked until input confirmation.",
                *draw_warnings,
                *analysis.warnings,
            ],
        }
    result = {
        "status": "ok",
        "reaction": analysis.to_dict(),
        "style": style,
        "quality": quality,
        "steps": [step.to_dict() for step in steps],
        "step_explanations": [_explanation_for_step(step, index) for index, step in enumerate(steps, start=1)],
        "warnings": [
            "Mechanism steps are rule-based hypotheses and need expert/literature validation.",
            "Electron-flow arrows are publication-vector SVG assets but schematic unless atom-exact endpoints/literature support are supplied.",
            *analysis.warnings,
        ],
        "evidence_boundary": {
            "tool_facts": "RDKit parsing/drawing plus local functional-group heuristics.",
            "rule_inferences": "Mechanism templates keyed from conservative reaction-class detection.",
            "llm_assumptions": "No paid reaction database or experimental validation was used.",
        },
    }
    if quality == "publication":
        result["publication_package"] = _publication_package(
            reaction_smiles=reaction_smiles,
            steps=steps,
            confirmed=structure_confirmed,
            classes=classes,
        )
    return result
