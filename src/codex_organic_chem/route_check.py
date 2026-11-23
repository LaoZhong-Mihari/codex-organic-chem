"""Step-by-step audit of a proposed synthetic route.

The audit does what a careful chemist does when re-reading a route before
showing it to anyone:

* every atom of every product is accounted for (formula delta per step, read
  against the stated conditions, and exact closure when reagents/byproducts
  are declared);
* every stereocentre is followed through each step by spatial arrangement
  under an atom mapping: retained, inverted, destroyed, created, or left
  undefined;
* functional and protecting groups that are present on both sides of a step
  are listed as bystanders and compared with what the reagents attack;
* the last product is compared with the target centre by centre.

Hard facts (unparseable structures, an unclosed declared balance, a final
product that is not the target) are ``blocking``. Everything that depends on a
heuristic (MCS atom mapping, reagent keyword rules) is ``review``: a prompt to
look again, never a verdict. Passing the audit establishes bookkeeping, not
that any step works in the flask.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

from .rdkit_tools import RDKIT_AVAILABLE
from .structure_compare import (
    StructureError,
    canonical,
    compare_structures,
    delta_text,
    double_bond_relation,
    formula_counter,
    hill_formula,
    parse,
    resolve_constrained_stereo,
    stereo_elements,
    tetrahedral_relation,
    undefined_stereo_count,
)

if RDKIT_AVAILABLE:  # pragma: no cover - exercised through tests
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFMCS
else:  # pragma: no cover
    Chem = None
    Descriptors = None
    rdFMCS = None


# --------------------------------------------------------------------------
# functional groups, protecting groups and latent liabilities

FUNCTIONAL_GROUPS: dict[str, str] = {
    # oxygen
    "alcohol_primary": "[CX4;H2,H3;!$(C[OX2][#6,#14]);!$(C([OX2H])[OX2])][OX2H]",
    "alcohol_secondary": "[CX4;H1;!$(C([OX2H])[OX2,NX3])][OX2H]",
    "alcohol_tertiary": "[CX4;H0;!$(C([OX2H])[OX2,NX3])][OX2H]",
    "phenol": "c[OX2H]",
    "hemiacetal": "[CX4]([OX2H])[OX2][#6]",
    "diol_1_2": "[OX2H][CX4][CX4][OX2H]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "enone": "[CX3]=[CX3][CX3]=[OX1]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6,#1][CX3](=O)[OX2H0][#6]",
    "lactone": "[#6][CX3;R](=O)[OX2;R][#6]",
    "anhydride": "[CX3](=O)[OX2][CX3](=O)",
    "acyl_halide": "[CX3](=O)[Cl,Br]",
    "epoxide": "[CX4]1[OX2][CX4]1",
    "enol_ether": "[CX3]=[CX3][OX2][#6]",
    "acetal_ketal": "[CX4;!$(C=O)]([OX2;!$(OC=O)][#6])[OX2;!$(OC=O)][#6]",
    "aryl_methyl_ether": "c[OX2][CH3]",
    # nitrogen
    "amine_primary": "[NX3;H2;!$(NC=[O,S,N]);!$(N[SX4]);!$(N[OX2])][#6]",
    "amine_secondary": "[NX3;H1;!$(NC=[O,S,N]);!$(N[SX4]);!$(N[OX2])]([#6])[#6]",
    "amine_tertiary": "[NX3;H0;!$(NC=[O,S,N]);!$(N[SX4]);!$(N[OX2]);!$(N=*)]([#6])([#6])[#6]",
    "amide": "[#6,#1][CX3](=O)[NX3]",
    "weinreb_amide": "[CX3](=O)[NX3]([CH3])[OX2][CH3]",
    "carbamate": "[NX3][CX3](=O)[OX2H0][#6]",
    "nitrile": "[CX2]#[NX1]",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "azide": "[#6][NX2]=[NX2+]=[NX1-]",
    "imine": "[CX3]=[NX2][#6,#1]",
    "oxime": "[CX3]=[NX2][OX2]",
    "isocyanate": "[NX2]=C=O",
    "sulfonamide": "[NX3][SX4](=O)(=O)[#6]",
    # carbon pi systems, halides
    "alkene": "[CX3]=[CX3]",
    "alkyne": "[CX2]#[CX2]",
    "terminal_alkyne": "[CX2H]#[CX2]",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "aryl_halide": "c[Cl,Br,I]",
    "vinyl_halide": "[CX3]=[CX3][Cl,Br,I]",
    "sulfonate_ester": "[#6][OX2][SX4](=O)(=O)[#6]",
    "triflate": "[OX2][SX4](=O)(=O)C(F)(F)F",
    # sulfur, boron, tin, phosphorus, silicon
    "thiol": "[SX2H]",
    "sulfide": "[#6][SX2][#6]",
    "sulfoxide": "[#6][SX3](=O)[#6]",
    "sulfone": "[#6][SX4](=O)(=O)[#6]",
    "thioacetal": "[CX4]([SX2])[SX2]",
    "boronic_acid_or_ester": "[#6][BX3]([OX2])[OX2]",
    "stannane": "[#6][Sn]",
    "phosphonate": "[PX4](=O)([OX2])([OX2])[#6]",
    "silyl_enol_ether": "[CX3]=[CX3][OX2][SiX4]",
}

PROTECTING_GROUPS: dict[str, str] = {
    "silyl_ether": "[#6][OX2][SiX4]",
    "TMS_ether": "[#6][OX2][Si]([CH3])([CH3])[CH3]",
    "TES_ether": "[#6][OX2][Si]([CH2][CH3])([CH2][CH3])[CH2][CH3]",
    "TBS_ether": "[#6][OX2][Si]([CH3])([CH3])C([CH3])([CH3])[CH3]",
    "TIPS_ether": "[#6][OX2][Si](C([CH3])[CH3])(C([CH3])[CH3])C([CH3])[CH3]",
    "TBDPS_ether": "[#6][OX2][Si](c1ccccc1)(c1ccccc1)C([CH3])([CH3])[CH3]",
    "benzyl_type_ether": "[#6;!$(C=O)][OX2][CH2]c1ccccc1",
    "PMB_ether": "[#6;!$(C=O)][OX2][CH2]c1ccc(O[CH3])cc1",
    "MOM_ether": "[#6][OX2][CH2][OX2][CH3]",
    "THP_ether": "[#6][OX2][CX4;H1]1[OX2][CH2][CH2][CH2][CH2]1",
    "acetonide": "[CX4;R]([CH3])([CH3])([OX2;R])[OX2;R]",
    "acetate_ester": "[CH3][CX3](=O)[OX2][#6]",
    "benzoate_ester": "c[CX3](=O)[OX2][CX4]",
    "pivalate_ester": "[CH3]C([CH3])([CH3])[CX3](=O)[OX2][#6]",
    "methyl_ester": "[#6][CX3](=O)[OX2][CH3]",
    "tert_butyl_ester": "[#6][CX3](=O)[OX2]C([CH3])([CH3])[CH3]",
    "benzyl_ester": "[#6][CX3](=O)[OX2][CH2]c1ccccc1",
    "Boc_carbamate": "[NX3][CX3](=O)[OX2]C([CH3])([CH3])[CH3]",
    "Cbz_carbamate": "[NX3][CX3](=O)[OX2][CH2]c1ccccc1",
    "Fmoc_carbamate": "[NX3][CX3](=O)[OX2][CH2][CH]1c2ccccc2-c2ccccc12",
    "N_benzyl": "[NX3;!$(NC=O)][CH2]c1ccccc1",
    "tosyl_sulfonamide": "[NX3][SX4](=O)(=O)c1ccc([CH3])cc1",
    "dithiane": "[CX4]1[SX2][CH2][CH2][CH2][SX2]1",
}

# Groups counted when asking "only some of several equal groups reacted?"
_SELECTIVITY_GROUPS = {
    "alcohol_primary",
    "alcohol_secondary",
    "alcohol_tertiary",
    "phenol",
    "aldehyde",
    "ketone",
    "ester",
    "carboxylic_acid",
    "amide",
    "alkene",
    "alkyne",
    "amine_primary",
    "amine_secondary",
    "epoxide",
    "alkyl_halide",
    "aryl_halide",
    "azide",
    "nitro",
    "silyl_ether",
    "benzyl_type_ether",
    "carbamate",
    "acetal_ketal",
    "nitrile",
}

_COMPILED: dict[str, Any] = {}


def _pattern(name: str, smarts: str) -> Any:
    if name not in _COMPILED:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:  # pragma: no cover - guarded by tests
            raise ValueError(f"Invalid SMARTS for {name}: {smarts}")
        _COMPILED[name] = patt
    return _COMPILED[name]


def _liabilities(mol: Any) -> Counter:
    """Stereocentres and motifs that are easy to lose under ordinary conditions."""
    found: Counter = Counter()
    acceptor = Chem.MolFromSmarts("[$([CX3]=[OX1]),$([CX2]#N),$([NX3+](=O)[O-]),$([SX4](=O)=O)]")
    acceptor_atoms = {m[0] for m in mol.GetSubstructMatches(acceptor)}
    alkene_atoms = {i for m in mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=[CX3]")) for i in m}
    for element in stereo_elements(mol):
        if element["kind"] != "tetrahedral" or not element["specified"]:
            continue
        atom = mol.GetAtomWithIdx(element["atom_index"])
        if atom.GetSymbol() != "C" or atom.GetTotalNumHs() == 0:
            continue
        neighbors = {n.GetIdx() for n in atom.GetNeighbors()}
        n_acceptors = len(neighbors & acceptor_atoms)
        if n_acceptors >= 2:
            found["stereocenter_between_two_acceptors"] += 1
        elif n_acceptors == 1:
            found["epimerizable_alpha_stereocenter"] += 1
        if neighbors & alkene_atoms:
            found["allylic_stereocenter"] += 1
    beta_hydroxy = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][OX2H]")
    found["beta_hydroxy_carbonyl"] += len(mol.GetSubstructMatches(beta_hydroxy))
    beta_leaving = Chem.MolFromSmarts("[CX3](=[OX1])[CX4;!H0][CX4][$([OX2;H0]),F,Cl,Br,I]")
    found["beta_leaving_group_carbonyl"] += len(mol.GetSubstructMatches(beta_leaving))
    return +found


def functional_group_counts(mol: Any) -> Counter:
    counts: Counter = Counter()
    for table in (FUNCTIONAL_GROUPS, PROTECTING_GROUPS):
        for name, smarts in table.items():
            matches = mol.GetSubstructMatches(_pattern(name, smarts))
            if matches:
                counts[name] = len(matches)
    counts.update(_liabilities(mol))
    return counts


def functional_groups(smiles: str) -> dict[str, Any]:
    """Functional groups, protecting groups and latent liabilities of one structure."""
    mol = parse(smiles)
    counts = functional_group_counts(mol)
    return {
        "status": "ok",
        "smiles": canonical(mol),
        "functional_groups": {k: v for k, v in sorted(counts.items()) if k in FUNCTIONAL_GROUPS},
        "protecting_groups": {k: v for k, v in sorted(counts.items()) if k in PROTECTING_GROUPS},
        "liabilities": {
            k: v for k, v in sorted(counts.items()) if k not in FUNCTIONAL_GROUPS and k not in PROTECTING_GROUPS
        },
    }


# --------------------------------------------------------------------------
# reagent rules: which bystander groups a reagent class is known to attack

_ReagentRule = tuple[str, str, dict[str, str]]

REAGENT_RULES: list[_ReagentRule] = [
    (
        "LiAlH4-type hydride",
        r"\blialh4\b|\blah\b|red-al|lithium alumin",
        {
            g: "reduced or displaced by strong aluminium hydrides"
            for g in (
                "ester lactone amide carbamate nitrile aldehyde ketone enone carboxylic_acid epoxide "
                "alkyl_halide sulfonate_ester azide nitro acyl_halide anhydride imine weinreb_amide"
            ).split()
        },
    ),
    (
        "borohydride",
        r"\bnabh4\b|\blibh4\b|selectride|\bzn\(bh4\)2\b|sodium borohydride",
        {
            "aldehyde": "reduced",
            "ketone": "reduced",
            "enone": "1,2- vs 1,4-reduction must be controlled",
            "acyl_halide": "reduced",
            "imine": "reduced",
        },
    ),
    (
        "DIBAL",
        r"\bdibal",
        {
            g: "reduced by DIBAL"
            for g in "ester lactone nitrile aldehyde ketone enone amide weinreb_amide epoxide".split()
        },
    ),
    (
        "borane",
        r"\bbh3\b|9-bbn|catecholborane|\b\(sia\)2bh\b|disiamylborane",
        {
            "alkene": "hydroborated",
            "alkyne": "hydroborated",
            "carboxylic_acid": "reduced by BH3",
            "amide": "reduced by BH3",
            "aldehyde": "reduced",
            "ketone": "reduced",
        },
    ),
    (
        "catalytic hydrogenation",
        r"(?<![a-z0-9])h2(?![a-z0-9])|hydrogenat|pd/c|pd\(oh\)2|pearlman|\bpto2\b|raney",
        {
            "alkene": "hydrogenated",
            "alkyne": "hydrogenated",
            "benzyl_type_ether": "hydrogenolysed",
            "benzyl_ester": "hydrogenolysed",
            "Cbz_carbamate": "hydrogenolysed",
            "N_benzyl": "hydrogenolysed (slower)",
            "nitro": "reduced to the amine",
            "azide": "reduced to the amine",
            "aryl_halide": "hydrodehalogenated (Br, I; Cl slower)",
            "imine": "reduced",
            "nitrile": "reduced (Raney Ni, forcing Pd)",
            "epoxide": "may be hydrogenolysed at benzylic/allylic positions",
            "allylic_stereocenter": "Pd/C isomerises alkenes; an adjacent stereocentre can epimerise",
        },
    ),
    (
        "fluoride",
        r"\btbaf\b|hf[·. ]*pyr|hf-pyr|\btasf\b|\bcsf\b|\bkf\b|\bhf\b",
        {
            "silyl_ether": "cleaved by fluoride (rate TMS > TES > TBS > TIPS > TBDPS)",
            "silyl_enol_ether": "cleaved by fluoride",
            "epimerizable_alpha_stereocenter": "TBAF is basic; alpha-stereocentres can epimerise",
            "beta_leaving_group_carbonyl": "basic TBAF can trigger beta-elimination",
            "acetate_ester": "acyl migration to a freed neighbouring hydroxyl is common",
        },
    ),
    (
        "protic or Lewis acid",
        r"\btfa\b|\bhcl\b|\bh2so4\b|\btsoh\b|p-tsoh|\bcsa\b|\bhbr\b|\bbf3|\btmsotf\b|\bticl4\b|"
        r"\balcl3\b|\bsncl4\b|\bh3o\+|acoh.*h2o|amberlyst|dowex|\bhclo4\b|\btfoh\b|\bmsoh\b|\bznbr2\b|"
        # HCl generated in situ: SOCl2 or AcCl in an alcohol
        r"\b(?:socl2|accl)\b.*\b(?:meoh|etoh|methanol|ethanol)\b|\b(?:meoh|etoh|methanol|ethanol)\b.*\b(?:socl2|accl)\b",
        {
            "Boc_carbamate": "cleaved by acid",
            "tert_butyl_ester": "cleaved by acid",
            "acetal_ketal": "hydrolysed or exchanged by acid",
            "acetonide": "cleaved by aqueous acid",
            "THP_ether": "cleaved by acid",
            "MOM_ether": "cleaved by strong/Lewis acid",
            "TMS_ether": "very acid-labile",
            "TES_ether": "acid-labile",
            "TBS_ether": "cleaved by aqueous or Lewis acid (primary fastest)",
            "PMB_ether": "cleaved by strong acid",
            "silyl_enol_ether": "hydrolysed",
            "enol_ether": "hydrolysed",
            "epoxide": "opened under acid",
            "alcohol_tertiary": "can ionise/eliminate under strong acid",
            "hemiacetal": "exchanges under acid",
        },
    ),
    (
        "mild acid catalyst (PPTS)",
        r"\bppts\b|pyridinium p-toluenesulfonate",
        {
            "TMS_ether": "cleaved even by PPTS",
            "TES_ether": "slowly cleaved by PPTS in alcohols",
            "THP_ether": "cleaved by PPTS in alcohols",
            "enol_ether": "hydrolysed/exchanged",
            "silyl_enol_ether": "hydrolysed",
            "acetal_ketal": "exchanges under the same conditions",
        },
    ),
    (
        "hydroxide/alkoxide (or carbonate in alcohol/water)",
        r"\blioh\b|\bnaoh\b|\bkoh\b|\bnaome\b|\bnaoet\b|\bba\(oh\)2\b|\bme3snoh\b|\blio?oh\b|"
        r"(k2co3|cs2co3|na2co3|nahco3).*(meoh|etoh|h2o|methanol|ethanol|water)|"
        r"(meoh|etoh|h2o|methanol|ethanol|water).*(k2co3|cs2co3|na2co3)",
        {
            "ester": "saponified or transesterified",
            "acetate_ester": "cleaved (K2CO3/MeOH is the standard removal)",
            "benzoate_ester": "cleaved",
            "methyl_ester": "saponified",
            "lactone": "opened",
            "Fmoc_carbamate": "cleaved by base",
            "acyl_halide": "hydrolysed",
            "anhydride": "hydrolysed",
            "epimerizable_alpha_stereocenter": "can epimerise under basic conditions",
            "stereocenter_between_two_acceptors": "epimerises readily under base",
            "beta_hydroxy_carbonyl": "retro-aldol or dehydration under base",
            "beta_leaving_group_carbonyl": "E1cB elimination under base",
            "alkyl_halide": "substitution/elimination",
            "sulfonate_ester": "substitution/elimination",
            "epoxide": "opened by hydroxide/alkoxide",
        },
    ),
    (
        "carbonate/amidine base",
        r"\bk2co3\b|\bcs2co3\b|\bna2co3\b|\bdbu\b|\bdbn\b|\btbd\b|\bk3po4\b",
        {
            "stereocenter_between_two_acceptors": "epimerises readily even with mild base",
            "beta_hydroxy_carbonyl": "retro-aldol or dehydration under base",
            "beta_leaving_group_carbonyl": "E1cB elimination under base",
            "Fmoc_carbamate": "cleaved by DBU/secondary amines",
        },
    ),
    (
        "strong amide/organolithium base",
        r"\blda\b|\blihmds\b|\bnahmds\b|\bkhmds\b|\bn-buli\b|\bs-buli\b|\bt-buli\b|\bbuli\b|\bmeli\b|"
        r"\bphli\b|\bnah\b|\bkh\b|\blitmp\b|\bnanh2\b|\bkotbu\b|\bt-buok\b|\bkot-bu\b",
        {
            "alcohol_primary": "deprotonated first (consumes base; alkoxide may react)",
            "alcohol_secondary": "deprotonated first (consumes base; alkoxide may react)",
            "alcohol_tertiary": "deprotonated first",
            "phenol": "deprotonated first",
            "carboxylic_acid": "deprotonated first",
            "amine_primary": "N-H deprotonated by RLi/NaH",
            "amide": "N-H or alpha C-H deprotonated",
            "terminal_alkyne": "deprotonated",
            "ketone": "enolised (or attacked by RLi)",
            "aldehyde": "attacked or enolised",
            "ester": "enolised, Claisen-type self-condensation, or attacked by RLi",
            "nitrile": "alpha-deprotonated or attacked",
            "alkyl_halide": "elimination or halogen-metal exchange",
            "aryl_halide": "halogen-metal exchange with RLi (Br, I)",
            "epoxide": "opened or isomerised",
            "epimerizable_alpha_stereocenter": "deprotonation erases the stereocentre unless re-protonated selectively",
            "stereocenter_between_two_acceptors": "deprotonated immediately",
            "beta_leaving_group_carbonyl": "enolate ejects the beta leaving group",
        },
    ),
    (
        "Grignard/organometallic nucleophile",
        r"mgbr|mgcl|mgi\b|grignard|\bcuprate\b|\br2culi\b|\bme2culi\b|\balme3\b|\bznet2\b|\bzn(?:me|et)2\b",
        {
            "alcohol_primary": "protonates the organometallic",
            "alcohol_secondary": "protonates the organometallic",
            "phenol": "protonates the organometallic",
            "carboxylic_acid": "protonates the organometallic",
            "amine_primary": "protonates the organometallic",
            "terminal_alkyne": "protonates the organometallic",
            "aldehyde": "adds",
            "ketone": "adds",
            "ester": "adds twice",
            "lactone": "adds",
            "nitrile": "adds",
            "epoxide": "opens",
            "acyl_halide": "adds",
            "weinreb_amide": "adds once",
            "enone": "1,2- vs 1,4-addition",
        },
    ),
    (
        "DDQ/CAN",
        r"\bddq\b|\bcan\b|ceric ammonium",
        {
            "PMB_ether": "cleaved oxidatively",
            "benzyl_type_ether": "electron-rich benzyl ethers are slowly cleaved",
            "enol_ether": "oxidised",
            "silyl_enol_ether": "oxidised to the enone",
        },
    ),
    (
        "peracid/epoxidising oxidant",
        r"m-?cpba|\bdmdo\b|peracetic|\btfpaa\b|\boxone\b|\bh2o2\b",
        {
            "alkene": "epoxidised (most electron-rich first)",
            "sulfide": "oxidised to sulfoxide/sulfone",
            "amine_tertiary": "N-oxidised",
            "amine_secondary": "N-oxidised",
            "ketone": "Baeyer-Villiger with peracids",
            "aldehyde": "oxidised to the acid/formate",
            "enol_ether": "epoxidised (Rubottom-type)",
            "silyl_enol_ether": "Rubottom oxidation",
            "boronic_acid_or_ester": "oxidised to the alcohol",
            "thioacetal": "oxidised",
        },
    ),
    (
        "OsO4/ozone/alkene cleavage",
        r"\boso4\b|\bo3\b|ozon|\bk2oso4|ad-mix|\bkmno4\b|\bruo4\b|\brucl3\b.*naio4",
        {
            "alkene": "every accessible alkene is dihydroxylated/cleaved",
            "alkyne": "oxidised (O3, RuO4)",
            "sulfide": "oxidised",
            "enol_ether": "cleaved",
            "diol_1_2": "cleaved if NaIO4 is present",
        },
    ),
    (
        "alcohol oxidant",
        r"\bdmp\b|dess-?martin|swern|\(cocl\)2.*dmso|\bpcc\b|\bpdc\b|\btpap\b|\bibx\b|jones|\bcro3\b|"
        r"tempo|parikh|so3[·. ]*pyr|\bnaocl\b|ley[- ]griffith",
        {
            "alcohol_primary": "oxidised",
            "alcohol_secondary": "oxidised",
            "diol_1_2": "may be cleaved (DMP/Pb(OAc)4-type) or over-oxidised",
            "sulfide": "oxidised by several of these reagents",
            "thiol": "oxidised",
            "aldehyde": "over-oxidised to the acid by Jones/aqueous oxidants",
            "amine_primary": "oxidised/condenses with the new carbonyl",
            "amine_secondary": "oxidised/condenses with the new carbonyl",
            "epimerizable_alpha_stereocenter": "a new alpha-stereocentre next to an aldehyde epimerises easily",
        },
    ),
    (
        "chlorite (Pinnick)",
        r"\bnaclo2\b|pinnick",
        {
            "alkene": "HOCl by-product adds to alkenes unless scavenged (2-methyl-2-butene)",
            "sulfide": "oxidised",
            "phenol": "ring-chlorinated by HOCl",
            "amine_primary": "oxidised/chlorinated",
            "amine_secondary": "oxidised/chlorinated",
        },
    ),
    (
        "MnO2",
        r"\bmno2\b",
        {"sulfide": "can be oxidised", "amine_primary": "benzylic/allylic amines are oxidised"},
    ),
    (
        "olefin metathesis",
        r"grubbs|hoveyda|\bhg-?ii\b|\bg-?ii\b|schrock|metathesis|\brcm\b|stewart",
        {
            "alkene": "every unhindered alkene can cross-metathesise or isomerise",
            "alkyne": "enyne metathesis",
            "amine_primary": "free amines poison Ru catalysts",
            "amine_secondary": "free amines poison Ru catalysts",
            "amine_tertiary": "basic amines slow Ru catalysts",
            "sulfide": "can coordinate Ru",
            "thiol": "poisons Ru",
        },
    ),
    (
        "Pd/Ni cross-coupling",
        r"pd\(pph3\)4|pd2\(dba\)3|pd\(oac\)2|pd\(dppf\)|pdcl2|suzuki|stille|negishi|sonogashira|heck|"
        r"buchwald|\bni\(cod\)2\b|kumada|tsuji|\bxphos\b|\bsphos\b",
        {
            "aryl_halide": "any C(sp2)-X can undergo oxidative addition (I > Br ~ OTf > Cl)",
            "vinyl_halide": "any C(sp2)-X can undergo oxidative addition",
            "triflate": "oxidative addition",
            "boronic_acid_or_ester": "transmetallates",
            "stannane": "transmetallates",
            "terminal_alkyne": "Sonogashira/Glaser coupling",
            "thiol": "poisons Pd",
            "amine_primary": "N-arylation or catalyst binding",
            "amine_secondary": "N-arylation",
            "alkene": "Heck insertion/isomerisation with Pd-H",
        },
    ),
    (
        "boron tribromide/strong Lewis acid dealkylation",
        r"\bbbr3\b|\bbcl3\b|\btmsi\b|\balcl3\b.*etsh|\bbi3\b",
        {
            g: "cleaved by BBr3-type reagents"
            for g in (
                "aryl_methyl_ether benzyl_type_ether PMB_ether MOM_ether acetal_ketal acetonide silyl_ether "
                "methyl_ester tert_butyl_ester Boc_carbamate THP_ether epoxide lactone"
            ).split()
        },
    ),
    (
        "Mitsunobu",
        r"\bdiad\b|\bdead\b|mitsunobu|\bdbad\b|\badd[p]\b",
        {
            "alcohol_primary": "every unhindered free alcohol is activated",
            "alcohol_secondary": "every unhindered free alcohol is activated (with inversion)",
            "phenol": "acts as the nucleophile",
            "carboxylic_acid": "acts as the nucleophile",
            "sulfonamide": "acts as the nucleophile",
            "azide": "reduced by PPh3 (Staudinger)",
        },
    ),
    (
        "phosphine",
        r"\bpph3\b|\bpbu3\b|\bpme3\b|staudinger",
        {"azide": "Staudinger reduction", "epoxide": "deoxygenated at high temperature"},
    ),
    (
        "acylation/silylation/alkylation of heteroatoms",
        r"\bac2o\b|\baccl\b|\bbzcl\b|\bpivcl\b|\btbscl\b|\btbsotf\b|\btescl\b|\btesotf\b|\btipscl\b|"
        r"\btipsotf\b|\btbdpscl\b|\btmscl\b|\btmsotf\b|\bmscl\b|\btscl\b|"
        r"\btf2o\b|\bbnbr\b|\bpmbcl\b|\bmomcl\b|\bmei\b|\bme2so4\b|\bbncl\b|\ballyl ?br",
        {
            "alcohol_primary": "reacts (fastest of the hydroxyls)",
            "alcohol_secondary": "reacts",
            "alcohol_tertiary": "reacts with triflates/forcing conditions",
            "phenol": "reacts",
            "amine_primary": "reacts faster than hydroxyl",
            "amine_secondary": "reacts faster than hydroxyl",
            "carboxylic_acid": "forms mixed anhydride/ester",
            "thiol": "reacts",
        },
    ),
    (
        "carbamate-forming reagent",
        r"\bboc2o\b|\bcbzcl\b|\bcbz-?osu\b|\bfmoc-?cl\b|\bfmoc-?osu\b|\balloc-?cl\b|\btroc-?cl\b",
        {
            "amine_primary": "every free amine is protected",
            "amine_secondary": "every free amine is protected",
            "phenol": "O-acylated when DMAP is present",
        },
    ),
    (
        "dissolving metal / single-electron reductant",
        r"\bsmi2\b|na/nh3|li/nh3|birch|\bzn\b.*(acoh|hcl|nh4cl)|\bfe\b.*(hcl|acoh|nh4cl)|\bsncl2\b|"
        r"\bna\(hg\)|\bmg\b.*meoh|\bli\b.*naphthalen|\blidbb\b",
        {
            "ketone": "reduced/ketyl formation",
            "aldehyde": "reduced/pinacol",
            "nitro": "reduced",
            "alkyl_halide": "reduced",
            "aryl_halide": "reduced",
            "sulfone": "desulfonylated",
            "benzyl_type_ether": "cleaved by Li/NH3",
            "alkyne": "reduced to the E-alkene by Na/NH3",
            "enone": "conjugate reduction",
            "epoxide": "opened reductively",
            "azide": "reduced",
            "tosyl_sulfonamide": "cleaved by SmI2/Na-naphthalenide",
        },
    ),
    (
        "radical conditions",
        r"\bbu3snh\b|\baibn\b|\bttmss\b|\(tms\)3sih|\bet3b\b|barton",
        {
            "alkyl_halide": "reduced/cyclises",
            "aryl_halide": "reduced (Br, I)",
            "vinyl_halide": "reduced",
            "alkene": "radical acceptor (5-exo/6-endo cyclisation)",
            "alkyne": "radical acceptor",
            "thiol": "H-atom donor",
            "nitro": "denitrated (tertiary/secondary)",
        },
    ),
    (
        "electrophilic halogen",
        r"\bnbs\b|\bncs\b|\bnis\b|\bbr2\b|\bi2\b|\bicl\b|\bselectfluor\b|\bcl2\b",
        {
            "alkene": "halofunctionalised",
            "alkyne": "halogenated",
            "phenol": "ring-halogenated",
            "enol_ether": "halogenated",
            "silyl_enol_ether": "alpha-halogenated",
            "amine_primary": "N-halogenated/oxidised",
            "amine_secondary": "N-halogenated/oxidised",
            "sulfide": "oxidised",
            "thiol": "oxidised to disulfide",
            "ketone": "alpha-halogenated under acid/base",
        },
    ),
    (
        "nucleophilic azide/cyanide",
        r"\bnan3\b|\btmsn3\b|\blin3\b|\bnacn\b|\bkcn\b|\btmscn\b|\bbu4ncn\b|\bet2alcn\b",
        {
            "alkyl_halide": "displaced",
            "sulfonate_ester": "displaced",
            "epoxide": "opened",
            "enone": "conjugate addition",
            "aldehyde": "cyanohydrin formation (cyanide)",
            "ketone": "cyanohydrin formation (cyanide)",
            "acyl_halide": "forms the acyl azide/cyanide",
        },
    ),
    (
        "thiolate",
        r"\bnaset\b|\bnasme\b|\bphsna\b|\bnasph\b|\betsna\b|\blisph\b|\bphsh\b.*(k2co3|et3n|nah)",
        {
            "alkyl_halide": "displaced",
            "sulfonate_ester": "displaced",
            "epoxide": "opened",
            "enone": "conjugate addition",
            "aryl_methyl_ether": "demethylated at high temperature",
            "methyl_ester": "demethylated (SN2 at methyl)",
            "aryl_halide": "SNAr on activated rings",
        },
    ),
    (
        "phosphorus ylide/olefination",
        r"wittig|ph3p=|\bhwe\b|horner|phosphonate.*(nah|lihmds|khmds|dbu)|tebbe|petasis|julia|peterson",
        {
            "aldehyde": "olefinated",
            "ketone": "olefinated (slower); basic ylides also enolise",
            "epimerizable_alpha_stereocenter": "basic ylides epimerise alpha-stereocentres",
            "beta_hydroxy_carbonyl": "base-mediated retro-aldol",
            "ester": "attacked by Tebbe/Petasis reagents",
            "lactone": "attacked by Tebbe/Petasis reagents",
        },
    ),
    (
        "diazomethane",
        r"\bch2n2\b|\btmschn2\b|diazomethane",
        {
            "carboxylic_acid": "methylated",
            "phenol": "methylated",
            "ketone": "homologated (slow)",
            "aldehyde": "homologated/epoxidised",
            "acyl_halide": "forms the diazoketone",
        },
    ),
    (
        "amide coupling/activating agent",
        r"\bhatu\b|\bhbtu\b|\bedc\b|\bedci\b|\bdcc\b|\bdic\b|\bpybop\b|\bt3p\b|\bcdi\b|\bbop\b|yamaguchi|"
        # SOCl2 in an alcohol is an HCl-mediated esterification, not an acid
        # chloride preparation: it has its own rule below
        r"^(?!.*\b(?:meoh|etoh|methanol|ethanol|i-?proh)\b).*\bsocl2\b|oxalyl chloride|\(cocl\)2(?!.*dmso)|shiina|\bmnba\b",
        {
            "carboxylic_acid": "every free acid is activated",
            "amine_primary": "acylated",
            "amine_secondary": "acylated",
            "alcohol_primary": "acylated (esterification conditions)",
            "alcohol_secondary": "acylated with DMAP/Yamaguchi",
            "phenol": "acylated",
            "epimerizable_alpha_stereocenter": "activated acids epimerise at the alpha-carbon (oxazolone/ketene)",
        },
    ),
    (
        "hydride transfer to iminium (reductive amination)",
        r"\bnabh3cn\b|\bnabh\(oac\)3\b|\bstab\b|picoline.?borane",
        {"aldehyde": "reductively aminated/reduced", "ketone": "reductively aminated", "imine": "reduced"},
    ),
    (
        "thermal conditions",
        r"\b(1[2-9]\d|[2-9]\d\d) ?°?c\b|reflux.*(xylene|toluene|dmf|dmso|dce|mesitylene)|sealed tube|microwave",
        {
            "Boc_carbamate": "thermolysed above ~150 C",
            "beta_hydroxy_carbonyl": "retro-aldol/dehydration on heating",
            "azide": "decomposes/cycloadds on heating",
            "epimerizable_alpha_stereocenter": "prolonged heating promotes epimerisation",
        },
    ),
]

_INVERTING_CONDITIONS = re.compile(
    r"\bdiad\b|\bdead\b|\bdbad\b|\baddp\b|mitsunobu|\bnan3\b|\blin3\b|\bdppa\b|\btmsn3\b|\bnacn\b|\bkcn\b|"
    r"\bcsoac\b|\bkoac\b|\bnaoac\b|\bbu4noac\b|\bko2\b|\bnai\b|\blibr\b|\bbu4ni\b|\btbai\b|\bnasph\b|"
    r"\bksac\b|\bkscn\b|phthalimide|\bsn2\b|\bnano2\b|appel|\bcbr4\b|\bpbr3\b|\bcsf\b.*\bdmf\b",
    re.IGNORECASE,
)

_RULES_COMPILED = [(name, re.compile(pattern, re.IGNORECASE), risks) for name, pattern, risks in REAGENT_RULES]


# --------------------------------------------------------------------------
# formula-delta dictionary (product minus the numbered reactants)

_NamedDelta = tuple[str, dict[str, int], str]


def _d(text: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for token in text.split():
        match = re.fullmatch(r"([+-]?\d*)([A-Z][a-z]?)", token)
        assert match, token
        raw = match.group(1)
        out[match.group(2)] = int(raw) if raw not in {"", "+", "-"} else (-1 if raw == "-" else 1)
    return out


# (label, delta, regex of reagent words that would justify it)
NAMED_DELTAS: list[_NamedDelta] = [
    ("no atoms gained or lost (isomerisation, rearrangement or addition)", {}, r""),
    ("+H2: reduction/hydrogenation", _d("2H"), r"h2\b|nabh4|lialh4|\blah\b|dibal|selectride|bh3|pd/c|pto2|raney|"
     r"hydrogen|smi2|nabh3cn|nabh\(oac\)3|red-al|luche|cbs|noyori|et3sih|diimide|lindlar|crabtree|nh3|zn|"
     r"stryker|l-selectride|k-selectride|ipr2?oh|hantzsch|n2h4|tsnhnh2|nbsh"),
    ("-H2: oxidation/dehydrogenation", _d("-2H"), r"dmp|dess|swern|pcc|pdc|tpap|tempo|mno2|ibx|jones|cro3|ddq|"
     r"so3|oxid|\(cocl\)2|naocl|o2\b|pd\(oac\)2|selenium|seo2|phse|saegusa|chloranil|baib|phi\(oac\)2|cu"),
    ("+O: epoxidation, hydroxylation or heteroatom oxidation", _d("O"), r"cpba|dmdo|h2o2|tbhp|t-buooh|vo\(acac\)2|"
     r"sharpless|oxone|seo2|davis|o2\b|naio4|mooph|peracid|shi |jacobsen|ti\(oipr\)4|baeyer|naclo2|pinnick|"
     r"oxaziridine|tfdo|rubottom"),
    ("-O: deoxygenation", _d("-O"), r"pph3|p\(oet\)3|zn|ti|sm|wcl6|barton|bu3snh|et3sih|deoxygen"),
    ("-H2O: dehydration, condensation or cyclisation", _d("-2H -O"), r"."),
    ("+H2O: hydration or hydrolysis", _d("2H O"), r"h2o|h3o|aq|hydro|lioh|naoh|koh|hcl|h2so4|oxymerc|bh3|9-bbn|"
     r"hg\(|wet|acoh|tfa|tsoh|ppts"),
    ("+H2O2: dihydroxylation", _d("2H 2O"), r"oso4|ad-mix|kmno4|k2oso4|nmo|upjohn|woodward|prevost|ruo4"),
    ("-CO2: decarboxylation", _d("-C -2O"), r"heat|°c|reflux|krapcho|licl|nacl|barton|decarbox|cu|ag|hcl|pd"),
    ("-N2: nitrogen extrusion", _d("-2N"), r"heat|rh2|cu|hv|light|wolff|curtius|°c|pd|staudinger|pph3|schmidt"),
    ("-N2 +H2: azide to amine", _d("-2N 2H"), r"pph3|h2\b|pd/c|lialh4|staudinger|zn|pme3|dtt|sncl2|lindlar"),
    ("-2O +2H: nitro to amine", _d("-2O 2H"), r"h2\b|fe\b|zn\b|sncl2|pd/c|raney|na2s2o4|in\b"),
    ("-HCl", _d("-H -Cl"), r"."),
    ("-HBr", _d("-H -Br"), r"."),
    ("-HI", _d("-H -I"), r"."),
    ("-HF", _d("-H -F"), r"."),
    ("H to Cl: chlorination", _d("-H Cl"), r"ncs|cl2|so2cl2|socl2|pcl5|tcca|c2cl6|chlor"),
    ("H to Br: bromination", _d("-H Br"), r"nbs|br2|cbr4|pyh.*br3|dbdmh|brom|cubr2"),
    ("H to I: iodination", _d("-H I"), r"nis|i2\b|icl|iod"),
    ("H to F: fluorination", _d("-H F"), r"selectfluor|nfsi|dast|deoxo|fluor"),
    ("+Br2: dibromination", _d("2Br"), r"br2|pyh.*br3|brom"),
    ("+Cl2: dichlorination", _d("2Cl"), r"cl2|so2cl2|chlor"),
    ("+HBr", _d("H Br"), r"hbr"),
    ("+HCl", _d("H Cl"), r"hcl"),
    ("OH to Cl", _d("-O -H Cl"), r"socl2|pcl|ccl4|ncs|mscl|\(cocl\)2|appel|pph3|tscl|ghosez"),
    ("OH to Br", _d("-O -H Br"), r"pbr3|cbr4|appel|nbs|hbr|pph3|socl?br2|libr"),
    ("OH to I", _d("-O -H I"), r"i2\b|appel|pph3|nai|hi\b|imidazole"),
    ("OH to N3", _d("-O -H 3N"), r"dppa|nan3|hn3|mitsunobu|diad|dead|tmsn3"),
    ("+CH2: methylation, methyl esterification or one-carbon homologation", _d("C 2H"), r"mei|me2so4|ch2n2|tmschn2|"
     r"meotf|ch2i2|simmons|zn|me3o|methyl|dmc|hcho|paraformaldehyde|corey|ch2br2|et2zn|memgbr|meli|arndt|"
     r"(?:socl2|h2so4|hcl|accl|tmscl|tsoh|amberlyst|dowex).*(?:\bmeoh\b|methanol)|(?:\bmeoh\b|methanol).*(?:socl2|h2so4|hcl|accl|tmscl|tsoh|amberlyst|dowex)"),
    ("-CH2: demethylation or methyl-ester cleavage", _d("-C -2H"), r"bbr3|lioh|naoh|koh|lii|naset|tmsi|bcl3|"
     r"me3snoh|hbr|alcl3|esterase|lipase|ba\(oh\)2|k2co3|nacn|licl"),
    ("-C2H4: ethyl-ester cleavage or metathesis (loss of ethylene)", _d("-2C -4H"), r"lioh|naoh|koh|grubbs|"
     r"hoveyda|metathesis|rcm|bbr3|schrock|stewart|k2co3|hcl"),
    ("+C2H4: ethylation or ethyl esterification", _d("2C 4H"), r"eti|etbr|et3o|etotf|et2so4|ethyl|etmgbr|"
     r"(?:socl2|h2so4|hcl|accl|tmscl|tsoh|amberlyst|dowex).*(?:\betoh\b|ethanol)|(?:\betoh\b|ethanol).*(?:socl2|h2so4|hcl|accl|tmscl|tsoh|amberlyst|dowex)"),
    ("+C3H4: allylation or acetonide formation", _d("3C 4H"), r"allyl|dimethoxypropane|acetone|2-methoxypropene|"
     r"dmp\b|me2c\(ome\)2"),
    ("+C2H2O: acetylation", _d("2C 2H O"), r"ac2o|accl|acetyl|vinyl acetate|ketene|isopropenyl acetate"),
    ("+C7H4O: benzoylation", _d("7C 4H O"), r"bzcl|bz2o|benzoyl|bzcn"),
    ("+C5H8O: pivaloylation or THP protection", _d("5C 8H O"), r"pivcl|piv2o|dhp|dihydropyran|pivaloyl"),
    ("+C5H8O2: Boc protection", _d("5C 8H 2O"), r"boc"),
    ("+C8H6O2: Cbz protection", _d("8C 6H 2O"), r"cbz|benzyl chloroformate|bnococl"),
    ("+C15H10O2: Fmoc protection", _d("15C 10H 2O"), r"fmoc"),
    ("+C7H6: benzylation", _d("7C 6H"), r"bnbr|bncl|benzyl|bnoc\(=?nh\)ccl3|trichloroacetimidate"),
    ("+C8H8O: PMB protection", _d("8C 8H O"), r"pmb|methoxybenzyl"),
    ("+C2H4O: MOM protection or ethylene ketal", _d("2C 4H O"), r"mom|ethylene glycol|hoch2ch2oh|ch2\(ome\)2|"
     r"\(ch2oh\)2|tmsoch2ch2otms"),
    ("+C19H14: tritylation", _d("19C 14H"), r"trcl|trityl"),
    ("+C3H8Si: TMS", _d("3C 8H Si"), r"tms"),
    ("+C6H14Si: TBS or TES silylation", _d("6C 14H Si"), r"tbs|tes|tbdms|et3si"),
    ("+C9H20Si: TIPS silylation", _d("9C 20H Si"), r"tips"),
    ("+C16H18Si: TBDPS silylation", _d("16C 18H Si"), r"tbdps|tbdpscl"),
    ("+CH2O2S: mesylation", _d("C 2H 2O S"), r"mscl|ms2o|mesyl"),
    ("+C7H6O2S: tosylation", _d("7C 6H 2O S"), r"tscl|ts2o|tosyl"),
    ("+CF3O2S -H: triflation", _d("C 3F 2O S -H"), r"tf2o|phntf2|comins|triflat|tfcl"),
    ("+NO2 -H: nitration", _d("N 2O -H"), r"hno3|no2bf4|nitrat|nano2|agno2"),
    ("+CH2 -O: carbonyl methylenation", _d("C 2H -O"), r"wittig|ph3p|tebbe|petasis|nysted|peterson|lombardo|"
     r"takai|julia|cp2ti"),
    ("+2H -O: C=O to CH2 (or C-OH to C-H)", _d("2H -O"), r"wolff|n2h4|clemmensen|zn|et3sih|barton|bu3snh|"
     r"raney|tshnnh2|nabh3cn|h2\b|pd/c|lialh4|bh3|alh3|mozingo"),
    ("-CH2 +O: alkene C=CH2 cleaved to C=O", _d("-C -2H O"), r"o3|ozon|oso4|naio4|ruo4|lemieux|kmno4"),
    ("-CH4O: loss of methanol (lactonisation, transesterification, elimination)", _d("-C -4H -O"), r"."),
    ("-C2H6O: loss of ethanol (lactonisation, Claisen/Dieckmann, transesterification)", _d("-2C -6H -O"), r"."),
    ("-C4H8: tert-butyl ester/ether cleavage", _d("-4C -8H"), r"tfa|hcl|h3po4|tmsotf|zn|acid|hco2h|sio2|tfoh"),
    ("-CH2O: methyl ester to aldehyde", _d("-C -2H -O"), r"dibal|red-al|lialh\(ot"),
    ("-CO: methyl ester to primary alcohol", _d("-C -O"), r"lialh4|\blah\b|dibal|libh4|red-al|nabh4|superhydride"),
    ("-C2H4O: ethyl ester to aldehyde", _d("-2C -4H -O"), r"dibal|red-al|lialh\(ot"),
    ("-C2H2O: ethyl ester to primary alcohol", _d("-2C -2H -O"), r"lialh4|\blah\b|dibal|libh4|red-al|nabh4"),
    ("+2H2: double reduction (alkyne to alkane, nitrile to amine, ester-level reductions)", _d("4H"),
     r"h2\b|pd/c|pto2|raney|lialh4|\blah\b|bh3|dibal|hydrogen"),
    ("+O -2H: CH2OH to CO2H, or CH2 to C=O", _d("O -2H"), r"jones|cro3|tempo|naclo2|pinnick|kmno4|ruo4|rucl3|"
     r"seo2|pdc|pcc|baib|oxid|tbhp|o2\b|naocl|h5io6"),
    ("+2O: ozonolysis keeping both carbonyls, or dioxygenation", _d("2O"), r"o3|ozon|oso4|naio4|ruo4|o2\b|kmno4"),
    ("-N +O +H: nitrile to aldehyde", _d("-N O H"), r"dibal|raney|h3o"),
    ("OMs to N3", _d("-C -3H -3O -S 3N"), r"n3|azide|dppa"),
    ("OTs to N3", _d("-7C -7H -3O -S 3N"), r"n3|azide"),
    ("OMs to I", _d("-C -3H -3O -S I"), r"nai|lii|ki\b|bu4ni|tbai"),
    ("OTs to I", _d("-7C -7H -3O -S I"), r"nai|lii|ki\b|bu4ni|tbai"),
    ("OMs to Br", _d("-C -3H -3O -S Br"), r"libr|nabr|bu4nbr|mgbr2"),
    ("OTs to Br", _d("-7C -7H -3O -S Br"), r"libr|nabr|bu4nbr|mgbr2"),
    ("OMs to CN", _d("-3H -3O -S N"), r"cn\b|cyanide"),
    ("OTs to CN", _d("-6C -7H -3O -S N"), r"cn\b|cyanide"),
    ("Br to N3", _d("-Br 3N"), r"n3|azide"),
    ("Cl to N3", _d("-Cl 3N"), r"n3|azide"),
    ("I to N3", _d("-I 3N"), r"n3|azide"),
    ("Br to CN", _d("-Br C N"), r"cn\b|cyanide"),
    ("Cl to CN", _d("-Cl C N"), r"cn\b|cyanide"),
    ("I to CN", _d("-I C N"), r"cn\b|cyanide"),
    ("Cl to I (Finkelstein)", _d("-Cl I"), r"nai|ki\b|lii"),
    ("Br to I (Finkelstein)", _d("-Br I"), r"nai|ki\b|lii"),
    ("+CO: carbonylation", _d("C O"), r"co\b|carbonyl|mo\(co\)6|pauson|co2\(co\)8|formyl|dmf"),
    ("+CO2: carboxylation", _d("C 2O"), r"co2|dry ice|carboxyl"),
    ("-B(OH)2 -Cl: Suzuki coupling", _d("-B -2O -2H -Cl"), r"pd|ni|suzuki"),
    ("-B(OH)2 -Br: Suzuki coupling", _d("-B -2O -2H -Br"), r"pd|ni|suzuki"),
    ("-B(OH)2 -I: Suzuki coupling", _d("-B -2O -2H -I"), r"pd|ni|suzuki"),
    ("-Bpin -Br: Suzuki coupling", _d("-B -2O -6C -12H -Br"), r"pd|ni|suzuki"),
    ("-Bpin -I: Suzuki coupling", _d("-B -2O -6C -12H -I"), r"pd|ni|suzuki"),
    ("-Bpin -Cl: Suzuki coupling", _d("-B -2O -6C -12H -Cl"), r"pd|ni|suzuki"),
    ("-MgBr +H: Grignard addition after work-up", _d("-Mg -Br H"), r"."),
    ("-MgCl +H: Grignard addition after work-up", _d("-Mg -Cl H"), r"."),
    ("-Li +H: organolithium addition after work-up", _d("-Li H"), r"."),
    ("-Ph3PO: Wittig olefination with a numbered ylide/phosphonium", _d("-18C -15H -P -O"), r"."),
    ("-(EtO)2PO2H: Horner-Wadsworth-Emmons", _d("-4C -11H -4O -P"), r"."),
    ("-(MeO)2PO2H: Horner-Wadsworth-Emmons", _d("-2C -7H -4O -P"), r"."),
    ("-SnBu3 -Br: Stille coupling", _d("-Sn -12C -27H -Br"), r"pd|stille|cu"),
    ("-SnBu3 -I: Stille coupling", _d("-Sn -12C -27H -I"), r"pd|stille|cu"),
    ("-SnMe3 -I: Stille coupling", _d("-Sn -3C -9H -I"), r"pd|stille|cu"),
]

_PROTECTION_WORDS = (
    "silylation",
    "protection",
    "acetylation",
    "benzoylation",
    "pivaloylation",
    "benzylation",
    "tritylation",
    "TMS",
    "mesylation",
    "tosylation",
    "triflation",
    "acetonide",
)


def _counter_key(delta: dict[str, int]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted((k, v) for k, v in delta.items() if v))


_DELTA_INDEX: dict[tuple, list[tuple[str, str]]] = {}
for _label, _delta, _words in NAMED_DELTAS:
    _DELTA_INDEX.setdefault(_counter_key(_delta), []).append((_label, _words))
    if _delta and any(word in _label for word in _PROTECTION_WORDS):
        reverse = {k: -v for k, v in _delta.items()}
        _DELTA_INDEX.setdefault(_counter_key(reverse), []).append(
            (f"reverse of [{_label}] (deprotection/removal)", r".")
        )


def interpret_delta(delta: dict[str, int], conditions: str = "") -> dict[str, Any]:
    """Name a formula delta and say whether the stated conditions account for it.

    A delta is only called explained when a reagent named in the conditions
    can cause it. Two simultaneous changes are accepted only if the conditions
    support both; otherwise a miscounted atom would always find some pair of
    textbook transformations that happens to add up.
    """

    def supported(words: str) -> bool:
        return not words or bool(re.search(words, conditions, re.IGNORECASE))

    key = _counter_key(delta)
    hits = list(_DELTA_INDEX.get(key, []))
    combos: list[str] = []
    if not hits and key:
        signed: list[tuple[str, dict[str, int], str]] = []
        for label, d, words in NAMED_DELTAS:
            if not d:
                continue
            signed.append((label, d, words))
            if any(word in label for word in _PROTECTION_WORDS):
                signed.append((f"reverse of [{label}]", {k: -v for k, v in d.items()}, r"."))
        target = Counter(dict(key))
        for i, (label_a, da, words_a) in enumerate(signed):
            if not supported(words_a):
                continue
            rest = Counter(target)
            rest.subtract(da)
            rest_key = _counter_key(dict(rest))
            for label_b, db, words_b in signed[i:]:
                if _counter_key(db) == rest_key and supported(words_b):
                    combos.append(f"{label_a}  AND  {label_b}")
            if len(combos) >= 3:
                break
    matched = [label for label, words in hits if supported(words)]
    if hits and matched:
        status = "consistent_with_conditions"
    elif hits:
        status = "named_but_conditions_do_not_mention_a_matching_reagent"
    elif combos:
        status = "explained_only_as_a_combination"
    else:
        status = "unexplained"
    return {
        "status": status,
        "matches": matched or [label for label, _ in hits],
        "possible_combinations": combos[:3],
    }


# --------------------------------------------------------------------------
# atom mapping between a reactant and the product


def _labelled(mol: Any, blocked: set[int] | None = None) -> Any:
    """Copy whose isotopes encode the element, so MCS can exclude blocked atoms."""
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetIsotope(900 + atom.GetIdx() if blocked and atom.GetIdx() in blocked else atom.GetAtomicNum())
        atom.SetAtomMapNum(0)
    return copy


def _map_number_mapping(reactant: Any, product: Any) -> dict[int, int] | None:
    r_maps = {a.GetAtomMapNum(): a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
    p_maps = {a.GetAtomMapNum(): a.GetIdx() for a in product.GetAtoms() if a.GetAtomMapNum()}
    shared = set(r_maps) & set(p_maps)
    heavy = sum(1 for a in reactant.GetAtoms() if a.GetAtomicNum() > 1)
    if not shared or len(shared) < 0.5 * heavy:
        return None
    return {r_maps[n]: p_maps[n] for n in shared}


def _with_stereo_hs(mol: Any) -> Any:
    centres = [e["atom_index"] for e in stereo_elements(mol) if e["kind"] == "tetrahedral"]
    centres = [i for i in centres if mol.GetAtomWithIdx(i).GetTotalNumHs() > 0]
    if not centres:
        return Chem.Mol(mol)
    return Chem.AddHs(mol, onlyOnAtoms=centres)


def _stereo_events(reactant: Any, product: Any, mapping: dict[int, int]) -> dict[str, Any]:
    """Follow every stereo element of the reactant into the product."""
    reactant_h, product_h = _with_stereo_hs(reactant), _with_stereo_hs(product)
    events: dict[str, list[dict[str, Any]]] = {
        "retained": [],
        "inverted": [],
        "epimerized_without_bond_change": [],
        "destroyed": [],
        "lost_definition": [],
        "ambiguous": [],
    }
    product_elements = {
        e["atom_index"]: e for e in stereo_elements(product) if e["kind"] == "tetrahedral"
    }
    for element in stereo_elements(reactant):
        if element["kind"] == "tetrahedral":
            if not element["specified"]:
                continue
            r_idx = element["atom_index"]
            entry = {
                "reactant_atom": r_idx,
                "element": element["element"],
                "substituents": element["substituents"],
                "cip_before": element["cip"],
            }
            p_idx = mapping.get(r_idx)
            if p_idx is None:
                events["ambiguous"].append({**entry, "reason": "atom not traced into the product"})
                continue
            entry["product_atom"] = p_idx
            target = product_elements.get(p_idx)
            if target is None:
                events["destroyed"].append(entry)
                continue
            entry["cip_after"] = target["cip"]
            if not target["specified"]:
                events["lost_definition"].append(entry)
                continue
            relation, detail = tetrahedral_relation(
                reactant_h.GetAtomWithIdx(r_idx),
                product_h.GetAtomWithIdx(p_idx),
                mapping,
                pair_unmapped=True,
            )
            entry.update(detail)
            if relation in {"same", "inverted"} and "replaced" not in detail:
                exchanged = _same_element_exchange(reactant, product, r_idx, p_idx, mapping)
                if exchanged:
                    entry["same_element_exchange_possible"] = exchanged
            if relation == "same":
                events["retained"].append(entry)
            elif relation == "inverted":
                # Whether a same-element exchange really happened depends on the
                # conditions; route_check re-buckets these with that knowledge.
                changed = "replaced" in entry or "same_element_exchange_possible" in entry
                events["inverted" if changed else "epimerized_without_bond_change"].append(entry)
            else:
                events["ambiguous"].append({**entry, "reason": detail.get("reason", relation)})
        else:
            if not element["specified"]:
                continue
            a1, a2 = element["atom_indices"]
            bond = reactant.GetBondBetweenAtoms(a1, a2)
            relation = double_bond_relation(bond, product, mapping)
            entry = {"reactant_bond": [a1, a2], "geometry_before": element["cip"]}
            if relation == "same":
                events["retained"].append(entry)
            elif relation == "inverted":
                events["epimerized_without_bond_change"].append({**entry, "note": "E/Z isomerised"})
            elif relation == "absent":
                events["destroyed"].append(entry)
            elif relation == "unspecified_b":
                events["lost_definition"].append(entry)
            else:
                events["ambiguous"].append({**entry, "reason": relation})
    return events


def _same_element_exchange(reactant: Any, product: Any, r_idx: int, p_idx: int, mapping: dict[int, int]) -> str | None:
    """Could a heteroatom at the centre have been swapped for another atom of the same element?

    In a Mitsunobu esterification or an SN2 by an O/N nucleophile on an O/N
    leaving group the element at the centre does not change, so an element-based
    mapping pairs the departing atom with the arriving one. The tell-tale is
    that this heteroatom carries a heavy substituent in the product that it did
    not have in the reactant (the rest of the nucleophile). O-acylation,
    sulfonylation and silylation look identical on the graph, which is why the
    caller decides with the reaction conditions.
    """
    for neighbor in reactant.GetAtomWithIdx(r_idx).GetNeighbors():
        if neighbor.GetAtomicNum() in (1, 6):
            continue
        image = mapping.get(neighbor.GetIdx())
        if image is None or product.GetBondBetweenAtoms(p_idx, image) is None:
            continue
        partner = product.GetAtomWithIdx(image)
        inherited = {mapping.get(n.GetIdx()) for n in neighbor.GetNeighbors()}
        new_heavy = [
            n for n in partner.GetNeighbors() if n.GetIdx() != p_idx and n.GetIdx() not in inherited
        ]
        if new_heavy:
            symbol = neighbor.GetSymbol()
            return f"{symbol}->{symbol} (same element; that {symbol} carries a new substituent in the product)"
    return None


def _event_cost(events: dict[str, Any]) -> tuple[int, int, int]:
    return (
        len(events["epimerized_without_bond_change"]) + len(events["ambiguous"]),
        len(events["inverted"]),
        len(events["lost_definition"]),
    )


def map_reactant_to_product(
    reactant: Any,
    product: Any,
    *,
    blocked_product_atoms: set[int] | None = None,
    timeout_s: int = 6,
) -> dict[str, Any]:
    """Atom mapping reactant -> product from atom maps, else maximum common substructure.

    Symmetric skeletons admit several equally large mappings; the one that
    implies the fewest unexplained stereo changes is chosen, so an inversion is
    reported only when no equivalent mapping avoids it.
    """
    direct = _map_number_mapping(reactant, product)
    if direct is not None:
        return {"method": "atom_maps", "mapping": direct, "events": _stereo_events(reactant, product, direct)}
    result = rdFMCS.FindMCS(
        [_labelled(reactant), _labelled(product, blocked_product_atoms)],
        atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        matchChiralTag=False,
        timeout=timeout_s,
    )
    if not result.numAtoms:
        return {"method": "mcs", "mapping": {}, "events": _stereo_events(reactant, product, {}), "timed_out": result.canceled}
    query = Chem.MolFromSmarts(result.smartsString)
    r_lab, p_lab = _labelled(reactant), _labelled(product, blocked_product_atoms)
    r_matches = r_lab.GetSubstructMatches(query, uniquify=False, maxMatches=24)
    p_matches = p_lab.GetSubstructMatches(query, uniquify=False, maxMatches=48)
    best: dict[str, Any] | None = None
    for r_match in r_matches[:6]:
        for p_match in p_matches:
            mapping = dict(zip(r_match, p_match))
            events = _stereo_events(reactant, product, mapping)
            cost = _event_cost(events)
            if best is None or cost < best["cost"]:
                best = {"mapping": mapping, "events": events, "cost": cost}
            if cost == (0, 0, 0):
                break
        if best and best["cost"] == (0, 0, 0):
            break
    assert best is not None
    return {
        "method": "mcs",
        "mapping": best["mapping"],
        "events": best["events"],
        "timed_out": bool(result.canceled),
    }


# --------------------------------------------------------------------------
# the audit


def _issue(severity: str, code: str, message: str, **extra: Any) -> dict[str, Any]:
    return {"severity": severity, "code": code, "message": message, **extra}


def _sum_formula(mols: list[Any]) -> Counter:
    total: Counter = Counter()
    for mol in mols:
        total.update(formula_counter(mol))
    return total


def _charge(mols: list[Any]) -> int:
    return sum(a.GetFormalCharge() for m in mols for a in m.GetAtoms())


def _delta(after: Counter, before: Counter) -> dict[str, int]:
    keys = set(after) | set(before)
    return {k: after.get(k, 0) - before.get(k, 0) for k in sorted(keys) if after.get(k, 0) != before.get(k, 0)}


class _Compounds:
    def __init__(self, raw: Any) -> None:
        self.smiles: dict[str, str] = {}
        self.names: dict[str, str] = {}
        self.mols: dict[str, Any] = {}
        self.errors: list[dict[str, Any]] = []
        for key, value in (raw or {}).items():
            key = str(key)
            if isinstance(value, dict):
                smiles = value.get("smiles")
                if value.get("name"):
                    self.names[key] = str(value["name"])
            else:
                smiles = value
            self._add(key, smiles)

    def _add(self, key: str, smiles: Any) -> None:
        try:
            self.mols[key] = parse(str(smiles), keep_maps=True)
            self.smiles[key] = str(smiles)
        except StructureError as exc:
            self.errors.append(_issue("blocking", "unparseable_structure", f"Compound {key}: {exc}", compound=key))

    def resolve(self, ref: Any) -> str | None:
        """Return the compound key for a reference; inline SMILES are registered."""
        ref = str(ref)
        if ref in self.mols:
            return ref
        if any(ch in ref for ch in "()=#[]@") or ref.isalpha() or Chem.MolFromSmiles(ref) is not None:
            try:
                key = f"inline:{ref}"
                if key not in self.mols:
                    self.mols[key] = parse(ref, keep_maps=True)
                    self.smiles[key] = ref
                return key
            except StructureError:
                return None
        return None


def _stripped(mol: Any) -> Any:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return copy


def route_check(route: dict[str, Any], *, mcs_timeout_s: int = 6) -> dict[str, Any]:
    """Audit a route. See the module docstring and ``route_check_example``."""
    if not RDKIT_AVAILABLE:
        return {"status": "unavailable", "warnings": ["RDKit is not available."]}
    if not isinstance(route, dict) or not isinstance(route.get("steps"), list) or not route["steps"]:
        return {
            "status": "error",
            "warnings": ["A route needs 'compounds' (id -> SMILES) and a non-empty 'steps' list."],
        }
    compounds = _Compounds(route.get("compounds"))
    issues: list[dict[str, Any]] = list(compounds.errors)
    step_reports: list[dict[str, Any]] = []
    produced_by: dict[str, int] = {}
    consumed: Counter = Counter()
    depth: dict[str, int] = {}
    group_cache: dict[str, Counter] = {}
    created_total = inverted_total = destroyed_total = 0

    def groups(key: str) -> Counter:
        if key not in group_cache:
            group_cache[key] = functional_group_counts(_stripped(compounds.mols[key]))
        return group_cache[key]

    for position, step in enumerate(route["steps"], start=1):
        label = step.get("step", position)
        report: dict[str, Any] = {"step": label, "conditions": step.get("conditions", ""), "issues": []}
        step_reports.append(report)
        local = report["issues"]
        conditions = str(step.get("conditions") or "")
        reactant_keys: list[str] = []
        for ref in step.get("reactants") or []:
            key = compounds.resolve(ref)
            if key is None:
                local.append(_issue("blocking", "unknown_reactant", f"Reactant {ref!r} is not a defined compound or SMILES."))
            else:
                reactant_keys.append(key)
        product_key = compounds.resolve(step.get("product"))
        if product_key is None:
            local.append(_issue("blocking", "unknown_product", f"Product {step.get('product')!r} is not a defined compound or SMILES."))
        if product_key is None or not reactant_keys:
            if not reactant_keys:
                local.append(_issue("blocking", "no_reactants", "The step lists no usable reactant."))
            continue
        report["reactants"] = [k.removeprefix("inline:") for k in reactant_keys]
        report["product"] = product_key.removeprefix("inline:")
        for key in reactant_keys:
            consumed[key] += 1
        if product_key in produced_by:
            local.append(_issue("review", "product_made_twice", f"Compound {product_key} is also the product of step {produced_by[product_key]}."))
        produced_by[product_key] = label
        depth[product_key] = 1 + max((depth.get(k, 0) for k in reactant_keys), default=0)
        reactant_mols = [compounds.mols[k] for k in reactant_keys]
        product_mol = compounds.mols[product_key]

        # --- atoms ---------------------------------------------------------
        extra_in, extra_out, declared = [], [], False
        for field, bucket in (("reagents", extra_in), ("byproducts", extra_out)):
            for smi in step.get(field) or []:
                declared = True
                try:
                    bucket.append(parse(str(smi)))
                except StructureError as exc:
                    local.append(_issue("blocking", "unparseable_structure", f"{field}: {exc}"))
        before = _sum_formula(reactant_mols)
        after = formula_counter(product_mol)
        delta = _delta(after, before)
        balance: dict[str, Any] = {
            "reactants_formula": hill_formula(before),
            "product_formula": hill_formula(after),
            "delta_product_minus_reactants": delta_text(delta),
        }
        if declared:
            residual = _delta(after + _sum_formula(extra_out), before + _sum_formula(extra_in))
            charge_residual = _charge([product_mol, *extra_out]) - _charge([*reactant_mols, *extra_in])
            balance["declared_reagents_and_byproducts"] = True
            if residual or charge_residual:
                balance["status"] = "does_not_close"
                balance["residual"] = delta_text(residual)
                if charge_residual:
                    balance["charge_residual"] = charge_residual
                complete = bool(step.get("reagents")) and bool(step.get("byproducts"))
                local.append(
                    _issue(
                        "blocking" if complete else "review",
                        "atom_balance_does_not_close",
                        f"With the declared reagents/byproducts the equation is off by {delta_text(residual) or '0'}"
                        + (f" and charge {charge_residual:+d}" if charge_residual else "")
                        + ". An atom was dropped, added or miscounted in a structure of this step.",
                    )
                )
            else:
                balance["status"] = "closed"
        else:
            interpretation = interpret_delta(delta, conditions)
            balance.update(interpretation)
            if interpretation["status"] == "unexplained":
                local.append(
                    _issue(
                        "review",
                        "unexplained_formula_change",
                        f"Product minus reactants = {delta_text(delta)}. This is not a recognised gain/loss; "
                        "confirm which reagent supplies and which byproduct removes these atoms, or fix the structure.",
                    )
                )
            elif interpretation["status"] == "named_but_conditions_do_not_mention_a_matching_reagent":
                local.append(
                    _issue(
                        "review",
                        "formula_change_not_matched_by_conditions",
                        f"Product minus reactants = {delta_text(delta)} ({'; '.join(interpretation['matches'][:2])}), "
                        f"but the conditions ({conditions or 'none given'}) name no reagent that does this. "
                        "Either the conditions are incomplete or a structure has the wrong atom count.",
                    )
                )
            elif interpretation["status"] == "explained_only_as_a_combination":
                local.append(
                    _issue(
                        "info",
                        "formula_change_is_a_combination",
                        f"Product minus reactants = {delta_text(delta)}: only explained as two changes at once "
                        f"({interpretation['possible_combinations'][0]}). Confirm both really happen in this step.",
                    )
                )
            charge_delta = _charge([product_mol]) - _charge(reactant_mols)
            if charge_delta:
                balance["charge_delta"] = charge_delta
        report["atom_balance"] = balance

        # --- mapping and stereo ---------------------------------------------
        covered: set[int] = set()
        stereo_report: dict[str, list[Any]] = {
            "retained": [],
            "inverted": [],
            "epimerized_without_bond_change": [],
            "destroyed": [],
            "lost_definition": [],
            "ambiguous": [],
        }
        methods = set()
        traced_centres: set[int] = set()
        lost_atoms: Counter = Counter()
        mapping_complete = True
        ordered = sorted(
            zip(reactant_keys, reactant_mols), key=lambda item: item[1].GetNumHeavyAtoms(), reverse=True
        )
        for key, reactant in ordered:
            try:
                mapped = map_reactant_to_product(
                    reactant, product_mol, blocked_product_atoms=covered or None, timeout_s=mcs_timeout_s
                )
            except Exception as exc:  # pragma: no cover - defensive: never fail the whole audit
                local.append(_issue("info", "mapping_failed", f"Atom mapping failed for {key}: {exc}"))
                continue
            methods.add(mapped["method"])
            if mapped.get("timed_out"):
                mapping_complete = False
                report.setdefault("notes", []).append(
                    f"MCS search for {key.removeprefix('inline:')} hit the time limit; mapping may be incomplete."
                )
            covered |= set(mapped["mapping"].values())
            lost_atoms.update(
                atom.GetSymbol()
                for atom in reactant.GetAtoms()
                if atom.GetAtomicNum() > 1 and atom.GetIdx() not in mapped["mapping"]
            )
            for bucket, entries in mapped["events"].items():
                for entry in entries:
                    stereo_report[bucket].append({"reactant": key.removeprefix("inline:"), **entry})
                    if "product_atom" in entry:
                        traced_centres.add(entry["product_atom"])
        inverting = bool(_INVERTING_CONDITIONS.search(conditions))
        kept_inverted = []
        for entry in stereo_report["inverted"]:
            if "replaced" in entry:
                kept_inverted.append(entry)
            elif inverting:
                kept_inverted.append({**entry, "replaced": entry["same_element_exchange_possible"]})
            else:
                stereo_report["derivatised_but_inverted"] = stereo_report.get("derivatised_but_inverted", []) + [entry]
        stereo_report["inverted"] = kept_inverted
        retained_despite_inverting = [
            entry
            for entry in stereo_report["retained"]
            if inverting and ("replaced" in entry or "same_element_exchange_possible" in entry)
        ]
        created, created_undefined = [], []
        for element in stereo_elements(product_mol):
            if element["kind"] != "tetrahedral" or element["atom_index"] in traced_centres:
                continue
            entry = {
                "product_atom": element["atom_index"],
                "element": element["element"],
                "substituents": element["substituents"],
                "cip": element["cip"],
            }
            (created if element["specified"] else created_undefined).append(entry)
        open_count, n_consistent = undefined_stereo_count(_stripped(product_mol))
        if not open_count:
            created_undefined = []
        heavy_product = product_mol.GetNumHeavyAtoms()
        report["mapping"] = {
            "method": "+".join(sorted(methods)) or "none",
            "product_heavy_atoms_traced_to_reactants": f"{len(covered)}/{heavy_product}",
            "untraced_product_atoms": dict(
                Counter(
                    product_mol.GetAtomWithIdx(i).GetSymbol()
                    for i in range(product_mol.GetNumAtoms())
                    if i not in covered
                )
            ),
        }
        gained_atoms = Counter(
            product_mol.GetAtomWithIdx(i).GetSymbol() for i in range(product_mol.GetNumAtoms()) if i not in covered
        )
        report["mapping"]["reactant_atoms_not_found_in_product"] = dict(lost_atoms)
        # Every product atom that cannot be traced to a numbered reactant must be
        # new according to the formula. If more atoms of an element are untraced
        # than the formula gains, an atom moved or was exchanged for one of the
        # same element: a rearrangement, a reagent-for-substrate swap, or a
        # substituent drawn in the wrong place (which a formula alone cannot see).
        declared_elements = {a.GetSymbol() for m in (*extra_in, *extra_out) for a in m.GetAtoms()}
        moved = {}
        for symbol, count in gained_atoms.items():
            if symbol == "H" or symbol in declared_elements:
                continue
            excess = count - max(0, delta.get(symbol, 0))
            if excess > 0:
                moved[symbol] = excess
        if moved and mapping_complete and "atom_maps" not in methods:
            local.append(
                _issue(
                    "review",
                    "atoms_relocated",
                    "The atom mapping cannot follow "
                    + ", ".join(f"{count} {symbol}" for symbol, count in sorted(moved.items()))
                    + " into the product, and the formula change does not account for them as new atoms. "
                    "Either the skeleton rearranges or a reagent atom replaces a substrate atom of the same "
                    "element in this step (declare reagents/byproducts as SMILES, or give atom-mapped SMILES), "
                    "or a chain length or substituent position differs from the reactant.",
                    elements=moved,
                )
            )
        report["stereo"] = {
            "n_retained": len(stereo_report["retained"]),
            "inverted_with_substituent_replaced": stereo_report["inverted"],
            "inverted_where_only_heteroatom_derivatised": stereo_report.get("derivatised_but_inverted", []),
            "epimerized_without_bond_change": stereo_report["epimerized_without_bond_change"],
            "destroyed": stereo_report["destroyed"],
            "lost_definition": stereo_report["lost_definition"],
            "created_defined": created,
            "created_undefined": created_undefined,
            "ambiguous": stereo_report["ambiguous"],
            "product_isomers_consistent_with_drawing": n_consistent,
        }
        created_total += len(created)
        inverted_total += len(stereo_report["inverted"])
        destroyed_total += len(stereo_report["destroyed"])
        if stereo_report.get("derivatised_but_inverted"):
            local.append(
                _issue(
                    "review",
                    "inversion_where_only_the_heteroatom_was_derivatised",
                    f"{len(stereo_report['derivatised_but_inverted'])} centre(s) are drawn inverted although the "
                    "conditions only derivatise the O/N attached to them (acylation, sulfonylation, silylation and "
                    "alkylation of the heteroatom do not break its bond to carbon). Either a wedge was flipped or "
                    "this is a Mitsunobu/SN2 step and the conditions should say so.",
                    elements=stereo_report["derivatised_but_inverted"],
                )
            )
        if retained_despite_inverting:
            local.append(
                _issue(
                    "review",
                    "retention_under_inverting_conditions",
                    f"{len(retained_despite_inverting)} centre(s) keep their configuration although the conditions "
                    "(Mitsunobu, azide/cyanide/halide/carboxylate displacement, Appel) substitute with inversion. "
                    "Check the wedge on the product, or explain the retention (double inversion, participation).",
                    elements=retained_despite_inverting,
                )
            )
        if stereo_report["epimerized_without_bond_change"]:
            local.append(
                _issue(
                    "review",
                    "configuration_changed_without_bond_change",
                    f"{len(stereo_report['epimerized_without_bond_change'])} stereo element(s) have the opposite "
                    "configuration in the product although none of their bonds changed. Unless this step is an "
                    "intended epimerisation/isomerisation, a wedge was flipped when the product was written.",
                    elements=stereo_report["epimerized_without_bond_change"],
                )
            )
        if stereo_report["inverted"]:
            local.append(
                _issue(
                    "info",
                    "inversion_at_substituted_centre",
                    f"{len(stereo_report['inverted'])} centre(s) are inverted where a substituent was replaced "
                    "(SN2/Mitsunobu-type). Confirm the mechanism really inverts; double inversion or SN1 would not.",
                    elements=stereo_report["inverted"],
                )
            )
        retained_replaced = [e for e in stereo_report["retained"] if "replaced" in e and not inverting]
        if retained_replaced:
            local.append(
                _issue(
                    "info",
                    "retention_at_substituted_centre",
                    f"{len(retained_replaced)} centre(s) keep their configuration although a substituent was "
                    "replaced at the centre. That needs a reason (double inversion, neighbouring-group "
                    "participation, the bond to the centre never breaking, or a retentive mechanism).",
                    elements=retained_replaced,
                )
            )
        if stereo_report["lost_definition"] or (open_count and created_undefined):
            local.append(
                _issue(
                    "review",
                    "undefined_stereo_in_product",
                    f"The product leaves {open_count} stereo element(s) undefined"
                    + (f" ({n_consistent} isomers fit the drawing)" if n_consistent else "")
                    + ". Draw the isomer that forms, or state that it is an inconsequential mixture removed later.",
                )
            )
        elif open_count:
            local.append(
                _issue(
                    "review",
                    "undefined_stereo_in_product",
                    f"The product has {open_count} undefined stereo element(s) inherited from its precursors.",
                )
            )

        # --- functional groups -----------------------------------------------
        before_groups: Counter = Counter()
        for key in reactant_keys:
            before_groups.update(groups(key))
        after_groups = groups(product_key)
        changed = {
            name: after_groups.get(name, 0) - before_groups.get(name, 0)
            for name in sorted(set(before_groups) | set(after_groups))
            if after_groups.get(name, 0) != before_groups.get(name, 0)
        }
        # A group present on both sides survived at least min(before, after) times.
        bystanders = {
            name: min(count, before_groups[name])
            for name, count in sorted(after_groups.items())
            if before_groups.get(name, 0) > 0
        }
        report["functional_groups"] = {"changed": changed, "bystanders_that_must_survive": bystanders}
        partially_converted: set[str] = set()
        for name in sorted(_SELECTIVITY_GROUPS):
            n_before, n_after = before_groups.get(name, 0), after_groups.get(name, 0)
            if n_before > n_after > 0:
                partially_converted.add(name)
                local.append(
                    _issue(
                        "review",
                        "selectivity_among_equal_groups",
                        f"{n_before - n_after} of {n_before} {name} group(s) reacted and {n_after} did not. "
                        "State what differentiates them (sterics, electronics, directing group, stoichiometry).",
                    )
                )
        checks = []
        for rule_name, pattern, risks in _RULES_COMPILED:
            if not pattern.search(conditions):
                continue
            for name, why in risks.items():
                # the selectivity question already covers a partly converted group
                if name in bystanders and name not in partially_converted:
                    checks.append({"reagent_class": rule_name, "group": name, "count": bystanders[name], "risk": why})
        if checks:
            report["compatibility_checks"] = checks
            summary = "; ".join(f"{c['group']} x{c['count']} vs {c['reagent_class']} ({c['risk']})" for c in checks[:6])
            local.append(
                _issue(
                    "review",
                    "bystander_group_vs_reagent",
                    "Groups drawn as unchanged that this reagent class is known to attack: "
                    + summary
                    + ("; ..." if len(checks) > 6 else "")
                    + ". Justify survival for this substrate or change the reagent/order.",
                )
            )
        if not conditions.strip():
            local.append(_issue("review", "no_conditions", "No conditions are given, so reagent compatibility was not checked."))

    # --- route level --------------------------------------------------------
    for report in step_reports:
        for item in report["issues"]:
            issues.append({"step": report["step"], **item})
    products = list(produced_by)
    final_key = None
    if products:
        unused = [k for k in products if consumed.get(k, 0) == 0]
        final_key = unused[-1] if unused else products[-1]
        for key in unused[:-1]:
            issues.append(
                _issue("review", "dead_end_intermediate", f"Compound {key} is made in step {produced_by[key]} but never used.")
            )
    starting = [k for k in compounds.mols if k not in produced_by and consumed.get(k, 0)]
    starting_report = []
    for key in starting:
        mol = _stripped(compounds.mols[key])
        centres = sum(1 for e in stereo_elements(mol) if e["kind"] == "tetrahedral" and e["specified"])
        heavy = mol.GetNumHeavyAtoms()
        entry = {
            "id": key.removeprefix("inline:"),
            "smiles": canonical(mol),
            "heavy_atoms": heavy,
            "mol_weight": round(Descriptors.MolWt(mol), 1),
            "defined_stereocenters": centres,
        }
        if key.startswith("inline:"):
            entry["note"] = "written inline in a step, not as a numbered compound"
        if centres >= 3 or heavy >= 25:
            entry["verify_availability"] = True
        starting_report.append(entry)
    flagged = [e["id"] for e in starting_report if e.get("verify_availability") and "note" not in e]
    if flagged:
        issues.append(
            _issue(
                "review",
                "starting_material_availability",
                f"Starting material(s) {', '.join(flagged)} are large or stereochemically rich. Check that they are "
                "really purchasable (and at what price) or add the steps that make them.",
            )
        )
    result: dict[str, Any] = {
        "n_steps": len(route["steps"]),
        "longest_linear_sequence": depth.get(final_key, 0) if final_key else 0,
        "final_product": final_key.removeprefix("inline:") if final_key else None,
        "starting_materials": starting_report,
        "stereo_summary": {
            "defined_centres_in_starting_materials": sum(e["defined_stereocenters"] for e in starting_report),
            "centres_created_in_steps": created_total,
            "centres_inverted_in_steps": inverted_total,
            "centres_destroyed_in_steps": destroyed_total,
        },
    }
    if final_key:
        final_groups = groups(final_key)
        left = {k: v for k, v in sorted(final_groups.items()) if k in PROTECTING_GROUPS}
        if left:
            result["protecting_groups_in_final_product"] = left
    target = route.get("target")
    if target and final_key:
        try:
            comparison = compare_structures(canonical(_stripped(compounds.mols[final_key])), str(target))
        except StructureError as exc:
            issues.append(_issue("blocking", "unparseable_structure", f"Target: {exc}"))
            comparison = None
        if comparison is not None:
            result["final_vs_target"] = comparison
            relation = comparison["relation"]
            if relation != "identical":
                detail = ""
                if relation == "different_formula":
                    detail = f" (target minus final product: {comparison['formula_delta_b_minus_a']})"
                elif relation == "diastereomers":
                    detail = f" ({comparison.get('n_differing', '?')} stereo element(s) differ)"
                elif relation == "same_constitution_stereo_undefined":
                    detail = (
                        f" (undefined elements: final {comparison['n_undefined_a']}, target "
                        f"{comparison['n_undefined_b']})"
                    )
                issues.append(
                    _issue(
                        "blocking",
                        "final_product_is_not_the_target",
                        f"The last product and the target are {relation.replace('_', ' ')}{detail}.",
                    )
                )
    elif not target:
        issues.append(_issue("review", "no_target", "No target was given, so the end point was not verified."))
    blocking = [i for i in issues if i["severity"] == "blocking"]
    review = [i for i in issues if i["severity"] == "review"]
    result.update(
        {
            "status": "blocked" if blocking else ("review" if review else "consistent"),
            "n_blocking": len(blocking),
            "n_review": len(review),
            "issues": issues,
            "steps": [
                {k: v for k, v in report.items() if k != "issues"} | {"issues": report["issues"]}
                for report in step_reports
            ],
            "scope": (
                "Bookkeeping audit: atoms, configuration, functional-group inventory and target identity. "
                "'review' items come from heuristic atom mapping and reagent keyword rules and may be false "
                "alarms. A consistent audit does not show that any step works experimentally."
            ),
        }
    )
    return result


def route_check_example() -> dict[str, Any]:
    """A small valid route showing every field of the route format."""
    return {
        "target": "C[C@H](N)CC(=O)OC",
        "compounds": {
            "1": {"smiles": "C[C@@H](O)CC(=O)OC", "name": "methyl (R)-3-hydroxybutanoate"},
            "2": {"smiles": "C[C@@H](OS(C)(=O)=O)CC(=O)OC", "name": "mesylate"},
            "3": {"smiles": "C[C@H](N=[N+]=[N-])CC(=O)OC", "name": "(S)-azide"},
            "4": {"smiles": "C[C@H](N)CC(=O)OC", "name": "methyl (S)-3-aminobutanoate"},
        },
        "steps": [
            {
                "step": 1,
                "reactants": ["1"],
                "product": "2",
                "conditions": "MsCl, Et3N, CH2Cl2, 0 °C",
                "why": "O-sulfonylation; the C-O bond is untouched, so the centre is retained.",
            },
            {
                "step": 2,
                "reactants": ["2"],
                "product": "3",
                "conditions": "NaN3, DMF, 60 °C",
                "reagents": ["[N-]=[N+]=[N-]"],
                "byproducts": ["CS(=O)(=O)[O-]"],
                "why": "SN2 displacement with inversion.",
            },
            {
                "step": 3,
                "reactants": ["3"],
                "product": "4",
                "conditions": "PPh3, THF/H2O, rt",
                "why": "Staudinger reduction; nitrogen stays on the same carbon.",
            },
        ],
    }


def route_to_figure_spec(
    route: dict[str, Any], *, steps_per_row: int | None = None, max_width: int = 1200
) -> dict[str, Any]:
    """Turn a route into a ``chem_route_figure`` spec so every compound gets drawn.

    Each compound is captioned with its number in bold and its name, if any,
    underneath. By default the steps form one row that the renderer wraps at
    arrows into rows of similar width (``max_width`` pixels); a step whose
    reactants do not include the previous product starts a new row. Pass
    ``steps_per_row`` for fixed rows instead, each starting with its first
    reactant.
    """
    compounds = _Compounds(route.get("compounds"))
    rows: list[dict[str, Any]] = []
    items: list[dict[str, Any]] = []
    previous_product: str | None = None
    count = 0

    def molecule(key: str) -> dict[str, Any]:
        return {
            "id": f"cmpd_{len(rows)}_{len(items)}",
            "label": key if not key.startswith("inline:") else "",
            "code": compounds.names.get(key, ""),
            "smiles": canonical(_stripped(compounds.mols[key])),
        }

    for step in route.get("steps") or []:
        reactant_keys = [k for k in (compounds.resolve(r) for r in step.get("reactants") or []) if k]
        product_key = compounds.resolve(step.get("product"))
        if not reactant_keys or product_key is None:
            continue
        chain_breaks = previous_product is not None and previous_product not in reactant_keys
        if (steps_per_row and count and count % steps_per_row == 0) or (items and chain_breaks):
            rows.append({"items": items})
            items = []
            previous_product = None
        main = [k for k in reactant_keys if k != previous_product]
        if previous_product is None or previous_product not in reactant_keys:
            first, rest = reactant_keys[0], reactant_keys[1:]
            items.append(molecule(first))
        else:
            rest = main
        for key in rest:
            items.append({"type": "plus"})
            items.append(molecule(key))
        above, below = split_conditions(str(step.get("conditions") or ""))
        items.append({"type": "arrow", "label": above, "sublabel": below})
        items.append(molecule(product_key))
        previous_product = product_key
        count += 1
    if items:
        rows.append({"items": items})
    spec: dict[str, Any] = {"title": route.get("title", ""), "style_preset": "acs", "rows": rows}
    if not steps_per_row:
        spec["max_width"] = max_width
    return spec


_SOLVENTS = {
    "ch2cl2", "dcm", "chcl3", "ccl4", "dce", "1,2-dce", "thf", "2-methf", "et2o", "ether", "diethyl ether",
    "dioxane", "1,4-dioxane", "dme", "meoh", "etoh", "iproh", "i-proh", "ipa", "tbuoh", "t-buoh", "n-buoh",
    "h2o", "water", "mecn", "ch3cn", "acetonitrile", "dmf", "dma", "dmac", "dmso", "nmp", "hmpa", "toluene",
    "phme", "benzene", "phh", "xylene", "xylenes", "hexane", "hexanes", "pentane", "heptane", "cyclohexane",
    "etoac", "acetone", "tfe", "hfip", "meno2", "ch3no2", "phcf3", "chlorobenzene", "phcl", "mtbe", "tbme",
    "cpme", "diglyme",
}
_SETTING = re.compile(
    r"°\s*C|\b(?:rt|r\.t\.|reflux|overnight|o/n|room temperature|heat|microwave)\b|Δ"
    r"|^-?\d+(?:\.\d+)?\s*(?:h|hr|hrs|min|d|days?|s|K)$"
    r"|^\d+(?:\.\d+)?\s*(?:atm|bar|psi)$|^(?:N2|Ar|argon|nitrogen|air)$",
    re.IGNORECASE,
)
_RATIO = re.compile(r"\s*\(?\d+(?:\.\d+)?\s*:\s*\d+(?:\.\d+)?\)?$")


def split_conditions(conditions: str) -> tuple[str, str]:
    """Split a conditions string into the text above and below a scheme arrow.

    Reagents and catalysts go above; solvents (including mixtures such as
    THF/H2O 3:1), temperatures, times and atmospheres go below, in their
    original order. Parts are separated at a comma followed by a space, so
    locants such as 2,2-dimethoxypropane stay whole.
    """
    plain = conditions.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
    above: list[str] = []
    below: list[str] = []
    for part in (p.strip() for p in re.split(r",\s+", plain.strip())):
        if not part:
            continue
        solvents = [s.strip().lower() for s in re.split(r"[/–]", _RATIO.sub("", part))]
        if _SETTING.search(part) or all(s in _SOLVENTS for s in solvents):
            below.append(part)
        else:
            above.append(part)
    return ", ".join(above), ", ".join(below)
