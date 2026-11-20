"""One import for free-form chemistry code.

    from codex_organic_chem import kit

    kit.compare("C[C@H](O)CC", "C[C@@H](O)CC")["relation"]      # 'enantiomers'
    kit.route_check(route)["issues"]                             # step-by-step audit
    kit.face(smiles, atom=8)["verdict"]                          # which face is open
    kit.cis_trans(smiles, 4, 6)["relationship"]                  # two ring substituents: cis or trans
    kit.substitution_products("CC1CCCCC1", "Cl")                 # count, do not guess

Every function returns plain dicts and raises ``kit.StructureError`` on
unparseable input. Use these together with ordinary RDKit code: write whatever
check the question needs instead of reasoning about a count, a formula or a
3D relationship in your head.
"""

from __future__ import annotations

from .model_kit import (
    atom_table,
    build_ensemble,
    dihedral_profile as dihedral,
    distance_profile as distance,
    embed_check,
    face_accessibility as face,
    isomer_energies,
    planarity,
    ring_relationship as cis_trans,
)
from .route_check import (
    functional_groups,
    interpret_delta,
    map_reactant_to_product,
    route_check,
    route_check_example,
    route_to_figure_spec,
)
from .structure_compare import (
    StructureError,
    compare_structures as compare,
    enumerate_stereoisomers,
    stereo_audit,
    substitution_products,
    symmetry_environments,
)

__all__ = [
    "StructureError",
    "atom_table",
    "build_ensemble",
    "cis_trans",
    "compare",
    "dihedral",
    "distance",
    "embed_check",
    "enumerate_stereoisomers",
    "face",
    "functional_groups",
    "interpret_delta",
    "isomer_energies",
    "map_reactant_to_product",
    "planarity",
    "route_check",
    "route_check_example",
    "route_to_figure_spec",
    "stereo_audit",
    "substitution_products",
    "symmetry_environments",
]
