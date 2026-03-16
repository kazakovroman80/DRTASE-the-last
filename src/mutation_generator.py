"""Mutation logic for symmetric homodimer mutations."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

from src.structure_loader import AtomRecord, StructureFile


@dataclass(frozen=True)
class MutationTarget:
    """One symmetric mutation to apply to chain A and chain B."""

    mutation: str
    chain_a_residue: int
    chain_b_residue: int
    description: str


# Allowed atom names for the target amino acids.
RESIDUE_ATOMS: Dict[str, List[str]] = {
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
}


def _atom_key(atom_name: str) -> str:
    return atom_name.strip()


def _build_residue_index(atoms: Iterable[AtomRecord]) -> Dict[Tuple[str, int], List[int]]:
    """Map (chain, residue_number) to atom list indices."""

    index: Dict[Tuple[str, int], List[int]] = {}
    for i, atom in enumerate(atoms):
        key = (atom.chain_id.strip(), atom.residue_number)
        index.setdefault(key, []).append(i)
    return index


def _mutate_tyr_to_phe(residue_atoms: List[AtomRecord]) -> List[AtomRecord]:
    """Mutate TYR-like residue atoms into PHE by dropping hydroxyl atom."""

    kept: List[AtomRecord] = []
    for atom in residue_atoms:
        atom_name = _atom_key(atom.atom_name)
        if atom_name not in RESIDUE_ATOMS["PHE"]:
            continue
        new_atom = deepcopy(atom)
        new_atom.residue_name = "PHE"
        kept.append(new_atom)
    return kept


def _mutate_asn_to_gln(residue_atoms: List[AtomRecord]) -> List[AtomRecord]:
    """Mutate ASN-like residue atoms into GLN.

    Strategy:
    - Keep backbone and CB coordinates unchanged.
    - Reuse ASN CG as GLN CD, OD1 as OE1, ND2 as NE2.
    - Create GLN CG midway between CB and CD to keep geometry simple.
    """

    by_name = {_atom_key(atom.atom_name): deepcopy(atom) for atom in residue_atoms}

    required = ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"]
    missing = [name for name in required if name not in by_name]
    if missing:
        raise ValueError(f"Cannot mutate ASN->GLN, missing atoms: {missing}")

    output: Dict[str, AtomRecord] = {}

    for name in ["N", "CA", "C", "O", "CB"]:
        atom = by_name[name]
        atom.residue_name = "GLN"
        output[name] = atom

    # Reuse existing coordinates for the terminal amide.
    cd = by_name["CG"]
    cd.residue_name = "GLN"
    cd.atom_name = " CD "
    output["CD"] = cd

    oe1 = by_name["OD1"]
    oe1.residue_name = "GLN"
    oe1.atom_name = " OE1"
    output["OE1"] = oe1

    ne2 = by_name["ND2"]
    ne2.residue_name = "GLN"
    ne2.atom_name = " NE2"
    output["NE2"] = ne2

    # New CG placed halfway between CB and CD.
    cb = output["CB"]
    cg = deepcopy(cb)
    cg.residue_name = "GLN"
    cg.atom_name = " CG "
    cg.x = (cb.x + cd.x) / 2.0
    cg.y = (cb.y + cd.y) / 2.0
    cg.z = (cb.z + cd.z) / 2.0
    output["CG"] = cg

    ordered: List[AtomRecord] = []
    for atom_name in RESIDUE_ATOMS["GLN"]:
        ordered.append(output[atom_name])
    return ordered


def mutate_residue_pair(
    structure: StructureFile,
    chain_a_residue: int,
    chain_b_residue: int,
    mutant_name: str,
) -> StructureFile:
    """Apply one symmetric mutation and return a new structure object."""

    new_structure = deepcopy(structure)
    residue_index = _build_residue_index(new_structure.atoms)

    targets = [("A", chain_a_residue), ("B", chain_b_residue)]

    replacements: Dict[Tuple[str, int], List[AtomRecord]] = {}
    for key in targets:
        if key not in residue_index:
            raise ValueError(f"Residue not found: chain {key[0]} residue {key[1]}")

        residue_atoms = [new_structure.atoms[i] for i in residue_index[key]]
        if mutant_name == "PHE":
            replacements[key] = _mutate_tyr_to_phe(residue_atoms)
        elif mutant_name == "GLN":
            replacements[key] = _mutate_asn_to_gln(residue_atoms)
        else:
            raise ValueError(f"Unsupported mutant residue: {mutant_name}")

    # Rebuild atom list with replacements inserted at original residue positions.
    replaced_keys = set(targets)
    output_atoms: List[AtomRecord] = []
    consumed = set()

    for atom in new_structure.atoms:
        key = (atom.chain_id.strip(), atom.residue_number)
        if key in replaced_keys:
            if key not in consumed:
                output_atoms.extend(replacements[key])
                consumed.add(key)
            continue
        output_atoms.append(atom)

    new_structure.atoms = output_atoms
    return new_structure
