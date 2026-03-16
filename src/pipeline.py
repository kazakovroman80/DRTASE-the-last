"""Pipeline entry points for generating candidate mutant structures."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List

from src.mutation_generator import MutationTarget, mutate_residue_pair
from src.structure_loader import load_pdb, write_pdb


INPUT_PDB = Path("data/input/mutant_dimer_background_minimized.pdb")
OUTPUT_STRUCTURES_DIR = Path("results/structures")
OUTPUT_TABLE_DIR = Path("results/tables")


def candidate_mutations() -> List[MutationTarget]:
    """Return the mutation list for this project."""

    return [
        MutationTarget("Y69F", 69, 269, "Symmetric Tyr->Phe mutation at the dimer interface region."),
        MutationTarget("N70Q", 70, 270, "Symmetric Asn->Gln mutation; preserves amide chemistry with longer side chain."),
        MutationTarget("N71Q", 71, 271, "Symmetric Asn->Gln mutation near position 70 with conservative polarity."),
        MutationTarget("Y86F", 86, 286, "Symmetric Tyr->Phe mutation removing hydroxyl while retaining aromatic ring."),
        MutationTarget("Y101F", 101, 301, "Symmetric Tyr->Phe mutation at position 101 in both protomers."),
    ]


def _target_residue_name(mutation_label: str) -> str:
    return "PHE" if mutation_label.endswith("F") else "GLN"


def run_pipeline() -> None:
    """Generate all mutant PDB files and the candidate summary table."""

    structure = load_pdb(INPUT_PDB)
    OUTPUT_STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for target in candidate_mutations():
        mutant_residue = _target_residue_name(target.mutation)
        mutated_structure = mutate_residue_pair(
            structure=structure,
            chain_a_residue=target.chain_a_residue,
            chain_b_residue=target.chain_b_residue,
            mutant_name=mutant_residue,
        )

        out_name = f"{target.mutation}_mutant_dimer.pdb"
        out_path = OUTPUT_STRUCTURES_DIR / out_name
        write_pdb(out_path, mutated_structure)

        rows.append(
            {
                "mutation": target.mutation,
                "chain_A_residue": target.chain_a_residue,
                "chain_B_residue": target.chain_b_residue,
                "description": target.description,
            }
        )

    table_path = OUTPUT_TABLE_DIR / "candidate_mutations.csv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mutation", "chain_A_residue", "chain_B_residue", "description"],
        )
        writer.writeheader()
        writer.writerows(rows)
