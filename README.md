# Homodimer Mutant Structure Pipeline

This project generates symmetric mutant models for a homodimeric enzyme from a minimized starting PDB file.

## Input

- `data/input/mutant_dimer_background_minimized.pdb`

The input dimer already contains background experimental mutations and is kept intact except for the requested candidate substitutions.

## Generated outputs

Running the pipeline creates:

- `results/structures/Y69F_mutant_dimer.pdb`
- `results/structures/N70Q_mutant_dimer.pdb`
- `results/structures/N71Q_mutant_dimer.pdb`
- `results/structures/Y86F_mutant_dimer.pdb`
- `results/structures/Y101F_mutant_dimer.pdb`
- `results/tables/candidate_mutations.csv`

## Project layout

```text
src/
    structure_loader.py
    mutation_generator.py
    pipeline.py
scripts/
    run_pipeline.py
requirements.txt
README.md
```

## How to run

```bash
python scripts/run_pipeline.py
```

## Notes on mutation behavior

- Mutations are applied symmetrically in chain A and chain B.
- Original residue numbering is preserved.
- Coordinates are retained for existing atoms whenever possible.
- `TYR -> PHE`: hydroxyl atom is removed.
- `ASN -> GLN`: existing atoms are reused where possible, and one additional side-chain atom (`CG`) is placed using a simple midpoint geometry.
