"""Utilities for reading and writing simple PDB structures.

The parser here intentionally focuses on ATOM/HETATM records and keeps the
original line formatting where possible. This keeps the code beginner-friendly
while still being robust enough for mutation workflows.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class AtomRecord:
    """Represents one ATOM/HETATM record in a PDB file."""

    record_type: str
    serial: int
    atom_name: str
    alt_loc: str
    residue_name: str
    chain_id: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float
    occupancy: float | None
    b_factor: float | None
    element: str
    charge: str

    @classmethod
    def from_pdb_line(cls, line: str) -> "AtomRecord":
        """Parse one fixed-width PDB ATOM/HETATM line."""

        def safe_float(text: str) -> float | None:
            text = text.strip()
            return float(text) if text else None

        return cls(
            record_type=line[0:6].strip(),
            serial=int(line[6:11]),
            atom_name=line[12:16],
            alt_loc=line[16:17],
            residue_name=line[17:20].strip(),
            chain_id=line[21:22],
            residue_number=int(line[22:26]),
            insertion_code=line[26:27],
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            occupancy=safe_float(line[54:60]) if len(line) >= 60 else None,
            b_factor=safe_float(line[60:66]) if len(line) >= 66 else None,
            element=line[76:78].strip() if len(line) >= 78 else "",
            charge=line[78:80].strip() if len(line) >= 80 else "",
        )

    def to_pdb_line(self) -> str:
        """Serialize the record back to a valid PDB ATOM/HETATM line."""

        occupancy = self.occupancy if self.occupancy is not None else 1.00
        b_factor = self.b_factor if self.b_factor is not None else 0.00

        # Derive element from atom name if missing.
        element = self.element.strip() or "".join(ch for ch in self.atom_name if ch.isalpha()).strip()[:2].rjust(2)

        return (
            f"{self.record_type:<6}{self.serial:>5d} "
            f"{self.atom_name:<4}{self.alt_loc:1}"
            f"{self.residue_name:>3} {self.chain_id:1}"
            f"{self.residue_number:>4d}{self.insertion_code:1}   "
            f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}"
            f"{occupancy:>6.2f}{b_factor:>6.2f}          "
            f"{element:>2}{self.charge:>2}"
        )


@dataclass
class StructureFile:
    """A parsed structure with atom records and trailing non-atom lines."""

    atoms: List[AtomRecord]
    trailing_lines: List[str]


def load_pdb(path: str | Path) -> StructureFile:
    """Load a PDB file into a simple structure container."""

    atoms: List[AtomRecord] = []
    trailing: List[str] = []

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(AtomRecord.from_pdb_line(line))
            elif line.startswith("END"):
                # Keep END at write time.
                continue
            else:
                trailing.append(line)

    return StructureFile(atoms=atoms, trailing_lines=trailing)


def write_pdb(path: str | Path, structure: StructureFile) -> None:
    """Write a structure back to disk in PDB format."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for serial, atom in enumerate(structure.atoms, start=1):
            atom.serial = serial
            handle.write(atom.to_pdb_line() + "\n")

        # Keep non-coordinate metadata except original CONECT records because
        # atom content can change during mutations.
        for line in structure.trailing_lines:
            if line.startswith("CONECT"):
                continue
            handle.write(line + "\n")

        handle.write("END\n")
