"""XYZ file reader/writer (single- and multi-frame)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core.molecule import Molecule
from ..core.results import ParseResult
from .errors import FileFormatError


def read(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    frames: list[Molecule] = []
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue
        try:
            natoms = int(stripped.split()[0])
        except ValueError as exc:
            raise FileFormatError(
                f"{path}: expected an atom count at line {i + 1}, got {stripped!r}"
            ) from exc
        if i + 1 + natoms >= len(lines) + 1:
            raise FileFormatError(f"{path}: truncated frame at line {i + 1}")
        title = lines[i + 1].strip() if i + 1 < len(lines) else ""
        symbols, coords = [], []
        for row in lines[i + 2:i + 2 + natoms]:
            parts = row.split()
            if len(parts) < 4:
                raise FileFormatError(f"{path}: malformed atom line {row!r}")
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        frames.append(Molecule(symbols, np.array(coords), title=title))
        i += 2 + natoms
    if not frames:
        raise FileFormatError(f"{path}: no frames found")
    return ParseResult(frames=frames, program="xyz", source=str(path))


def write(path, molecule: Molecule) -> None:
    write_frames(path, [molecule])


def write_frames(path, molecules: list[Molecule]) -> None:
    with open(path, "w") as fh:
        for mol in molecules:
            fh.write(f"{mol.natoms}\n{mol.title or mol.formula()}\n")
            for sym, (x, y, z) in zip(mol.symbols, mol.coords):
                fh.write(f"{sym:<3s} {x:15.8f} {y:15.8f} {z:15.8f}\n")
