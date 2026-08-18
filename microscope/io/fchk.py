"""Gaussian formatted checkpoint (.fchk) reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core import elements
from ..core.molecule import Molecule
from ..core.results import ParseResult
from .errors import FileFormatError

BOHR_TO_ANGSTROM = 0.52917721092


def _read_array(lines: list[str], i: int, count: int, cast) -> tuple[list, int]:
    values: list = []
    i += 1
    while i < len(lines) and len(values) < count:
        values.extend(cast(tok) for tok in lines[i].split())
        i += 1
    return values[:count], i


def read(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    n = len(lines)
    zs: list[int] = []
    coords: np.ndarray | None = None
    charge, mult = 0, 1

    i = 0
    while i < n:
        line = lines[i]
        if line.startswith("Atomic numbers"):
            count = int(line.split("N=")[1])
            zs, i = _read_array(lines, i, count, int)
            continue
        if line.startswith("Current cartesian coordinates"):
            count = int(line.split("N=")[1])
            vals, i = _read_array(lines, i, count, float)
            coords = np.array(vals).reshape(-1, 3) * BOHR_TO_ANGSTROM
            continue
        if line.startswith("Charge") and "I" in line:
            charge = int(line.split()[-1])
        elif line.startswith("Multiplicity") and "I" in line:
            mult = int(line.split()[-1])
        i += 1

    if not zs or coords is None:
        raise FileFormatError(f"{path}: no geometry found in fchk file")
    if len(zs) != len(coords):
        raise FileFormatError(f"{path}: atom count mismatch in fchk file")

    title = lines[0].strip() if lines else ""
    mol = Molecule(
        symbols=[elements.SYMBOLS[z] if 0 < z < len(elements.SYMBOLS) else "X" for z in zs],
        coords=coords,
        charge=charge,
        multiplicity=mult,
        title=title,
    )
    return ParseResult(frames=[mol], program="Gaussian fchk", source=str(path))
