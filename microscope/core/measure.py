"""Naming what a selection of atoms measures.

A distance, angle or torsion, and the label a chemist would write for it. It
is the geometry and the wording only — where the text is drawn is the
viewer's business — so it can be checked without a display.
"""

from __future__ import annotations

from . import geometry
from .molecule import Molecule

# What a selection of this many atoms measures, and the symbol for it.
KINDS = {2: "d", 3: "∠", 4: "φ"}


def value(molecule: Molecule, atoms: list[int]) -> str:
    """The measurement between 2, 3 or 4 atoms, formatted with its unit."""
    pts = [molecule.coords[i] for i in atoms]
    if len(atoms) == 2:
        return f"{geometry.distance(*pts):.3f} Å"
    if len(atoms) == 3:
        return f"{geometry.angle(*pts):.1f}°"
    return f"{geometry.dihedral(*pts):.1f}°"


def tag(molecule: Molecule, atom: int) -> str:
    """How an atom is named in the viewer: element plus its 1-based number."""
    return f"{molecule.symbols[atom]}{atom + 1}"


def describe(molecule: Molecule | None, atoms: list[int]) -> str:
    """The status-bar line for a selection, empty when there is nothing to say.

    One atom names itself, two to four measure something, and more than that
    is a count — there is no single number for five atoms.
    """
    if molecule is None or not atoms:
        return ""
    if len(atoms) == 1:
        return f"{tag(molecule, atoms[0])} selected"
    if len(atoms) > 4:
        return f"{len(atoms)} atoms selected"
    names = "–".join(tag(molecule, i) for i in atoms)
    return f"{KINDS[len(atoms)]}({names}) = {value(molecule, atoms)}"
