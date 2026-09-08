"""MDL molfile (V2000) writer.

Every chemical drawing program reads a molfile, ChemDraw included, so this is
the interchange format for handing a perceived Lewis structure to an editor.
Two flavours are written from the same block builder: the plain 3-D structure
(``write``), and the flat drawing exactly as it appears on screen
(``write_2d``), so a figure can be finished in ChemDraw at the angle it was
composed at here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core.molecule import Molecule

TARGET_BOND_LENGTH = 1.5      # what drawing programs expect of 2-D coordinates
_CHARGE_CODES = {3: 1, 2: 2, 1: 3, -1: 5, -2: 6, -3: 7}   # legacy atom-block field


def write(path, molecule: Molecule) -> None:
    """Write the 3-D structure, single bonds only (no orders are perceived)."""
    bonds = molecule.bonds
    if bonds is None:
        bonds = molecule.perceive_bonds()
    orders = np.ones(len(bonds), dtype=int)
    Path(path).write_text(_block(
        molecule.symbols, molecule.coords, np.asarray(bonds).reshape(-1, 2),
        orders, title=molecule.title or Path(path).stem))


def write_2d(path, structure, coords2d, include=None, title: str = "") -> None:
    """Write the flat drawing: *coords2d* are the screen positions of every
    atom and *include* lists the ones actually drawn (the rest are hydrogens
    folded into a label, and are written as implicit hydrogen counts)."""
    symbols, coords, bonds, orders, hydrogens = _flatten(
        structure, coords2d, include)
    name = title or structure.molecule.title or Path(path).stem
    if structure.bracketed:
        # V2000 charges belong to atoms; one that sits on the species as a
        # whole has nowhere to go but the title line
        name = f"{name} ({_species_note(structure)})"
    Path(path).write_text(_block(symbols, coords, bonds, orders,
                                 charges=_charges(structure, include),
                                 radicals=_radicals(structure, include),
                                 hydrogens=hydrogens, title=name))


def _species_note(structure) -> str:
    parts = []
    if structure.net_charge:
        parts.append(f"charge {structure.net_charge:+d}")
    if structure.net_radicals:
        parts.append(f"{structure.net_radicals} unpaired electron"
                     + ("s" if structure.net_radicals > 1 else ""))
    return ", ".join(parts) + " on the species"


def _flatten(structure, coords2d, include):
    """Screen coordinates -> molfile 2-D coordinates for the drawn atoms."""
    include = list(range(structure.natoms)) if include is None else list(include)
    position = {a: k for k, a in enumerate(include)}
    coords2d = np.asarray(coords2d, dtype=float).reshape(-1, 2)[include]

    bonds, orders = [], []
    for (i, j), order in zip(structure.bonds, structure.orders):
        if int(i) in position and int(j) in position:
            bonds.append((position[int(i)], position[int(j)]))
            orders.append(int(order))
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)

    scale = _scale(coords2d, bonds)
    xy = (coords2d - coords2d.mean(axis=0)) * scale
    coords = np.column_stack([xy[:, 0], -xy[:, 1], np.zeros(len(xy))])  # y is up here

    symbols = [structure.molecule.symbols[a] for a in include]
    hydrogens = [sum(1 for j in structure.neighbors[a]
                     if structure.molecule.symbols[j] == "H" and j not in position)
                 for a in include]
    return symbols, coords, bonds, np.array(orders, dtype=int), hydrogens


def _scale(coords2d, bonds) -> float:
    if not len(bonds):
        return 1.0
    lengths = np.linalg.norm(coords2d[bonds[:, 0]] - coords2d[bonds[:, 1]], axis=1)
    typical = float(np.median(lengths))
    return TARGET_BOND_LENGTH / typical if typical > 1e-9 else 1.0


def _charges(structure, include) -> list[int]:
    include = range(structure.natoms) if include is None else include
    return [int(structure.charges[a]) for a in include]


def _radicals(structure, include) -> list[int]:
    include = range(structure.natoms) if include is None else include
    return [int(structure.radicals[a]) for a in include]


def _block(symbols, coords, bonds, orders, charges=None, radicals=None,
           hydrogens=None, title: str = "") -> str:
    lines = [title[:80], "  microscope", ""]
    lines.append(f"{len(symbols):3d}{len(bonds):3d}"
                 "  0  0  0  0  0  0  0  0999 V2000")
    for k, symbol in enumerate(symbols):
        x, y, z = coords[k]
        charge = _CHARGE_CODES.get(charges[k], 0) if charges else 0
        # the valence field pins the implicit hydrogen count (15 means none)
        valence = 0
        if hydrogens is not None:
            bonded = int(sum(int(o) for b, o in zip(bonds, orders)
                             if k in (int(b[0]), int(b[1]))))
            valence = min(14, bonded + hydrogens[k]) or 15
        # dd ccc sss hhh bbb vvv then six unused fields, as V2000 lays them out
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3s}"
                     f" 0{charge:3d}  0  0  0{valence:3d}" + "  0" * 6)
    for (i, j), order in zip(bonds, orders):
        lines.append(f"{i + 1:3d}{j + 1:3d}{int(order):3d}  0  0  0  0")
    lines.extend(_property_lines("CHG", charges))
    lines.extend(_property_lines("RAD", radicals, code=lambda n: 2 if n == 1 else 3))
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def _property_lines(tag: str, values, code=lambda n: n) -> list[str]:
    """M  CHG / M  RAD lines, eight entries to a line as the format requires."""
    if not values:
        return []
    marked = [(k, code(v)) for k, v in enumerate(values) if v]
    lines = []
    for start in range(0, len(marked), 8):
        chunk = marked[start:start + 8]
        entries = "".join(f"{k + 1:4d}{v:4d}" for k, v in chunk)
        lines.append(f"M  {tag}{len(chunk):3d}{entries}")
    return lines
