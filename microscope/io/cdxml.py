"""ChemDraw CDXML writer.

CDXML is ChemDraw's documented XML format, so a structure written here opens
natively in ChemDraw with the double bonds, charges and implicit hydrogens
already in place — ready to be rearranged, labelled or dropped into a scheme.

The coordinates written are the ones on screen: CDXML measures in points with
y running down the page, which is the same convention as the viewport, so the
drawing arrives in ChemDraw at the angle it was composed at. Only the scale
changes, to put the median bond at ChemDraw's standard 30-point length.
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import quoteattr

import numpy as np

from .. import __version__

BOND_LENGTH = 30.0            # ChemDraw's default, in points
MARGIN = 40.0                 # keeps the drawing off the corner of the page

_HEADER = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">\n'
)


def write(path, structure, coords2d, include=None) -> None:
    """Write the drawn part of *structure* at its screen positions *coords2d*.

    *include* lists the atoms actually drawn; hydrogens left out of it are
    written as implicit hydrogen counts on the atom they belong to, which is
    how ChemDraw itself stores a skeletal structure.

    Only what ChemDraw is certain to understand goes in the file — element,
    position, charge, implicit hydrogens and bond order. Unpaired electrons
    are not written; they survive in the MDL molfile instead.
    """
    include = list(range(structure.natoms)) if include is None else list(include)
    drawn = set(include)
    points = _page_coordinates(structure, coords2d, include)
    symbols = structure.molecule.symbols
    numbers = structure.molecule.atomic_numbers

    body = []
    for node_id, atom in enumerate(include, start=1):
        x, y = points[node_id - 1]
        attributes = [f'id="{node_id}"', f'p="{x:.2f} {y:.2f}"',
                      f'Element="{int(numbers[atom])}"']
        implicit = sum(1 for j in structure.neighbors[atom]
                       if symbols[j] == "H" and j not in drawn)
        if implicit:
            attributes.append(f'NumHydrogens="{implicit}"')
        charge = int(structure.charges[atom])
        if charge:
            attributes.append(f'Charge="{charge}"')
        body.append("     <n " + " ".join(attributes) + "/>")

    node_of = {atom: node_id for node_id, atom in enumerate(include, start=1)}
    bond_id = len(include)
    for (i, j), order in zip(structure.bonds, structure.orders):
        i, j = int(i), int(j)
        if i not in drawn or j not in drawn:
            continue
        bond_id += 1
        body.append(f'     <b id="{bond_id}" B="{node_of[i]}" E="{node_of[j]}"'
                    f' Order="{int(order)}"/>')

    text = (_HEADER
            + f'<CDXML CreationProgram={quoteattr("microscope " + __version__)}'
              f' BondLength="{BOND_LENGTH:.0f}">\n'
              ' <page HeightPages="1" WidthPages="1">\n'
              '  <fragment>\n'
            + "\n".join(body)
            + "\n  </fragment>\n </page>\n</CDXML>\n")
    Path(path).write_text(text)


def _page_coordinates(structure, coords2d, include) -> np.ndarray:
    """Screen pixels -> CDXML points, with the standard bond length."""
    points = np.asarray(coords2d, dtype=float).reshape(-1, 2)[include]
    position = {a: k for k, a in enumerate(include)}
    lengths = [float(np.linalg.norm(points[position[int(i)]] - points[position[int(j)]]))
               for i, j in structure.bonds
               if int(i) in position and int(j) in position]
    typical = float(np.median(lengths)) if lengths else 0.0
    scale = BOND_LENGTH / typical if typical > 1e-9 else 1.0
    points = points * scale
    return points - points.min(axis=0) + MARGIN
