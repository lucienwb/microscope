"""Non-covalent contact detection (hydrogen bonds)."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from . import geometry
from .molecule import Molecule

# Elements that can act as hydrogen-bond donors/acceptors.
HBOND_ELEMENTS = {7, 8, 9, 16, 17}


def find_hbonds(molecule: Molecule, max_h_acceptor: float = 2.6,
                min_angle: float = 120.0) -> list[tuple[int, int, float]]:
    """Detect D–H···A hydrogen bonds.

    Criteria: H covalently bonded to a donor (N/O/F/S/Cl), an acceptor atom of
    the same element set not bonded to H, d(H···A) <= *max_h_acceptor* Å and
    angle(D–H···A) >= *min_angle* degrees.

    Returns a list of (h_index, acceptor_index, distance).
    """
    bonds = molecule.bonds
    if bonds is None:
        bonds = molecule.perceive_bonds()
    zs = molecule.atomic_numbers
    coords = molecule.coords
    adjacency: dict[int, set[int]] = defaultdict(set)
    for a, b in bonds:
        adjacency[int(a)].add(int(b))
        adjacency[int(b)].add(int(a))

    acceptors = np.where(np.isin(zs, list(HBOND_ELEMENTS)))[0]
    hbonds: list[tuple[int, int, float]] = []
    for h in np.where(zs == 1)[0]:
        donors = [d for d in adjacency[int(h)] if zs[d] in HBOND_ELEMENTS]
        if not donors:
            continue
        donor = donors[0]
        for acc in acceptors:
            acc = int(acc)
            if acc == donor or acc in adjacency[int(h)]:
                continue
            dist = float(np.linalg.norm(coords[acc] - coords[h]))
            if not 1.2 <= dist <= max_h_acceptor:
                continue
            angle = geometry.angle(coords[donor], coords[h], coords[acc])
            if angle >= min_angle:
                hbonds.append((int(h), acc, dist))
    return hbonds
