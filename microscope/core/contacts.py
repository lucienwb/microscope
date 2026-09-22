"""Non-covalent contact detection (hydrogen bonds)."""

from __future__ import annotations

import numpy as np

from .molecule import Molecule
from .neighbors import close_pairs

# Elements that can act as hydrogen-bond donors/acceptors.
HBOND_ELEMENTS = {7, 8, 9, 16, 17}


def find_hbonds(molecule: Molecule, max_h_acceptor: float = 2.6,
                min_angle: float = 120.0) -> list[tuple[int, int, float]]:
    """Detect D–H···A hydrogen bonds.

    Criteria: H covalently bonded to a donor (N/O/F/S/Cl), an acceptor atom of
    the same element set not bonded to H, d(H···A) <= *max_h_acceptor* Å and
    angle(D–H···A) >= *min_angle* degrees.

    Returns a list of (h_index, acceptor_index, distance), by hydrogen then
    acceptor. The candidates come from one neighbour search over the
    hydrogens and polar atoms together, not from every hydrogen against every
    acceptor: on a protein with its hydrogens that was minutes.
    """
    bonds = molecule.bonds
    if bonds is None:
        bonds = molecule.perceive_bonds()
    zs = molecule.atomic_numbers
    coords = molecule.coords
    n = len(zs)
    polar = np.isin(zs, list(HBOND_ELEMENTS))
    bonds = np.asarray(bonds, dtype=np.int64).reshape(-1, 2)

    # each hydrogen's donor: the polar atom it is bonded to (the lowest-numbered,
    # should it somehow have two)
    ends = np.concatenate([bonds, bonds[:, ::-1]])
    ends = ends[(zs[ends[:, 0]] == 1) & polar[ends[:, 1]]]
    donor = np.full(n, -1, dtype=np.int64)
    ends = ends[np.lexsort((-ends[:, 1], ends[:, 0]))]     # lowest partner written last
    donor[ends[:, 0]] = ends[:, 1]
    hydrogen = donor >= 0
    if not hydrogen.any() or not polar.any():
        return []

    i, j, dist = close_pairs(coords, max_h_acceptor, among=np.flatnonzero(hydrogen | polar))
    swap = ~hydrogen[i]                          # put the hydrogen first
    h, acc = np.where(swap, j, i), np.where(swap, i, j)
    keep = hydrogen[h] & polar[acc] & (acc != donor[h]) & (dist >= 1.2)
    h, acc, dist = h[keep], acc[keep], dist[keep]
    # an acceptor bonded to the hydrogen itself is not an acceptor
    low, high = np.minimum(bonds[:, 0], bonds[:, 1]), np.maximum(bonds[:, 0], bonds[:, 1])
    bonded = np.isin(np.minimum(h, acc) * n + np.maximum(h, acc), low * n + high)
    h, acc, dist = h[~bonded], acc[~bonded], dist[~bonded]
    to_donor = coords[donor[h]] - coords[h]
    to_acceptor = coords[acc] - coords[h]
    cosine = np.einsum("ij,ij->i", to_donor, to_acceptor) / (
        np.linalg.norm(to_donor, axis=1) * np.linalg.norm(to_acceptor, axis=1))
    wide = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))) >= min_angle
    h, acc, dist = h[wide], acc[wide], dist[wide]
    order = np.lexsort((acc, h))
    return [(int(a), int(b), float(d)) for a, b, d in zip(h[order], acc[order], dist[order])]
