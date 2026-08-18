"""Structure editing: distance/angle/dihedral adjustment, atom deletion.

All operations are fragment-aware: adjusting an internal coordinate moves the
whole group of atoms attached beyond the adjusted bond (GaussView behavior).
If the coordinate lies in a ring, only the single end atom is moved instead.
Functions return a new coordinates array (or Molecule) and never modify input.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .geometry import angle as _angle
from .geometry import dihedral as _dihedral
from .geometry import rotation_matrix
from .molecule import Molecule


class EditError(ValueError):
    pass


def _adjacency(bonds) -> dict[int, set[int]]:
    adj: dict[int, set[int]] = defaultdict(set)
    if bonds is not None:
        for a, b in bonds:
            adj[int(a)].add(int(b))
            adj[int(b)].add(int(a))
    return adj


def fragment_beyond(bonds, anchor: int, moving: int) -> list[int] | None:
    """Atoms reachable from *moving* without crossing the anchor-moving bond.

    Returns None when *anchor* is reachable anyway, i.e. the bond is part of
    a ring and there is no clean fragment to move.
    """
    adj = _adjacency(bonds)
    seen = {moving}
    stack = [moving]
    while stack:
        a = stack.pop()
        for b in adj[a]:
            if {a, b} == {anchor, moving}:
                continue
            if b not in seen:
                seen.add(b)
                stack.append(b)
    if anchor in seen:
        return None
    return sorted(seen)


def _moving_set(mol: Molecule, anchor: int, moving: int,
                move_fragment: bool) -> list[int]:
    if move_fragment:
        frag = fragment_beyond(mol.bonds, anchor, moving)
        if frag is not None:
            return frag
    return [moving]


def set_distance(mol: Molecule, i: int, j: int, value: float,
                 move_fragment: bool = True) -> np.ndarray:
    """Coordinates with d(i, j) set to *value*, translating j's fragment."""
    if value <= 0:
        raise EditError("distance must be positive")
    coords = mol.coords.copy()
    d = coords[j] - coords[i]
    norm = float(np.linalg.norm(d))
    if norm < 1e-6:
        raise EditError("atoms coincide; cannot define a direction")
    frag = _moving_set(mol, i, j, move_fragment)
    coords[frag] += (value - norm) * d / norm
    return coords


def set_angle(mol: Molecule, i: int, j: int, k: int, value: float,
              move_fragment: bool = True) -> np.ndarray:
    """Coordinates with angle i-j-k set to *value* degrees, rotating k's side."""
    if not 0.0 < value < 180.0:
        raise EditError("angle must be between 0 and 180 degrees")
    coords = mol.coords.copy()
    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]
    axis = np.cross(v1, v2)
    if np.linalg.norm(axis) < 1e-8:  # collinear: any perpendicular axis works
        helper = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(v1, helper)) > 0.9 * np.linalg.norm(v1):
            helper = np.array([0.0, 1.0, 0.0])
        axis = np.cross(v1, helper)
    frag = _moving_set(mol, j, k, move_fragment)
    pivot = coords[j]

    def rotated(delta_deg: float) -> np.ndarray:
        rot = rotation_matrix(axis, np.radians(delta_deg))
        out = coords.copy()
        out[frag] = (rot @ (out[frag] - pivot).T).T + pivot
        return out

    current = _angle(coords[i], coords[j], coords[k])
    result = rotated(value - current)
    if abs(_angle(result[i], result[j], result[k]) - value) > 1e-4:
        result = rotated(current - value)
    return result


def _wrap_deg(a: float) -> float:
    return (a + 180.0) % 360.0 - 180.0


def set_dihedral(mol: Molecule, i: int, j: int, k: int, l: int, value: float,
                 move_fragment: bool = True) -> np.ndarray:
    """Coordinates with dihedral i-j-k-l set to *value* degrees.

    Rotates everything attached through k about the j-k axis; if that bond is
    in a ring, only atom l is rotated.
    """
    coords = mol.coords.copy()
    axis = coords[k] - coords[j]
    if np.linalg.norm(axis) < 1e-6:
        raise EditError("axis atoms coincide")
    if move_fragment:
        frag = fragment_beyond(mol.bonds, j, k)
        if frag is None:
            frag = [l]
    else:
        frag = [l]
    pivot = coords[j]

    def rotated(delta_deg: float) -> np.ndarray:
        rot = rotation_matrix(axis, np.radians(delta_deg))
        out = coords.copy()
        out[frag] = (rot @ (out[frag] - pivot).T).T + pivot
        return out

    current = _dihedral(coords[i], coords[j], coords[k], coords[l])
    result = rotated(value - current)
    achieved = _dihedral(result[i], result[j], result[k], result[l])
    if abs(_wrap_deg(achieved - value)) > 1e-3:
        result = rotated(current - value)
    return result


def translate_atoms(mol: Molecule, indices, delta) -> np.ndarray:
    """Coordinates with *indices* shifted rigidly by the 3-vector *delta*."""
    idx = _unique(indices)
    if not idx:
        raise EditError("no atoms to move")
    coords = mol.coords.copy()
    coords[idx] += np.asarray(delta, dtype=float)
    return coords


def rotate_atoms(mol: Molecule, indices, axis, angle_deg: float,
                 pivot=None) -> np.ndarray:
    """Coordinates with *indices* rotated about *axis* through *pivot*.

    *pivot* defaults to the centroid of the moved atoms, which is what the
    interactive manipulator wants (the selection turns in place).
    """
    idx = _unique(indices)
    if not idx:
        raise EditError("no atoms to rotate")
    axis = np.asarray(axis, dtype=float)
    if np.linalg.norm(axis) < 1e-9:
        raise EditError("degenerate rotation axis")
    coords = mol.coords.copy()
    pivot = (coords[idx].mean(axis=0) if pivot is None
             else np.asarray(pivot, dtype=float))
    rot = rotation_matrix(axis, np.radians(angle_deg))
    coords[idx] = (rot @ (coords[idx] - pivot).T).T + pivot
    return coords


def connected_fragment(bonds, seeds) -> list[int]:
    """Every atom bonded, directly or indirectly, to any atom in *seeds*."""
    adj = _adjacency(bonds)
    seen: set[int] = set()
    stack = [int(s) for s in seeds]
    while stack:
        a = stack.pop()
        if a in seen:
            continue
        seen.add(a)
        stack.extend(b for b in adj[a] if b not in seen)
    return sorted(seen)


def _unique(indices) -> list[int]:
    """Indices without duplicates, order irrelevant (used as a fancy index)."""
    return sorted({int(i) for i in indices})


def delete_atoms(mol: Molecule, indices) -> Molecule:
    """New Molecule with *indices* removed; bonds are re-perceived."""
    drop = set(int(x) for x in indices)
    keep = [x for x in range(mol.natoms) if x not in drop]
    if not keep:
        raise EditError("cannot delete every atom")
    new = Molecule(
        symbols=[mol.symbols[x] for x in keep],
        coords=mol.coords[keep],
        charge=mol.charge,
        multiplicity=mol.multiplicity,
        title=mol.title,
    )
    new.perceive_bonds()
    return new
