"""Structure editing tests."""

import math

import numpy as np
import pytest

from microscope.core import editing, geometry
from microscope.core.molecule import Molecule


def ethane() -> Molecule:
    """Staggered ethane: C0-C1 along x, tetrahedral hydrogens."""
    c, s = -1.0 / 3.0, math.sqrt(8.0) / 3.0
    coords = [[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]]
    symbols = ["C", "C"]
    for phi_deg in (0, 120, 240):                 # H on C0
        phi = math.radians(phi_deg)
        coords.append([1.09 * c, 1.09 * s * math.cos(phi), 1.09 * s * math.sin(phi)])
        symbols.append("H")
    for phi_deg in (60, 180, 300):                # H on C1 (staggered)
        phi = math.radians(phi_deg)
        coords.append([1.54 - 1.09 * c, 1.09 * s * math.cos(phi), 1.09 * s * math.sin(phi)])
        symbols.append("H")
    mol = Molecule(symbols, np.array(coords))
    mol.perceive_bonds()
    return mol


def water() -> Molecule:
    mol = Molecule(["O", "H", "H"],
                   np.array([[0.0, 0.0, 0.117],
                             [0.0, 0.757, -0.469],
                             [0.0, -0.757, -0.469]]))
    mol.perceive_bonds()
    return mol


def test_ethane_topology():
    mol = ethane()
    assert len(mol.bonds) == 7
    frag = editing.fragment_beyond(mol.bonds, 0, 1)
    assert frag == [1, 5, 6, 7]                    # C1 and its hydrogens


def test_set_distance_moves_fragment():
    mol = ethane()
    old_ch = geometry.distance(mol.coords[1], mol.coords[5])
    coords = editing.set_distance(mol, 0, 1, 2.0)
    assert abs(np.linalg.norm(coords[1] - coords[0]) - 2.0) < 1e-9
    # C1's hydrogens moved along, keeping their local geometry
    assert abs(np.linalg.norm(coords[5] - coords[1]) - old_ch) < 1e-9
    # C0 side untouched
    np.testing.assert_allclose(coords[2], mol.coords[2])


def test_set_angle():
    mol = water()
    coords = editing.set_angle(mol, 1, 0, 2, 90.0)
    assert abs(geometry.angle(coords[1], coords[0], coords[2]) - 90.0) < 1e-6
    # O-H bond lengths preserved
    for h in (1, 2):
        assert abs(np.linalg.norm(coords[h] - coords[0])
                   - np.linalg.norm(mol.coords[h] - mol.coords[0])) < 1e-9


def test_set_dihedral():
    mol = ethane()
    d0 = geometry.dihedral(mol.coords[2], mol.coords[0], mol.coords[1], mol.coords[5])
    assert abs(abs(d0) - 60.0) < 1.0               # staggered
    coords = editing.set_dihedral(mol, 2, 0, 1, 5, 0.0)   # eclipse it
    d1 = geometry.dihedral(coords[2], coords[0], coords[1], coords[5])
    assert abs(d1) < 1e-6
    # the other C1 hydrogens rotated by the same amount (still tetrahedral)
    ang = geometry.angle(coords[5], coords[1], coords[6])
    assert abs(ang - geometry.angle(mol.coords[5], mol.coords[1], mol.coords[6])) < 1e-6
    # C0 hydrogens untouched
    np.testing.assert_allclose(coords[2], mol.coords[2])


def test_ring_bond_falls_back_to_single_atom():
    # benzene ring: adjusting a ring bond moves only the end atom
    rc, rh = 1.397, 2.481
    symbols, coords = [], []
    for k in range(6):
        a = math.radians(60 * k)
        symbols += ["C", "H"]
        coords += [[rc * math.cos(a), rc * math.sin(a), 0.0],
                   [rh * math.cos(a), rh * math.sin(a), 0.0]]
    mol = Molecule(symbols, np.array(coords))
    mol.perceive_bonds()
    assert editing.fragment_beyond(mol.bonds, 0, 2) is None
    coords2 = editing.set_distance(mol, 0, 2, 1.6)
    moved = np.where(np.abs(coords2 - mol.coords).max(axis=1) > 1e-9)[0]
    assert list(moved) == [2]


def test_delete_atoms():
    mol = water()
    new = editing.delete_atoms(mol, [2])
    assert new.formula() == "HO"
    assert new.natoms == 2
    assert len(new.bonds) == 1
    with pytest.raises(editing.EditError):
        editing.delete_atoms(mol, [0, 1, 2])
