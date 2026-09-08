"""Molecules the Lewis tests are built from.

Idealized geometries rather than files: the point of each is one bit of
chemistry — a ring that must alternate, a boron that must not become an anion
— and building them here keeps that visible in the test that uses it.
"""

import numpy as np

from microscope.core import lewis
from microscope.core.molecule import Molecule
from microscope.render import lewis2d
from microscope.render.camera import OrthoCamera

W, H = 1000.0, 800.0

def _molecule(symbols, coords, charge=0, multiplicity=1):
    mol = Molecule(symbols, np.array(coords, dtype=float), charge=charge,
                   multiplicity=multiplicity)
    mol.perceive_bonds()
    return mol


def _benzene():
    angles = np.arange(6) * np.pi / 3.0
    radius = 1.39                       # ring bond 1.39 A gives this radius too
    ring = np.column_stack([radius * np.cos(angles), radius * np.sin(angles),
                            np.zeros(6)])
    return _molecule(["C"] * 6 + ["H"] * 6,
                     np.vstack([ring, ring * (1.0 + 1.09 / radius)]))


def _acene(rings):
    """A linear polyacene: *rings* hexagons fused in a row, hydrogens included."""
    a = 1.40
    step = a * np.sqrt(3.0)
    carbons = [(x * step, y * a, 0.0)
               for x in np.arange(rings + 1) - 0.5 for y in (0.5, -0.5)]
    carbons += [(k * step, y * a, 0.0)
                for k in range(rings) for y in (1.0, -1.0)]
    carbons = np.array(carbons)

    hydrogens = []
    for c in carbons:                       # one H on every two-coordinate C
        deltas = carbons - c
        near = [d for d in deltas if 1e-6 < np.linalg.norm(d) < 1.6 * a]
        if len(near) != 2:
            continue
        out = -sum(d / np.linalg.norm(d) for d in near)
        hydrogens.append(c + 1.09 * out / np.linalg.norm(out))

    return _molecule(["C"] * len(carbons) + ["H"] * len(hydrogens),
                     np.vstack([carbons, hydrogens]))


def _water():
    return _molecule(["O", "H", "H"],
                     [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])


def _formaldehyde():
    return _molecule(["C", "O", "H", "H"],
                     [[0.0, 0.0, 0.0], [1.21, 0.0, 0.0],
                      [-0.57, 0.94, 0.0], [-0.57, -0.94, 0.0]])


def _nitrate():
    angles = np.arange(3) * 2.0 * np.pi / 3.0
    oxygens = np.column_stack([1.26 * np.cos(angles), 1.26 * np.sin(angles),
                               np.zeros(3)])
    return _molecule(["N", "O", "O", "O"], np.vstack([[0, 0, 0], oxygens]),
                     charge=-1)


def _camera(mol, half_height=None):
    camera = OrthoCamera()
    camera.fit(*mol.bounding_sphere())
    if half_height:
        camera.half_height = half_height
    return camera


def _orders(structure):
    counts = {1: 0, 2: 0, 3: 0}
    for order in structure.orders:
        counts[int(order)] += 1
    return counts


def _drawing(mol, options=None):
    """A structure plus its screen positions, as the writers take them."""
    options = options or lewis2d.LewisOptions()
    structure = lewis.perceive(mol)
    hidden = lewis2d.hidden_hydrogens(structure, options)
    include = [a for a in range(mol.natoms) if not hidden[a]]
    points = _camera(mol).project(mol.coords, W, H)
    return structure, points, include

def _dicyanobenzene():
    """A benzene ring with two para nitriles, all at reference lengths."""
    angles = np.arange(6) * np.pi / 3.0
    ring = np.column_stack([1.39 * np.cos(angles), 1.39 * np.sin(angles),
                            np.zeros(6)])
    out = ring / 1.39                                   # radial unit vectors
    symbols = ["C"] * 6
    coords = [ring]
    for k in range(6):
        if k in (0, 3):                                 # para positions: -C#N
            symbols += ["C", "N"]
            coords.append(np.array([ring[k] + 1.43 * out[k],
                                    ring[k] + (1.43 + 1.16) * out[k]]))
        else:
            symbols.append("H")
            coords.append(np.array([ring[k] + 1.09 * out[k]]))
    return _molecule(symbols, np.vstack(coords))
