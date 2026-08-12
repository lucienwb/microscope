"""Build GPU instance buffers from a Molecule + Style."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.molecule import Molecule
from .styles import Style


@dataclass
class SceneBuffers:
    spheres: np.ndarray     # (N, 7)  float32: center xyz, radius, color rgb
    cylinders: np.ndarray   # (M, 13) float32: A xyz, B xyz, radius, colorA rgb, colorB rgb


def build_scene(molecule: Molecule, style: Style) -> SceneBuffers:
    zs = molecule.atomic_numbers
    coords = molecule.coords.astype(np.float32)
    radii = np.array([style.atom_radius(z) for z in zs], dtype=np.float32)
    colors = np.array([style.atom_color(z) for z in zs], dtype=np.float32)

    spheres = np.hstack([coords, radii[:, None], colors]).astype(np.float32)

    bonds = molecule.bonds
    if bonds is None or len(bonds) == 0:
        cylinders = np.zeros((0, 13), dtype=np.float32)
    else:
        i, j = bonds[:, 0], bonds[:, 1]
        r = np.full((len(bonds), 1), style.bond_radius, dtype=np.float32)
        cylinders = np.hstack([coords[i], coords[j], r, colors[i], colors[j]]).astype(np.float32)

    return SceneBuffers(spheres=np.ascontiguousarray(spheres),
                        cylinders=np.ascontiguousarray(cylinders))
