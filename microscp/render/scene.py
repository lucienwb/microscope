"""Build GPU instance buffers from a Molecule + Style."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.contacts import find_hbonds
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

    if style.show_hbonds:
        dashes = _hbond_dashes(molecule, style)
        if len(dashes):
            cylinders = np.vstack([cylinders, dashes])

    return SceneBuffers(spheres=np.ascontiguousarray(spheres),
                        cylinders=np.ascontiguousarray(cylinders))


def _hbond_dashes(molecule: Molecule, style: Style) -> np.ndarray:
    """Short cylinder segments forming dashed H···acceptor lines."""
    rows = []
    color = np.asarray(style.hbond_color, dtype=np.float32)
    for h, acc, _dist in find_hbonds(molecule):
        p1 = molecule.coords[h].astype(np.float32)
        p2 = molecule.coords[acc].astype(np.float32)
        length = float(np.linalg.norm(p2 - p1))
        if length < 1e-6:
            continue
        u = (p2 - p1) / length
        # keep dashes clear of the atom spheres at both ends
        s = 0.24
        end = length - 0.26
        while s < end:
            e = min(s + style.hbond_dash, end)
            rows.append(np.concatenate([
                p1 + u * s, p1 + u * e, [style.hbond_radius], color, color]))
            s = e + style.hbond_gap
    if not rows:
        return np.zeros((0, 13), dtype=np.float32)
    return np.array(rows, dtype=np.float32)
