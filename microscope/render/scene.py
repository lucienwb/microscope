"""Build GPU instance buffers from a Molecule + Style."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.contacts import find_hbonds
from ..core.isosurface import isosurface_mesh
from ..core.molecule import Molecule
from ..core.volume import VolumeData
from .styles import Style


# per-atom representations
REP_BALL = 0     # ball-and-stick (the default "CPK" look)
REP_STICK = 1    # bonds only, atoms shrunk to smooth tube joints
REP_LINE = 2     # thin wireframe
LINE_RADIUS = 0.032          # Angstrom


@dataclass
class SceneBuffers:
    spheres: np.ndarray     # (N, 8)  float32: center xyz, radius, color rgb, quadrant flag
    cylinders: np.ndarray   # (M, 13) float32: A xyz, B xyz, radius, colorA rgb, colorB rgb


def build_scene(molecule: Molecule, style: Style,
                reps: np.ndarray | None = None) -> SceneBuffers:
    """*reps* holds one REP_* code per atom (None = all REP_BALL); a bond is
    drawn at the thinnest representation of its two atoms."""
    zs = np.asarray(molecule.atomic_numbers)
    coords = molecule.coords.astype(np.float32)
    radii = np.array([style.atom_radius(z) for z in zs], dtype=np.float32)
    colors = np.array([style.atom_color(z) for z in zs], dtype=np.float32)

    if reps is None or len(reps) != len(zs):
        reps = np.zeros(len(zs), dtype=int)
    reps = np.asarray(reps, dtype=int)
    radii[reps == REP_STICK] = style.bond_radius
    radii[reps == REP_LINE] = LINE_RADIUS
    if style.bond_color is not None:
        # uniform-bond styles (Houkmol): tube joints take the bond color
        colors[reps != REP_BALL] = np.asarray(style.bond_color, dtype=np.float32)

    quad = np.zeros((len(zs), 1), dtype=np.float32)
    if style.quadrant_color is not None:
        quad[(zs > 1) & (reps == REP_BALL)] = 1.0   # seam lines on heavy atoms
    spheres = np.hstack([coords, radii[:, None], colors, quad]).astype(np.float32)

    bonds = molecule.bonds
    if bonds is None or len(bonds) == 0:
        cylinders = np.zeros((0, 13), dtype=np.float32)
    else:
        i, j = bonds[:, 0], bonds[:, 1]
        r = np.full((len(bonds), 1), style.bond_radius, dtype=np.float32)
        r[(reps[i] == REP_LINE) | (reps[j] == REP_LINE)] = LINE_RADIUS
        if style.bond_color is not None:
            ca = np.tile(np.asarray(style.bond_color, dtype=np.float32), (len(bonds), 1))
            cb = ca
        else:
            ca, cb = colors[i], colors[j]
        cylinders = np.hstack([coords[i], coords[j], r, ca, cb]).astype(np.float32)

    if style.show_hbonds:
        dashes = _hbond_dashes(molecule, style)
        if len(dashes):
            cylinders = np.vstack([cylinders, dashes])

    return SceneBuffers(spheres=np.ascontiguousarray(spheres),
                        cylinders=np.ascontiguousarray(cylinders))


def build_surface_meshes(volume: VolumeData, isovalue: float, style: Style) -> list:
    """Isosurface meshes for the renderer: (vertices, normals, rgba) tuples.

    Orbital-like signed data gets both the +isovalue and -isovalue lobes in
    the style's two surface colors; densities get the positive surface only.
    """
    level = abs(float(isovalue))
    if level == 0.0:
        return []
    meshes = []
    verts, normals = isosurface_mesh(volume, level)
    if len(verts):
        meshes.append((verts, normals, (*style.surface_positive, style.surface_opacity)))
    if volume.is_signed:
        verts, normals = isosurface_mesh(volume, -level)
        if len(verts):
            meshes.append((verts, normals,
                           (*style.surface_negative, style.surface_opacity)))
    return meshes


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
