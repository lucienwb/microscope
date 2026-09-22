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
    # looked up once per element, not once per atom: a protein is five elements
    # and a hundred thousand atoms
    elements, which = np.unique(zs, return_inverse=True)
    radii = np.array([style.atom_radius(z) for z in elements], dtype=np.float32)[which]
    colors = np.array([style.atom_color(z) for z in elements],
                      dtype=np.float32).reshape(-1, 3)[which]

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
    """Short cylinder segments forming dashed H···acceptor lines, every dash of
    every H-bond at once."""
    found = find_hbonds(molecule)
    if not found:
        return np.zeros((0, 13), dtype=np.float32)
    pairs = np.array([(h, a) for h, a, _ in found], dtype=np.int64)
    p1 = molecule.coords[pairs[:, 0]].astype(np.float32)
    p2 = molecule.coords[pairs[:, 1]].astype(np.float32)
    length = np.linalg.norm(p2 - p1, axis=1)
    usable = length >= 1e-6
    p1, p2, length = p1[usable], p2[usable], length[usable]
    u = (p2 - p1) / length[:, None]
    # dashes start 0.24 A out and stop 0.26 A short, clear of the atom spheres
    period = style.hbond_dash + style.hbond_gap
    end = length - 0.26
    count = np.maximum(np.ceil((end - 0.24) / period), 0).astype(np.int64)
    which = np.repeat(np.arange(len(length)), count)
    k = np.arange(len(which)) - np.repeat(np.cumsum(count) - count, count)
    start = 0.24 + k * period
    stop = np.minimum(start + style.hbond_dash, end[which])
    ends = np.column_stack([p1[which] + u[which] * start[:, None],
                            p1[which] + u[which] * stop[:, None]])
    color = np.asarray(style.hbond_color, dtype=np.float32)
    return np.hstack([ends, np.full((len(which), 1), style.hbond_radius),
                      np.tile(color, (len(which), 2))]).astype(np.float32)
