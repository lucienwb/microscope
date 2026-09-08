"""Isosurface extraction from volumetric data — pure numpy, no external mesher.

Uses marching tetrahedra (each grid cell split into 6 tetrahedra around the
0-6 diagonal): the per-tetrahedron case table is tiny and unambiguous, unlike
the 256-entry marching-cubes tables. Vertex normals come from the trilinearly
interpolated gradient of the scalar field, which gives smooth shading without
any mesh connectivity bookkeeping.
"""

from __future__ import annotations

import numpy as np

from .volume import VolumeData

_CORNERS = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                     (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], dtype=np.int64)

# Six tetrahedra per cell, all sharing the 0-6 body diagonal.
_TETS = ((0, 5, 1, 6), (0, 1, 2, 6), (0, 2, 3, 6),
         (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6))

# Triangles per tetrahedron sign case (bit i set = local corner i above the
# level). Each triangle is 3 edges (a, b) whose crossing points are joined;
# complementary cases share the same crossed edges.
_ONE = ((((0, 1), (0, 2), (0, 3)),),          # corner 0 above
        (((1, 0), (1, 3), (1, 2)),),          # corner 1
        (((2, 0), (2, 1), (2, 3)),),          # corner 2
        (((3, 0), (3, 2), (3, 1)),))          # corner 3
_TWO = {
    0b0011: (((0, 2), (0, 3), (1, 3)), ((0, 2), (1, 3), (1, 2))),
    0b0101: (((0, 1), (0, 3), (2, 3)), ((0, 1), (2, 3), (2, 1))),
    0b1001: (((0, 1), (0, 2), (2, 3)), ((0, 1), (2, 3), (1, 3))),
}
_CASE_TRIS = {}
for _i in range(4):
    _CASE_TRIS[1 << _i] = _ONE[_i]
    _CASE_TRIS[0b1111 ^ (1 << _i)] = _ONE[_i]
for _code, _tris in _TWO.items():
    _CASE_TRIS[_code] = _tris
    _CASE_TRIS[0b1111 ^ _code] = _tris


def marching_tetrahedra(values: np.ndarray, level: float) -> np.ndarray:
    """Extract the *level* isosurface of a 3D grid as a triangle soup.

    Returns an array of shape (ntriangles, 3, 3): triangle corners in
    fractional grid-index coordinates.
    """
    values = np.asarray(values, dtype=np.float32)
    if values.ndim != 3 or min(values.shape) < 2:
        return np.zeros((0, 3, 3), dtype=np.float32)
    nx, ny, nz = values.shape
    corner_vals = [values[dx:nx - 1 + dx, dy:ny - 1 + dy, dz:nz - 1 + dz]
                   for dx, dy, dz in _CORNERS]
    straddling = ((np.minimum.reduce(corner_vals) <= level)
                  & (level < np.maximum.reduce(corner_vals)))
    # Only the shell of cells the surface passes through is kept, and their
    # indices fit an int32 many times over; the corner positions are built one
    # tetrahedron at a time rather than held as an (M, 8, 3) array, which on a
    # 128^3 grid is the difference between 150 MB and 40.
    cells = np.argwhere(straddling).astype(np.int32)
    if len(cells) == 0:
        return np.zeros((0, 3, 3), dtype=np.float32)

    cell_vals = np.stack([values[cells[:, 0] + dx, cells[:, 1] + dy, cells[:, 2] + dz]
                          for dx, dy, dz in _CORNERS], axis=1)          # (M, 8)
    base = cells.astype(np.float32)

    triangles = []
    for tet in _TETS:
        v = cell_vals[:, tet]                                           # (M, 4)
        corners = _CORNERS[list(tet)].astype(np.float32)                # (4, 3)
        p = base[:, None, :] + corners[None, :, :]                      # (M, 4, 3)
        case = ((v[:, 0] > level) * 1 + (v[:, 1] > level) * 2
                + (v[:, 2] > level) * 4 + (v[:, 3] > level) * 8)
        for code, tris in _CASE_TRIS.items():
            sel = case == code
            if not sel.any():
                continue
            vs, ps = v[sel], p[sel]
            for edges in tris:
                pts = []
                for a, b in edges:
                    t = (level - vs[:, a]) / (vs[:, b] - vs[:, a])
                    pts.append(ps[:, a] + t[:, None] * (ps[:, b] - ps[:, a]))
                triangles.append(np.stack(pts, axis=1))
    if not triangles:
        return np.zeros((0, 3, 3), dtype=np.float32)
    return np.concatenate(triangles, axis=0)


def _trilinear(grid: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Sample *grid* at fractional index positions *pts* (N, 3).

    The corner indices are int32: a grid axis longer than two billion points
    is not a thing, and at half a million sample points the index arrays are
    otherwise the largest temporaries in the whole isosurface.
    """
    n = np.array(grid.shape)
    p = np.clip(pts, 0.0, n - 1.0)
    i0 = np.minimum(p.astype(np.int32), (n - 2).astype(np.int32))
    f = p - i0
    x0, y0, z0 = i0[:, 0], i0[:, 1], i0[:, 2]
    fx, fy, fz = f[:, 0], f[:, 1], f[:, 2]
    out = np.zeros(len(pts), dtype=np.float32)
    for dx in (0, 1):
        wx = fx if dx else 1.0 - fx
        for dy in (0, 1):
            wy = fy if dy else 1.0 - fy
            for dz in (0, 1):
                wz = fz if dz else 1.0 - fz
                out += wx * wy * wz * grid[x0 + dx, y0 + dy, z0 + dz]
    return out


def isosurface_mesh(volume: VolumeData, level: float) -> tuple[np.ndarray, np.ndarray]:
    """Extract the *level* isosurface of *volume* in world coordinates.

    Returns ``(vertices, normals)`` float32 arrays of shape (N, 3) where
    consecutive triples of vertices form triangles. Normals are unit vectors
    pointing away from the enclosed (above-level) region.
    """
    tris = marching_tetrahedra(volume.values, level)
    verts = tris.reshape(-1, 3)
    if len(verts) == 0:
        empty = np.zeros((0, 3), dtype=np.float32)
        return empty, empty.copy()
    # one gradient component at a time: all three together is three more
    # copies of the whole grid, and only one is needed at once
    field = np.asarray(volume.values, dtype=np.float32)
    g = np.stack([_trilinear(np.gradient(field, axis=k), verts)
                  for k in range(3)], axis=1)
    normals = -(g @ np.linalg.inv(volume.axes).T)   # index-space -> world
    norm = np.linalg.norm(normals, axis=1)
    bad = norm < 1e-12
    normals[bad] = (0.0, 0.0, 1.0)
    norm[bad] = 1.0
    normals /= norm[:, None]
    world = volume.grid_to_world(verts)
    return world.astype(np.float32), normals.astype(np.float32)
