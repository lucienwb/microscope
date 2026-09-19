"""Isosurface extraction from volumetric data — pure numpy, no external mesher.

Uses marching tetrahedra (each grid cell split into 6 tetrahedra around the
0-6 diagonal): the per-tetrahedron case table is tiny and unambiguous, unlike
the 256-entry marching-cubes tables. Vertex normals come from the gradient of
the scalar field, interpolated along the grid edge each vertex sits on, which
gives smooth shading without any mesh connectivity bookkeeping.
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
    return _march(values, level, normals=False)[0]


def _straddling(values: np.ndarray, level: float) -> np.ndarray:
    """Cells with corners on both sides of the level. The running min and max
    are kept in place: reducing a list of the eight corner views first stacks
    them, eight copies of the grid, which on an orbital grid is 250 MB."""
    nx, ny, nz = values.shape
    views = [values[dx:nx - 1 + dx, dy:ny - 1 + dy, dz:nz - 1 + dz] for dx, dy, dz in _CORNERS]
    low, high = views[0].copy(), views[0].copy()
    for view in views[1:]:
        np.minimum(low, view, out=low)
        np.maximum(high, view, out=high)
    return (low <= level) & (level < high)


def _corner_gradients(values: np.ndarray, corners: list) -> np.ndarray:
    """(M, 8, 3) central-difference gradient at each corner of the kept cells,
    one-sided at the edges of the grid, as np.gradient does - but only at the
    points the surface needs, not over the whole grid."""
    n = values.shape
    out = np.empty((len(corners[0][0]), 8, 3), dtype=np.float32)
    for c, index in enumerate(corners):
        for k in range(3):
            up, down = list(index), list(index)
            up[k] = np.minimum(index[k] + 1, n[k] - 1)
            down[k] = np.maximum(index[k] - 1, 0)
            out[:, c, k] = ((values[tuple(up)] - values[tuple(down)])
                            / (up[k] - down[k]).astype(np.float32))
    return out


def _march(values: np.ndarray, level: float, normals: bool):
    """Triangles, and with *normals* the field's gradient at each vertex:
    every vertex lies on an edge between two grid points, so its gradient is
    the one at those two points, interpolated as the vertex is."""
    values = np.asarray(values, dtype=np.float32)
    empty = np.zeros((0, 3, 3), dtype=np.float32)
    if values.ndim != 3 or min(values.shape) < 2:
        return empty, empty.copy()
    # Only the shell of cells the surface passes through is kept, and their
    # indices fit an int32 many times over; the corner positions are built one
    # tetrahedron at a time rather than held as an (M, 8, 3) array.
    cells = np.argwhere(_straddling(values, level)).astype(np.int32)
    if len(cells) == 0:
        return empty, empty.copy()
    corners = [(cells[:, 0] + dx, cells[:, 1] + dy, cells[:, 2] + dz)
               for dx, dy, dz in _CORNERS]
    cell_vals = np.stack([values[index] for index in corners], axis=1)      # (M, 8)
    cell_grads = _corner_gradients(values, corners) if normals else None
    base = cells.astype(np.float32)

    triangles, gradients = [], []
    for tet in _TETS:
        v = cell_vals[:, tet]                                           # (M, 4)
        offsets = _CORNERS[list(tet)].astype(np.float32)                # (4, 3)
        p = base[:, None, :] + offsets[None, :, :]                      # (M, 4, 3)
        g = cell_grads[:, tet] if normals else None                     # (M, 4, 3)
        case = ((v[:, 0] > level) * 1 + (v[:, 1] > level) * 2
                + (v[:, 2] > level) * 4 + (v[:, 3] > level) * 8)
        for code, tris in _CASE_TRIS.items():
            sel = case == code
            if not sel.any():
                continue
            vs, ps = v[sel], p[sel]
            gs = g[sel] if normals else None
            for edges in tris:
                pts, grads = [], []
                for a, b in edges:
                    t = ((level - vs[:, a]) / (vs[:, b] - vs[:, a]))[:, None]
                    pts.append(ps[:, a] + t * (ps[:, b] - ps[:, a]))
                    if normals:
                        grads.append(gs[:, a] + t * (gs[:, b] - gs[:, a]))
                triangles.append(np.stack(pts, axis=1))
                if normals:
                    gradients.append(np.stack(grads, axis=1))
    if not triangles:
        return empty, empty.copy()
    return (np.concatenate(triangles, axis=0),
            np.concatenate(gradients, axis=0) if normals else empty)


def isosurface_mesh(volume: VolumeData, level: float) -> tuple[np.ndarray, np.ndarray]:
    """Extract the *level* isosurface of *volume* in world coordinates.

    Returns ``(vertices, normals)`` float32 arrays of shape (N, 3) where
    consecutive triples of vertices form triangles. Normals are unit vectors
    pointing away from the enclosed (above-level) region.
    """
    tris, grads = _march(volume.values, level, normals=True)
    verts = tris.reshape(-1, 3)
    if len(verts) == 0:
        empty = np.zeros((0, 3), dtype=np.float32)
        return empty, empty.copy()
    g = grads.reshape(-1, 3)
    normals = -(g @ np.linalg.inv(volume.axes).T)   # index-space -> world
    norm = np.linalg.norm(normals, axis=1)
    bad = norm < 1e-12
    normals[bad] = (0.0, 0.0, 1.0)
    norm[bad] = 1.0
    normals /= norm[:, None]
    world = volume.grid_to_world(verts)
    return world.astype(np.float32), normals.astype(np.float32)
