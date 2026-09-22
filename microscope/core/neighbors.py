"""Pairs of points closer than a cutoff, found through a grid of cells.

Bond perception and hydrogen-bond detection both ask this of up to a hundred
thousand atoms. With cells one cutoff wide only neighbouring cells can hold a
close pair, so the points are sorted by cell and, for each of the thirteen
forward neighbour offsets and the cell itself, every point's partners are
found with one searchsorted over the whole structure. There is no Python loop
over cells: on a protein that loop was most of the time bond perception took
(1.3 s down to 0.2 for 98,000 atoms).
"""

from __future__ import annotations

import numpy as np

# half of the 27-cell neighbourhood, plus the cell itself: each pair of cells once
_FORWARD = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
            if (dx, dy, dz) > (0, 0, 0) or (dx, dy, dz) == (0, 0, 0)]


def close_pairs(coords: np.ndarray, cutoff: float,
                among: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(i, j, distance)`` for every pair with ``i < j`` no further apart than
    *cutoff*, sorted by ``i`` then ``j``. *among* restricts the search to those
    point indices; the returned indices are still into *coords*."""
    coords = np.asarray(coords, dtype=float)
    index = np.arange(len(coords)) if among is None else np.asarray(among, dtype=np.int64)
    empty = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
    if len(index) < 2 or cutoff <= 0:
        return empty
    points = coords[index]
    # cells numbered from 1 in a grid one cell wider on every side, so that a
    # neighbour offset never wraps round into the next row
    cells = np.floor((points - points.min(axis=0)) / cutoff).astype(np.int64) + 1
    dims = cells.max(axis=0) + 2
    key = (cells[:, 0] * dims[1] + cells[:, 1]) * dims[2] + cells[:, 2]
    order = np.argsort(key, kind="stable")
    sorted_key = key[order]
    found_i, found_j, found_d = [], [], []
    for dx, dy, dz in _FORWARD:
        shift = (dx * dims[1] + dy) * dims[2] + dz
        target = key + shift
        low = np.searchsorted(sorted_key, target, "left")
        count = np.searchsorted(sorted_key, target, "right") - low
        if not count.any():
            continue
        a = np.repeat(np.arange(len(key)), count)
        # each point's partners are a run of the sorted order starting at `low`
        runs = np.arange(len(a)) - np.repeat(np.cumsum(count) - count, count)
        b = order[np.repeat(low, count) + runs]
        if shift == 0:
            keep = a < b                      # within one cell, each pair once
            a, b = a[keep], b[keep]
        d = np.linalg.norm(points[a] - points[b], axis=1)
        near = d <= cutoff
        a, b = index[a[near]], index[b[near]]
        found_i.append(np.minimum(a, b))
        found_j.append(np.maximum(a, b))
        found_d.append(d[near])
    if not found_i:
        return empty
    i, j, d = np.concatenate(found_i), np.concatenate(found_j), np.concatenate(found_d)
    ordered = np.lexsort((j, i))
    return i[ordered], j[ordered], d[ordered]
