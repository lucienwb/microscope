"""Geometric measurements and transforms."""

from __future__ import annotations

import numpy as np


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(b) - np.asarray(a)))


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle a-b-c in degrees."""
    u = np.asarray(a) - np.asarray(b)
    v = np.asarray(c) - np.asarray(b)
    cosang = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))


def dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    """Dihedral a-b-c-d in degrees, in (-180, 180]."""
    a, b, c, d = (np.asarray(p, dtype=float) for p in (a, b, c, d))
    b0, b1, b2 = b - a, c - b, d - c
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(y, x)))


def angle_arc_points(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                     radius: float | None = None, segments: int = 24) -> np.ndarray:
    """Arc marking the angle a-b-c: points in the a,b,c plane, centered at b,
    sweeping from the b->a direction to the b->c direction (the interior angle).
    Returns (segments+1, 3) world coordinates, or (0, 3) if degenerate."""
    a, b, c = (np.asarray(p, dtype=float) for p in (a, b, c))
    u, v = a - b, c - b
    lu, lv = np.linalg.norm(u), np.linalg.norm(v)
    if lu < 1e-9 or lv < 1e-9:
        return np.zeros((0, 3))
    e1, ec = u / lu, v / lv
    n = np.cross(e1, ec)
    nn = np.linalg.norm(n)
    if nn < 1e-9:                       # collinear: no plane to draw in
        return np.zeros((0, 3))
    e2 = np.cross(n / nn, e1)           # in-plane, perpendicular to e1, on c's side
    theta = np.arccos(np.clip(np.dot(e1, ec), -1.0, 1.0))
    if radius is None:
        radius = float(np.clip(0.33 * min(lu, lv), 0.15, 0.55))
    phi = np.linspace(0.0, theta, segments + 1)[:, None]
    return b + radius * (np.cos(phi) * e1 + np.sin(phi) * e2)


def _dihedral_frame(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                    d: np.ndarray) -> tuple | None:
    """(u, pa, pd, la, ld) of the a-b-c-d dihedral: unit bond axis and the
    outer-bond components perpendicular to it; None if degenerate."""
    axis = c - b
    lax = np.linalg.norm(axis)
    if lax < 1e-9:
        return None
    u = axis / lax
    pa = (a - b) - np.dot(a - b, u) * u
    pd = (d - c) - np.dot(d - c, u) * u
    la, ld = np.linalg.norm(pa), np.linalg.norm(pd)
    if la < 1e-9 or ld < 1e-9:          # an outer atom sits on the bond axis
        return None
    return u, pa, pd, la, ld


def dihedral_arc_points(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray,
                        radius: float | None = None, segments: int = 32) -> np.ndarray:
    """Rotation arrow for the dihedral a-b-c-d: an arc around the central b-c
    bond, at its midpoint, in the plane perpendicular to the bond, sweeping
    from the a side to the d side (Newman-projection style; follow the points
    to know the twist direction). Its endpoints lie on the arms returned by
    dihedral_arm_points. Returns (segments+1, 3), or (0, 3) if degenerate or
    the dihedral is ~0."""
    a, b, c, d = (np.asarray(p, dtype=float) for p in (a, b, c, d))
    frame = _dihedral_frame(a, b, c, d)
    if frame is None:
        return np.zeros((0, 3))
    u, pa, pd, la, ld = frame
    e1, f = pa / la, pd / ld
    phi = np.arctan2(np.dot(np.cross(e1, f), u), np.dot(e1, f))
    if abs(phi) < 1e-6:
        return np.zeros((0, 3))
    e2 = np.cross(u, e1)
    if radius is None:
        radius = float(np.clip(0.45 * min(la, ld), 0.25, 0.8))
    s = np.linspace(0.0, phi, segments + 1)[:, None]
    mid = 0.5 * (b + c)
    return mid + radius * (np.cos(s) * e1 + np.sin(s) * e2)


def dihedral_arm_points(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                        d: np.ndarray) -> np.ndarray:
    """The two straight arms the dihedral arc spans between: segments from the
    b-c bond midpoint to the projections of a and d onto the plane of the arc.
    Viewed down the bond they overlay the a-b and c-d bonds, so the arc
    visibly connects the two sides of the dihedral. Returns (2, 2, 3) as
    [[mid, a side], [mid, d side]], or (0, 2, 3) if degenerate."""
    a, b, c, d = (np.asarray(p, dtype=float) for p in (a, b, c, d))
    frame = _dihedral_frame(a, b, c, d)
    if frame is None:
        return np.zeros((0, 2, 3))
    u, pa, pd, la, ld = frame
    mid = 0.5 * (b + c)
    return np.array([[mid, mid + pa], [mid, mid + pd]])


def rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation matrix about *axis* by *angle_rad*."""
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return np.eye(3)
    x, y, z = axis / norm
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    t = 1.0 - c
    return np.array([
        [t * x * x + c,     t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c,     t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ])
