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
