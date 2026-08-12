"""Orthographic trackball camera."""

from __future__ import annotations

import math

import numpy as np

from ..core.geometry import rotation_matrix


class OrthoCamera:
    def __init__(self):
        self.center = np.zeros(3)
        self.rotation = np.eye(3)      # world -> view
        self.half_height = 10.0
        self.distance = 30.0
        self.near = 0.5
        self.far = 100.0
        self.scene_radius = 10.0

    def fit(self, center: np.ndarray, radius: float) -> None:
        radius = max(float(radius), 1.0)
        self.center = np.asarray(center, dtype=float).copy()
        self.scene_radius = radius
        self.half_height = radius * 1.1 + 0.5
        self.distance = radius * 3.0 + 5.0
        self.near = 0.5
        self.far = self.distance + radius * 3.0 + 5.0

    def reset_orientation(self) -> None:
        self.rotation = np.eye(3)

    def view_matrix(self) -> np.ndarray:
        V = np.eye(4)
        V[:3, :3] = self.rotation
        V[:3, 3] = -self.rotation @ self.center
        V[2, 3] -= self.distance
        return V

    def proj_matrix(self, aspect: float) -> np.ndarray:
        h = self.half_height
        w = h * max(aspect, 1e-6)
        n, f = self.near, self.far
        P = np.zeros((4, 4))
        P[0, 0] = 1.0 / w
        P[1, 1] = 1.0 / h
        P[2, 2] = -2.0 / (f - n)
        P[2, 3] = -(f + n) / (f - n)
        P[3, 3] = 1.0
        return P

    def rotate_drag(self, dx: float, dy: float, speed: float = 0.008) -> None:
        angle = speed * math.hypot(dx, dy)
        if angle < 1e-9:
            return
        axis = np.array([dy, dx, 0.0])
        self.rotation = rotation_matrix(axis, angle) @ self.rotation

    def pan_drag(self, dx: float, dy: float, viewport_height_px: int) -> None:
        world_per_px = 2.0 * self.half_height / max(viewport_height_px, 1)
        right = self.rotation[0]   # world direction of the view x axis
        up = self.rotation[1]
        self.center = self.center - right * (dx * world_per_px) + up * (dy * world_per_px)

    def zoom(self, steps: float) -> None:
        self.half_height *= 0.9 ** steps
        lo = self.scene_radius * 0.02 + 1e-3
        hi = self.scene_radius * 20.0 + 10.0
        self.half_height = float(np.clip(self.half_height, lo, hi))

    def project(self, points: np.ndarray, width: float, height: float) -> np.ndarray:
        """Project world points to widget pixel coordinates (origin top-left)."""
        points = np.atleast_2d(points)
        M = self.proj_matrix(width / max(height, 1)) @ self.view_matrix()
        homo = np.hstack([points, np.ones((len(points), 1))])
        ndc = (M @ homo.T).T
        sx = (ndc[:, 0] + 1.0) * 0.5 * width
        sy = (1.0 - ndc[:, 1]) * 0.5 * height
        return np.column_stack([sx, sy])


def orientation_along(direction: np.ndarray) -> np.ndarray | None:
    """World->view rotation whose view z axis (toward the viewer) is *direction*.

    Used to look straight down a bond. Returns None for a degenerate direction.
    """
    z = np.asarray(direction, dtype=float)
    norm = np.linalg.norm(z)
    if norm < 1e-9:
        return None
    z = z / norm
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = np.cross(helper, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.array([x, y, z])


def orientation_from_plane(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray | None:
    """World->view rotation that puts the plane of points a, b, c in the screen
    plane, with the a->b direction horizontal. Returns None if collinear."""
    a, b, c = (np.asarray(p, dtype=float) for p in (a, b, c))
    n = np.cross(b - a, c - a)
    nn = np.linalg.norm(n)
    if nn < 1e-9:
        return None
    n = n / nn
    x = b - a
    x = x - np.dot(x, n) * n
    xn = np.linalg.norm(x)
    if xn < 1e-9:
        return None
    x = x / xn
    y = np.cross(n, x)
    return np.array([x, y, n])
