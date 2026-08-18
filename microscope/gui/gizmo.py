"""Interactive translate/rotate manipulator ("gizmo") for the selected atoms.

Modelled on the LightCone web tool: as soon as atoms are selected, three
coloured axis arrows and three rotation rings appear at the centroid of the
selection. Dragging an arrow slides the selection along that axis, dragging a
ring turns it about that axis, and dragging the centre dot moves it freely in
the screen plane.

The camera is orthographic, so a world direction projects to a fixed screen
direction: hit-testing is plain 2-D pixel geometry and the drag amounts come
out of a least-squares projection of the mouse movement onto the projected
axis. Drawing is QPainter overlay work, which also gives LightCone's
"always on top of the molecule" look for free.

Handles are identified by a ``(kind, index)`` pair: ``("axis", k)`` and
``("ring", k)`` for k = 0, 1, 2 (world X, Y, Z) and ``("center", -1)``.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPen, QPolygonF

from ..render.camera import OrthoCamera

AXIS_NAMES = ("X", "Y", "Z")
AXIS_VECTORS = np.eye(3)
AXIS_COLORS = (QColor(214, 57, 57), QColor(40, 158, 60), QColor(56, 108, 214))
RING_COLORS = (QColor(226, 124, 124), QColor(110, 190, 122), QColor(122, 156, 226))
CENTER_COLOR = QColor(70, 70, 70)

AXIS_LENGTH_PX = 92.0     # on-screen arrow length, so the gizmo never rescales
RING_RADIUS_FRAC = 0.62   # ring radius as a fraction of the arrow length
PICK_TOL_PX = 8.0
CENTER_RADIUS_PX = 6.0
MIN_AXIS_PX = 14.0        # axis pointing at the viewer: too foreshortened to use

CENTER = ("center", -1)


# --------------------------------------------------------------------- geometry

def selection_center(coords: np.ndarray, indices) -> np.ndarray | None:
    """Centroid of the selected atoms — where the manipulator sits."""
    idx = [int(i) for i in indices]
    if not idx:
        return None
    return np.asarray(coords)[idx].mean(axis=0)


def world_length(camera: OrthoCamera, height_px: int,
                 pixels: float = AXIS_LENGTH_PX) -> float:
    """World length that covers *pixels* on screen at the current zoom."""
    return pixels * 2.0 * camera.half_height / max(height_px, 1)


def ring_points(center: np.ndarray, k: int, radius: float,
                n: int = 72) -> np.ndarray:
    """Closed circle of radius *radius* around *center*, normal to axis *k*."""
    u = AXIS_VECTORS[(k + 1) % 3]
    v = AXIS_VECTORS[(k + 2) % 3]
    t = np.linspace(0.0, 2.0 * np.pi, n + 1)
    return center + radius * (np.outer(np.cos(t), u) + np.outer(np.sin(t), v))


def _segment_distance(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    d = b - a
    dd = float(d @ d)
    if dd < 1e-12:
        return float(np.linalg.norm(p - a))
    t = float(np.clip((p - a) @ d / dd, 0.0, 1.0))
    return float(np.linalg.norm(p - (a + t * d)))


def _polyline_distance(p: np.ndarray, pts: np.ndarray) -> float:
    seg = pts[1:] - pts[:-1]
    dd = np.einsum("ij,ij->i", seg, seg)
    dd[dd < 1e-12] = 1e-12
    t = np.clip(np.einsum("ij,ij->i", p - pts[:-1], seg) / dd, 0.0, 1.0)
    closest = pts[:-1] + t[:, None] * seg
    return float(np.linalg.norm(closest - p, axis=1).min())


def _screen(camera: OrthoCamera, points, width: float, height: float) -> np.ndarray:
    return camera.project(np.atleast_2d(points), width, height)


def hit_test(camera: OrthoCamera, width: float, height: float,
             center: np.ndarray, mouse, length: float | None = None):
    """Handle under the mouse, or None. Arrows win over rings over the ring
    interior, matching the drawing order (arrows sit on top)."""
    if center is None:
        return None
    if length is None:
        length = world_length(camera, int(height))
    p = np.asarray(mouse, dtype=float)
    center2d = _screen(camera, center, width, height)[0]

    if float(np.linalg.norm(p - center2d)) <= CENTER_RADIUS_PX + 2.0:
        return CENTER

    tips = _screen(camera, [center + AXIS_VECTORS[k] * length for k in range(3)],
                   width, height)
    best, best_d = None, PICK_TOL_PX
    for k in range(3):
        if float(np.linalg.norm(tips[k] - center2d)) < MIN_AXIS_PX:
            continue          # edge-on axis: unusable, and it would eat clicks
        d = _segment_distance(p, center2d, tips[k])
        if d < best_d:
            best, best_d = ("axis", k), d
    if best is not None:
        return best

    radius = length * RING_RADIUS_FRAC
    for k in range(3):
        ring = _screen(camera, ring_points(center, k, radius), width, height)
        d = _polyline_distance(p, ring)
        if d < best_d:
            best, best_d = ("ring", k), d
    return best


# --------------------------------------------------------------------- drag math

def axis_translation(camera: OrthoCamera, width: float, height: float,
                     center: np.ndarray, k: int, delta_px) -> float:
    """World distance to slide along axis *k* for a mouse movement of
    *delta_px*, by projecting the movement onto the axis' screen direction."""
    pts = _screen(camera, [center, center + AXIS_VECTORS[k]], width, height)
    s = pts[1] - pts[0]                     # pixels per world unit along the axis
    denom = float(s @ s)
    if denom < 1e-6:
        return 0.0
    return float(np.asarray(delta_px, dtype=float) @ s / denom)


def plane_translation(camera: OrthoCamera, height: float, delta_px) -> np.ndarray:
    """World displacement for a free drag in the screen plane."""
    dx, dy = (float(v) for v in delta_px)
    per_px = 2.0 * camera.half_height / max(height, 1)
    right, up = camera.rotation[0], camera.rotation[1]
    return right * (dx * per_px) - up * (dy * per_px)


def screen_angle(camera: OrthoCamera, width: float, height: float,
                 center: np.ndarray, mouse) -> float:
    """Angle of the mouse around the projected centre, in radians, CCW on
    screen (screen y grows downward, so it is flipped here)."""
    center2d = _screen(camera, center, width, height)[0]
    d = np.asarray(mouse, dtype=float) - center2d
    return float(np.arctan2(-d[1], d[0]))


def rotation_sign(camera: OrthoCamera, k: int) -> float:
    """+1 if axis *k* points towards the viewer, so a right-hand-rule rotation
    looks counter-clockwise on screen; -1 when it points away."""
    toward_viewer = camera.rotation[2]      # world direction of the view +z axis
    return 1.0 if float(AXIS_VECTORS[k] @ toward_viewer) >= 0.0 else -1.0


def wrap_angle(delta: float) -> float:
    """Angle difference folded into (-pi, pi] so dragging past ±180° works."""
    return (delta + np.pi) % (2.0 * np.pi) - np.pi


# --------------------------------------------------------------------- drawing

def _arrow_head(painter, tip: np.ndarray, direction: np.ndarray,
                size: float, color: QColor) -> None:
    n = float(np.linalg.norm(direction))
    if n < 1e-6:
        return
    t = direction / n
    perp = np.array([-t[1], t[0]])
    b1 = tip - t * size + perp * size * 0.42
    b2 = tip - t * size - perp * size * 0.42
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    painter.drawPolygon(QPolygonF([QPointF(*tip), QPointF(*b1), QPointF(*b2)]))
    painter.restore()


def draw(painter, camera: OrthoCamera, width: float, height: float,
         center: np.ndarray, hovered=None, active=None,
         length: float | None = None) -> None:
    """Paint the manipulator. *hovered*/*active* thicken the handle under the
    mouse and the one being dragged (LightCone highlights the same way)."""
    if center is None:
        return
    if length is None:
        length = world_length(camera, int(height))
    lit = active or hovered
    center2d = _screen(camera, center, width, height)[0]
    radius = length * RING_RADIUS_FRAC

    for k in range(3):                                  # rings first: behind
        ring = _screen(camera, ring_points(center, k, radius), width, height)
        pen = QPen(RING_COLORS[k])
        pen.setWidthF(3.0 if lit == ("ring", k) else 1.6)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPolyline(QPolygonF([QPointF(x, y) for x, y in ring]))

    font = QFont()
    font.setPixelSize(12)
    font.setBold(True)
    painter.setFont(font)
    tips = _screen(camera, [center + AXIS_VECTORS[k] * length for k in range(3)],
                   width, height)
    for k in range(3):
        d = tips[k] - center2d
        if float(np.linalg.norm(d)) < MIN_AXIS_PX:
            continue
        wide = lit == ("axis", k)
        pen = QPen(AXIS_COLORS[k])
        pen.setWidthF(3.4 if wide else 2.2)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        head = 11.0 if wide else 9.0
        shaft = tips[k] - d / np.linalg.norm(d) * head
        painter.drawLine(QPointF(*center2d), QPointF(*shaft))
        _arrow_head(painter, tips[k], d, head, AXIS_COLORS[k])
        label = center2d + d * (1.0 + 17.0 / max(float(np.linalg.norm(d)), 1e-6))
        painter.setPen(QPen(AXIS_COLORS[k]))
        painter.drawText(QPointF(label[0] - 4.0, label[1] + 4.0), AXIS_NAMES[k])

    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(CENTER_COLOR if lit == CENTER else QColor(120, 120, 120))
    r = CENTER_RADIUS_PX + (1.5 if lit == CENTER else 0.0)
    painter.drawEllipse(QPointF(*center2d), r, r)
