"""Measurement/label annotations drawn with QPainter.

Shared between the live viewport overlay and Export Image, so figures carry
the same publication-style annotations as the screen:

- distances: value written parallel to the bond, centered on it
- angles: arc in the a-b-c plane with the value at the bisector
- dihedrals: rotation arrow around the central bond (Newman style)
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF

from ..core import geometry, measure
from ..core.molecule import Molecule
from ..render.camera import OrthoCamera

TEXT_COLOR = QColor(25, 25, 25)
HALO_COLOR = QColor(255, 255, 255)
PIN_COLOR = QColor(95, 95, 95)


def measurement_value(mol: Molecule, idxs: list[int]) -> str:
    """Kept as a name the painter can call; the wording lives in core.measure."""
    return measure.value(mol, idxs)


def overlay_metrics(height_px: int, half_height: float, scale: float = 1.0) -> dict:
    """Sizes that track the zoom level so text stays legible at any scale.

    *scale* raises the legibility clamps for high-resolution export renders.
    """
    px_per_ang = height_px / (2.0 * half_height)
    return {
        "label_px": int(np.clip(0.17 * px_per_ang, 9 * scale, 34 * scale)),
        "value_px": int(np.clip(0.20 * px_per_ang, 10 * scale, 40 * scale)),
        "marker_r": float(np.clip(0.14 * px_per_ang, 8 * scale, 60 * scale)),
        "offset": int(np.clip(0.09 * px_per_ang, 5 * scale, 24 * scale)),
        "line_w": float(np.clip(0.018 * px_per_ang, 1.3 * scale, 3.5 * scale)),
    }


def draw_halo_text(painter: QPainter, x: float, y: float, text: str,
                   centered: bool = False, angle_deg: float = 0.0,
                   dy: float = 0.0) -> None:
    """Text with a white halo; optionally centered on (x, y) and rotated.

    *dy* offsets the baseline in the rotated frame (negative = above the
    anchor), used to float a label just off the line it annotates."""
    painter.save()
    painter.translate(x, y)
    if angle_deg:
        painter.rotate(angle_deg)
    dx = -painter.fontMetrics().horizontalAdvance(text) / 2.0 if centered else 0.0
    painter.setPen(QPen(HALO_COLOR))
    for ox, oy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        painter.drawText(QPointF(dx + ox, dy + oy), text)
    painter.setPen(QPen(TEXT_COLOR))
    painter.drawText(QPointF(dx, dy), text)
    painter.restore()


def draw_atom_labels(painter: QPainter, mol: Molecule, pts: np.ndarray,
                     metrics: dict, label_mode: str) -> None:
    font = QFont()
    font.setPixelSize(metrics["label_px"])
    painter.setFont(font)
    fm = painter.fontMetrics()
    baseline_dy = (fm.ascent() - fm.descent()) / 2.0
    for i in range(mol.natoms):
        sym = mol.symbols[i]
        if label_mode == "element":
            text = sym
        elif label_mode == "number":
            text = str(i + 1)
        else:
            text = f"{sym}{i + 1}"
        # centered on the atom, not beside it
        draw_halo_text(painter, pts[i, 0], pts[i, 1] + baseline_dy, text,
                       centered=True)


def _polyline(painter: QPainter, pts: np.ndarray) -> None:
    painter.drawPolyline(QPolygonF([QPointF(x, y) for x, y in pts]))


def _project_arc(painter: QPainter, camera: OrthoCamera, arc3d: np.ndarray,
                 width: float, height: float, metrics: dict,
                 color: QColor) -> np.ndarray | None:
    if len(arc3d) == 0:
        return None
    arc = camera.project(arc3d, width, height)
    pen = QPen(color)
    pen.setWidthF(metrics["line_w"])
    painter.setPen(pen)
    _polyline(painter, arc)
    return arc


def _arrow_head(painter: QPainter, arc: np.ndarray, size: float,
                color: QColor) -> None:
    """Small filled triangle at the end of a projected arc."""
    tip = arc[-1]
    back = arc[-4] if len(arc) >= 4 else arc[0]
    t = tip - back
    n = np.linalg.norm(t)
    if n < 1e-6:
        return
    t = t / n
    p = np.array([-t[1], t[0]])
    b1 = tip - t * size + p * size * 0.45
    b2 = tip - t * size - p * size * 0.45
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    painter.drawPolygon(QPolygonF([QPointF(*tip), QPointF(*b1), QPointF(*b2)]))
    painter.restore()


def _label_outward(painter: QPainter, anchor2d: np.ndarray, center2d: np.ndarray,
                   text: str, metrics: dict) -> None:
    """Centered text near *anchor2d*, pushed away from *center2d* far enough
    that the text box (and its halo) clears the anchor point entirely."""
    d = anchor2d - center2d
    n = np.linalg.norm(d)
    d = d / n if n > 1e-6 else np.array([0.0, -1.0])
    fm = painter.fontMetrics()
    half_w = fm.horizontalAdvance(text) / 2.0
    half_h = fm.height() / 2.0
    push = metrics["offset"] * 0.35 + abs(d[0]) * half_w + abs(d[1]) * half_h
    pos = anchor2d + d * push
    baseline = pos[1] + (fm.ascent() - fm.descent()) / 2.0
    draw_halo_text(painter, pos[0], baseline, text, centered=True)


def draw_measurement(painter: QPainter, mol: Molecule, camera: OrthoCamera,
                     width: float, height: float, idxs: list[int],
                     metrics: dict, color: QColor = PIN_COLOR) -> None:
    """One pinned/selected measurement: guide lines plus the styled value."""
    coords = mol.coords[list(idxs)]
    pts = camera.project(coords, width, height)
    pen = QPen(color)
    pen.setWidthF(metrics["line_w"])
    pen.setStyle(Qt.PenStyle.DashLine)
    painter.setPen(pen)
    _polyline(painter, pts)

    font = QFont()
    font.setPixelSize(metrics["value_px"])
    painter.setFont(font)
    text = measurement_value(mol, idxs)

    if len(idxs) == 2:
        mid = (pts[0] + pts[1]) / 2.0
        d = pts[1] - pts[0]
        ang = np.degrees(np.arctan2(d[1], d[0])) if np.linalg.norm(d) > 1e-6 else 0.0
        if ang > 90.0:                  # keep the text upright
            ang -= 180.0
        elif ang <= -90.0:
            ang += 180.0
        draw_halo_text(painter, mid[0], mid[1], text, centered=True,
                       angle_deg=ang, dy=-metrics["offset"] * 0.9)
        return

    if len(idxs) == 3:
        arc3d = geometry.angle_arc_points(*coords)
        arc = _project_arc(painter, camera, arc3d, width, height, metrics, color)
        vertex = pts[1]
        anchor = arc[len(arc) // 2] if arc is not None else vertex
        _label_outward(painter, anchor, vertex, text, metrics)
        return

    arc3d = geometry.dihedral_arc_points(*coords)
    center = (pts[1] + pts[2]) / 2.0
    if len(arc3d):
        # solid arms from the bond midpoint to each side, so the arc visibly
        # connects the a-b and c-d bonds (exact overlay when viewed down b-c)
        for seg in geometry.dihedral_arm_points(*coords):
            _project_arc(painter, camera, seg, width, height, metrics, color)
    arc = _project_arc(painter, camera, arc3d, width, height, metrics, color)
    if arc is not None:
        _arrow_head(painter, arc, max(7.0, metrics["offset"] * 1.3), color)
        _label_outward(painter, arc[len(arc) // 2], center, text, metrics)
    else:
        _label_outward(painter, center + np.array([0.0, -metrics["offset"]]),
                       center, text, metrics)


AXIS_COLORS = (QColor(214, 57, 57), QColor(40, 158, 60), QColor(56, 108, 214))
AXIS_NAMES = ("X", "Y", "Z")


def draw_axis_indicator(painter: QPainter, camera: OrthoCamera,
                        width: float, height: float, scale: float = 1.0) -> None:
    """Corner triad showing how the world x/y/z axes currently point.

    Purely an orientation cue: it sits in a fixed screen corner at a fixed
    size, and only its directions follow the camera. Axes pointing away from
    the viewer are drawn first and faded, so the triad reads as 3-D.
    """
    arm = 34.0 * scale
    # kept clear of the bottom-left measurement HUD in the live viewport
    origin = np.array([arm + 22.0 * scale, height - arm - 40.0 * scale])
    font = QFont()
    font.setPixelSize(int(13 * scale))
    font.setBold(True)
    painter.setFont(font)
    # rotation maps world -> view, so column k is world axis k seen in view
    # space; view x/y are screen right/up and view z points at the viewer
    view_dirs = camera.rotation
    order = np.argsort(view_dirs[2])          # back to front
    for k in order:
        d = view_dirs[:, k]
        tip = origin + np.array([d[0], -d[1]]) * arm
        toward = float(d[2])
        color = QColor(AXIS_COLORS[k])
        if toward < 0.0:                      # pointing away: fade it back
            color.setAlpha(150)
        pen = QPen(color)
        pen.setWidthF(2.4 * scale)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(QPointF(*origin), QPointF(*tip))
        flat = float(np.hypot(d[0], d[1]))
        if flat > 0.12:                       # too edge-on to label readably
            label = origin + np.array([d[0], -d[1]]) * (arm + 9.0 * scale)
            painter.setPen(QPen(color))
            fm = painter.fontMetrics()
            painter.drawText(
                QPointF(label[0] - fm.horizontalAdvance(AXIS_NAMES[k]) / 2.0,
                        label[1] + (fm.ascent() - fm.descent()) / 2.0),
                AXIS_NAMES[k])
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor(90, 90, 90))
    painter.drawEllipse(QPointF(*origin), 2.6 * scale, 2.6 * scale)
    painter.setBrush(Qt.BrushStyle.NoBrush)


def draw_annotations(painter: QPainter, mol: Molecule, camera: OrthoCamera,
                     width: float, height: float, pinned: list[list[int]],
                     label_mode: str = "none", scale: float = 1.0,
                     axes: bool = False) -> None:
    """Everything Export Image needs on top of the raster render."""
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    metrics = overlay_metrics(int(height), camera.half_height, scale)
    if label_mode != "none":
        pts = camera.project(mol.coords, width, height)
        draw_atom_labels(painter, mol, pts, metrics, label_mode)
    for idxs in pinned:
        draw_measurement(painter, mol, camera, width, height, idxs, metrics)
    if axes:
        draw_axis_indicator(painter, camera, width, height, scale)
