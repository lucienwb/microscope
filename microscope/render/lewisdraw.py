"""Painting a Lewis structure with QPainter.

Shared by the interactive Lewis view, Export Image and `scope -s --lewis`, so
a saved figure is exactly what was on screen. Drawing needs real font metrics
— a label's box is what the bonds stop short of — so the layout is computed
here with the painter's own font and then drawn.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QMarginsF, QPointF, QRect, QSize, QSizeF, Qt
from PySide6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontMetricsF,
    QImage,
    QPageSize,
    QPainter,
    QPdfWriter,
    QPen,
)

from ..core.lewis import LewisStructure
from .camera import OrthoCamera
from .formats import VECTOR_SUFFIXES
from .lewis2d import LewisLayout, LewisOptions, layout

# ChemDraw sets structures in Arial; Qt substitutes on machines without it.
FONT_FAMILIES = ("Arial", "Helvetica", "Helvetica Neue", "DejaVu Sans")
BOND_COLOR = QColor(0, 0, 0)
SUBSCRIPT_SCALE = 0.68
SUBSCRIPT_DROP = 0.24       # how far a subscript sits below the baseline
CHARGE_SCALE = 0.72


def _fonts(font_px: float) -> tuple[QFont, QFont]:
    main = QFont()
    main.setFamilies(list(FONT_FAMILIES))
    main.setPixelSize(max(1, int(round(font_px))))
    sub = QFont(main)
    sub.setPixelSize(max(1, int(round(font_px * SUBSCRIPT_SCALE))))
    return main, sub


def measure(runs, font_px: float) -> tuple[float, float]:
    """Pixel extent of a label, measured with the font it will be drawn in."""
    main, sub = _fonts(font_px)
    normal, small = QFontMetricsF(main), QFontMetricsF(sub)
    width = sum((small if is_sub else normal).horizontalAdvance(text)
                for text, is_sub in runs)
    return width, normal.capHeight()


def _qcolor(rgb) -> QColor:
    return QColor.fromRgbF(*(float(c) for c in rgb))


def draw_layout(painter: QPainter, plan: LewisLayout) -> None:
    """Draw a finished layout: bonds, then labels over them, then electrons."""
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)

    pen = QPen(BOND_COLOR)
    pen.setWidthF(plan.line_width)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    painter.setPen(pen)
    for bond in plan.bonds:
        for start, end in bond.lines:
            painter.drawLine(QPointF(*start), QPointF(*end))

    main, sub = _fonts(plan.font_px)
    normal, small = QFontMetricsF(main), QFontMetricsF(sub)
    charge_font = QFont(main)
    charge_font.setPixelSize(max(1, int(round(plan.font_px * CHARGE_SCALE))))
    dot_radius = plan.line_width * 1.35

    if plan.brackets is not None:
        _draw_brackets(painter, plan.brackets, plan.line_width, charge_font)

    for atom in plan.atoms:
        painter.setPen(QPen(_qcolor(atom.color)))
        if atom.labelled:
            width = sum((small if is_sub else normal).horizontalAdvance(text)
                        for text, is_sub in atom.runs)
            x = float(atom.pos[0]) - width / 2.0
            baseline = float(atom.pos[1]) + normal.capHeight() / 2.0
            for text, is_sub in atom.runs:
                painter.setFont(sub if is_sub else main)
                dy = plan.font_px * SUBSCRIPT_DROP if is_sub else 0.0
                painter.drawText(QPointF(x, baseline + dy), text)
                x += (small if is_sub else normal).horizontalAdvance(text)
        if atom.charge_text and atom.charge_pos is not None:
            painter.setFont(charge_font)
            metrics = QFontMetricsF(charge_font)
            painter.drawText(
                QPointF(float(atom.charge_pos[0])
                        - metrics.horizontalAdvance(atom.charge_text) / 2.0,
                        float(atom.charge_pos[1]) + metrics.capHeight() / 2.0),
                atom.charge_text)
        if atom.dots:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(_qcolor(atom.color)))
            for dot in atom.dots:
                painter.drawEllipse(QPointF(*dot), dot_radius, dot_radius)
            painter.setBrush(Qt.BrushStyle.NoBrush)


def _draw_brackets(painter: QPainter, box, line_width: float,
                   font: QFont) -> None:
    """Square brackets round the drawing, with the species charge outside."""
    pen = QPen(BOND_COLOR)
    pen.setWidthF(line_width)
    pen.setCapStyle(Qt.PenCapStyle.SquareCap)
    painter.setPen(pen)
    for x, arm in ((box.left, box.arm), (box.right, -box.arm)):
        painter.drawLine(QPointF(x, box.top), QPointF(x, box.bottom))
        painter.drawLine(QPointF(x, box.top), QPointF(x + arm, box.top))
        painter.drawLine(QPointF(x, box.bottom), QPointF(x + arm, box.bottom))

    if not box.runs or box.label_pos is None:
        return
    painter.setFont(font)
    metrics = QFontMetricsF(font)
    text = "".join(t for t, _ in box.runs)
    painter.drawText(
        QPointF(float(box.label_pos[0]) - metrics.horizontalAdvance(text) / 2.0,
                float(box.label_pos[1]) + metrics.capHeight() / 2.0), text)


def draw_lewis(painter: QPainter, structure: LewisStructure, camera: OrthoCamera,
               width: float, height: float,
               options: LewisOptions | None = None) -> LewisLayout:
    """Lay the structure out for this canvas and paint it; returns the layout
    so callers can report its warnings."""
    plan = layout(structure, camera, width, height, options, measure=measure)
    draw_layout(painter, plan)
    return plan


def render_lewis_image(structure: LewisStructure, camera: OrthoCamera,
                       width: int, height: int,
                       options: LewisOptions | None = None,
                       transparent: bool = True,
                       background=(1.0, 1.0, 1.0)) -> QImage:
    """Draw the structure into a QImage. No OpenGL is involved — a Lewis
    figure is line art, so it is painted straight onto the image."""
    image = QImage(int(width), int(height), QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(Qt.GlobalColor.transparent if transparent else _qcolor(background))
    painter = QPainter(image)
    try:
        draw_lewis(painter, structure, camera, int(width), int(height), options)
    finally:
        painter.end()
    return image


# A Lewis structure is line art, so it is worth keeping as line art: SVG for
# an illustration program, PDF for a manuscript.


def write_lewis_vector(path, structure: LewisStructure, camera: OrthoCamera,
                       width: int, height: int,
                       options: LewisOptions | None = None,
                       background=None) -> None:
    """Write the drawing to SVG or PDF at *width* x *height* points."""
    suffix = Path(path).suffix.lower()
    if suffix == ".svg":
        device = _svg_device(path, width, height)
    elif suffix == ".pdf":
        device = _pdf_device(path, width, height)
    else:
        raise ValueError(f"cannot write {suffix or 'that'} as vector art "
                         f"(use {' or '.join(VECTOR_SUFFIXES)})")
    painter = QPainter(device)
    try:
        if background is not None:
            painter.fillRect(0, 0, int(width), int(height), _qcolor(background))
        draw_lewis(painter, structure, camera, int(width), int(height), options)
    finally:
        painter.end()


def _svg_device(path, width: int, height: int):
    from PySide6.QtSvg import QSvgGenerator  # QtSvg is a separate module

    generator = QSvgGenerator()
    generator.setFileName(str(path))
    generator.setSize(QSize(int(width), int(height)))
    generator.setViewBox(QRect(0, 0, int(width), int(height)))
    generator.setTitle("Lewis structure")
    generator.setDescription("Drawn by microscope")
    return generator


def _pdf_device(path, width: int, height: int):
    writer = QPdfWriter(str(path))
    writer.setResolution(72)                          # one unit is one point
    writer.setPageSize(QPageSize(QSizeF(width, height), QPageSize.Unit.Point,
                                 "figure", QPageSize.SizeMatchPolicy.ExactMatch))
    writer.setPageMargins(QMarginsF(0, 0, 0, 0))
    return writer
