"""Paint an orbital energy-level diagram with QPainter.

The dialog draws it next to the orbital list, and the same function writes it
to PNG, SVG or PDF - it is line art, like a Lewis structure, so it stays line
art. Levels are placed by render.levels; this only turns them into pixels.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QImage, QPainter, QPen

from .formats import VECTOR_SUFFIXES
from .imagefile import vector_device, write_image
from .levels import Level, LevelDiagram, side_by_side

INK = QColor(40, 40, 40)
EMPTY = QColor(150, 150, 150)       # virtual orbitals
PICKED = QColor(224, 130, 20)       # the orbital on show, as the spectra mark theirs
LINE = 34.0                         # length of a lone level, px
MARGIN = (58.0, 18.0, 6.0, 22.0)    # left (energy axis), top, right, bottom
CENTRE = 0.3                        # levels sit this far into their column; names
                                    # and the gap take the rest, so columns never meet
ARROW = 6.0                         # half the height of an electron's arrow, px


def draw_levels(painter: QPainter, diagram: LevelDiagram, width: float, height: float,
                picked: tuple[str, int] | None = None,
                font_px: float = 12.0) -> dict[tuple[str, int], QRectF]:
    """Draw *diagram* into a *width* x *height* area; returns where each level
    is, for clicks to be matched against."""
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    font = QFont()
    font.setPixelSize(max(int(font_px), 6))
    painter.setFont(font)
    metrics = QFontMetricsF(font)
    left, top, right, bottom = MARGIN
    if not diagram.levels:
        painter.setPen(EMPTY)
        painter.drawText(QRectF(0, 0, width, height), Qt.AlignmentFlag.AlignCenter,
                         "no orbital energies in this file")
        return {}

    low, high = diagram.span
    pad = 0.06 * max(high - low, 1.0)
    low, high = low - pad, high + pad
    header = metrics.height() + 4 if diagram.columns > 1 else 0.0
    top += header

    def y_of(energy: float) -> float:
        return top + (high - energy) / (high - low) * (height - top - bottom)

    _axis(painter, metrics, left, top, height - bottom, low, high, y_of)
    column_width = (width - left - right) / diagram.columns
    boxes: dict[tuple[str, int], QRectF] = {}
    for column in range(diagram.columns):
        middle = left + (column + CENTRE) * column_width
        if diagram.columns > 1:
            painter.setPen(INK)
            painter.drawText(QRectF(middle - 30, top - header, 60, header),
                             Qt.AlignmentFlag.AlignCenter, "α" if column == 0 else "β")
        levels = [lv for lv in diagram.levels if lv.column == column]
        heights = [y_of(lv.energy) for lv in levels]
        placed = side_by_side(heights, 2 * ARROW + 3)
        ends: dict[int, float] = {}                 # right end of each level's group
        for lv, y, (slot, partners) in zip(levels, heights, placed):
            box, end = _level(painter, lv, middle, y, slot, partners, column_width,
                              picked)
            boxes[(lv.spin, lv.index)] = box
            ends[lv.index] = end
        _frontier(painter, metrics, levels, diagram, ends, y_of,
                  column_end=left + (column + 1) * column_width)
    return boxes


def _axis(painter, metrics, x, y0, y1, low, high, y_of) -> None:
    painter.setPen(QPen(INK, 1.0))
    painter.drawLine(QPointF(x, y0), QPointF(x, y1))
    step = _nice_step((high - low) / 6)
    for value in np.arange(np.ceil(low / step) * step, high, step):
        y = y_of(value)
        painter.drawLine(QPointF(x - 4, y), QPointF(x, y))
        text = f"{value:g}".replace("-", "−")
        painter.drawText(QRectF(0, y - 9, x - 7, 18),
                         Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, text)
    painter.save()
    painter.translate(10, (y0 + y1) / 2)
    painter.rotate(-90)
    painter.drawText(QRectF(-60, -10, 120, 20), Qt.AlignmentFlag.AlignCenter, "E (eV)")
    painter.restore()


def _nice_step(rough: float) -> float:
    power = 10 ** np.floor(np.log10(max(rough, 1e-6)))
    return float(next(f * power for f in (1, 2, 5, 10) if f * power >= rough))


def _level(painter, lv: Level, middle: float, y: float, slot: int, partners: int,
           column_width: float, picked) -> tuple[QRectF, float]:
    """One level with its electrons, its place among partners that share the
    height; returns where it is and where its group ends on the right."""
    length = min(LINE, 0.45 * column_width / partners)
    gap = 0.25 * length
    total = partners * length + (partners - 1) * gap
    x0 = middle - total / 2 + slot * (length + gap)
    chosen = picked == (lv.spin, lv.index)
    colour = PICKED if chosen else INK if lv.electrons else EMPTY
    pen = QPen(colour, 3.0 if chosen else 2.0)
    pen.setCapStyle(Qt.PenCapStyle.FlatCap)
    painter.setPen(pen)
    painter.drawLine(QPointF(x0, y), QPointF(x0 + length, y))
    # electrons as arrows straddling the level, as a textbook draws them
    ups = ([True, False] if lv.electrons == 2 else [lv.spin != "beta"] if lv.electrons
           else [])
    painter.setPen(QPen(colour, 1.2))
    for n, up in enumerate(ups):
        x = x0 + length * (n + 1) / (len(ups) + 1)
        tip, tail = (y - ARROW, y + ARROW) if up else (y + ARROW, y - ARROW)
        head = 3.0 if up else -3.0
        painter.drawLine(QPointF(x, tail), QPointF(x, tip))
        painter.drawLine(QPointF(x, tip), QPointF(x - 2.5, tip + head))
        painter.drawLine(QPointF(x, tip), QPointF(x + 2.5, tip + head))
    return QRectF(x0 - 2, y - 7, length + 4, 14), middle + total / 2


def _frontier(painter, metrics, levels: list[Level], diagram: LevelDiagram,
              ends: dict[int, float], y_of, column_end: float) -> None:
    """HOMO and LUMO named beside their groups, and the gap between them."""
    painter.setPen(INK)
    for lv in levels:
        if lv.name in ("HOMO", "LUMO"):
            painter.drawText(QRectF(ends[lv.index] + 6, y_of(lv.energy) - 9, 60, 18),
                             Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                             lv.name)
    spin = levels[0].spin if levels else ""
    homo = next((lv for lv in levels if lv.name == "HOMO"), None)
    lumo = next((lv for lv in levels if lv.name == "LUMO"), None)
    if spin in diagram.gaps and homo is not None and lumo is not None:
        x = max(ends[homo.index], ends[lumo.index]) + 46
        y_homo, y_lumo = y_of(homo.energy), y_of(lumo.energy)
        painter.setPen(QPen(EMPTY, 1.0, Qt.PenStyle.DashLine))
        painter.drawLine(QPointF(x, y_homo), QPointF(x, y_lumo))
        painter.setPen(INK)
        text = f"{diagram.gaps[spin]:.2f} eV"
        wide = metrics.horizontalAdvance(text) + 2
        # beside the line, on whichever side the column still has room
        x_text = x + 5 if x + 5 + wide <= column_end else x - 5 - wide
        painter.drawText(QRectF(x_text, (y_homo + y_lumo) / 2 - 9, wide, 18),
                         Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, text)


def write_levels(path, diagram: LevelDiagram, width: int | None = None,
                 height: int = 480, picked: tuple[str, int] | None = None) -> None:
    """Write the diagram as SVG or PDF (points), or as a PNG/TIFF at 3x."""
    width = width or (280 if diagram.columns == 1 else 460)
    if Path(path).suffix.lower() in VECTOR_SUFFIXES:
        device = vector_device(path, width, height, "Orbital energy levels")
        painter = QPainter(device)          # the device must outlive the painter
        try:
            painter.fillRect(0, 0, width, height, Qt.GlobalColor.white)
            draw_levels(painter, diagram, width, height, picked)
        finally:
            painter.end()
        return
    scale = 3
    image = QImage(width * scale, height * scale, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(Qt.GlobalColor.white)
    painter = QPainter(image)
    try:
        painter.scale(scale, scale)
        draw_levels(painter, diagram, width, height, picked)
    finally:
        painter.end()
    write_image(image, str(path))
