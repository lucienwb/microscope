"""The orbital energy-level diagram beside the orbital list: click a level to
draw that orbital. Drawn on white whatever the window theme, since it is a
preview of the figure it exports."""

from __future__ import annotations

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QPainter
from PySide6.QtWidgets import QWidget

from ..render.levels import LevelDiagram
from ..render.levelsdraw import draw_levels


class LevelsView(QWidget):
    picked = Signal(str, int)            # spin, 0-based orbital index

    def __init__(self, parent=None):
        super().__init__(parent)
        self.diagram = LevelDiagram()
        self.chosen: tuple[str, int] | None = None
        self._boxes: dict = {}
        self.setMinimumSize(250, 320)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("Click a level to draw that orbital")

    def show_diagram(self, diagram: LevelDiagram, chosen: tuple[str, int] | None) -> None:
        self.diagram, self.chosen = diagram, chosen
        # alpha and beta side by side, each with its names and gap beside it
        self.setMinimumWidth(250 if diagram.columns == 1 else 420)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        try:
            painter.fillRect(self.rect(), Qt.GlobalColor.white)
            self._boxes = draw_levels(painter, self.diagram, self.width(), self.height(),
                                      self.chosen)
        finally:
            painter.end()

    def mousePressEvent(self, event):
        at = QPointF(event.position())
        hit = next((key for key, box in self._boxes.items() if box.contains(at)), None)
        if hit is None:          # a near miss still counts: levels are thin
            near = [(abs(box.center().y() - at.y()), key) for key, box in self._boxes.items()
                    if box.left() - 6 <= at.x() <= box.right() + 6]
            hit = min(near)[1] if near and min(near)[0] < 10 else None
        if hit is not None:
            self.picked.emit(*hit)
