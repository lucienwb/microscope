"""The flat, ChemDraw-style view of the structure.

A read-only companion to the 3-D viewport: it shares the same camera object,
so the orientation carries over both ways — line up the molecule in 3-D,
switch here to see it as a Lewis structure, turn it a little more until the
ring is face-on and the substituents are clear, then save the picture or hand
it to ChemDraw. Nothing here edits the molecule; picking the angle is the
whole job.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter
from PySide6.QtWidgets import QWidget

from ..core import lewis
from ..core.molecule import Molecule
from ..render.camera import OrthoCamera
from ..render.lewis2d import LewisOptions, hidden_hydrogens
from ..render.lewisdraw import draw_lewis

PAGE_COLOR = QColor(255, 255, 255)
HINT_COLOR = QColor(140, 140, 140)
NOTE_COLOR = QColor(150, 110, 40)
SPIN_SPEED = 0.4        # degrees of in-plane rotation per pixel dragged


class LewisView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.camera = OrthoCamera()
        self.molecule: Molecule | None = None
        self.structure: lewis.LewisStructure | None = None
        self.options = LewisOptions()
        self._last_pos = None
        self.setAutoFillBackground(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(400, 300)

    # ------------------------------------------------------------------ data

    def set_camera(self, camera: OrthoCamera) -> None:
        """Share the 3-D viewport's camera so both views turn together."""
        self.camera = camera

    def set_molecule(self, molecule: Molecule | None) -> None:
        self.molecule = molecule
        self.structure = None if molecule is None else lewis.perceive(molecule)
        self.update()

    def set_option(self, name: str, value: bool) -> None:
        setattr(self.options, name, bool(value))
        self.update()

    def drawn_atoms(self) -> list[int]:
        """Indices of the atoms that appear in the drawing — the hydrogens
        folded into a label are left out, and out of an exported file too."""
        if self.structure is None:
            return []
        hidden = hidden_hydrogens(self.structure, self.options)
        return [a for a in range(self.structure.natoms) if not hidden[a]]

    def projected(self, width: int | None = None,
                  height: int | None = None) -> np.ndarray:
        """Where every atom sits on the page, in pixels."""
        if self.molecule is None:
            return np.zeros((0, 2))
        return self.camera.project(self.molecule.coords,
                                   self.width() if width is None else width,
                                   self.height() if height is None else height)

    # ------------------------------------------------------------------ paint

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), PAGE_COLOR)
        if self.structure is None:
            painter.setPen(HINT_COLOR)
            font = QFont()
            font.setPointSize(13)
            painter.setFont(font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "Open a file to draw its Lewis structure")
            painter.end()
            return
        plan = draw_lewis(painter, self.structure, self.camera,
                          self.width(), self.height(), self.options)
        self._draw_notes(painter, plan)
        painter.end()

    def _draw_notes(self, painter: QPainter, plan) -> None:
        """Anything the drawing cannot say for itself, in the corner.

        Overlapping atoms and a charge that disagrees with the file are both
        things the reader has to know about, and both are answered by turning
        the structure or distrusting the perception — so they belong on
        screen, and only on screen: exported figures never carry them.
        """
        if not plan.warnings:
            return
        font = QFont()
        font.setPointSize(10)
        painter.setFont(font)
        painter.setPen(NOTE_COLOR)
        line = painter.fontMetrics().height() + 2
        y = self.height() - 10 - line * (len(plan.warnings) - 1)
        for warning in plan.warnings:
            painter.drawText(QPointF(12, y), f"⚠ {warning}")
            y += line

    # ------------------------------------------------------------------ input

    def mousePressEvent(self, event):
        self._last_pos = event.position()

    def mouseMoveEvent(self, event):
        if self._last_pos is None or not event.buttons():
            return
        dx = event.position().x() - self._last_pos.x()
        dy = event.position().y() - self._last_pos.y()
        if event.buttons() & Qt.MouseButton.LeftButton:
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                self.camera.spin(-dx * SPIN_SPEED)   # turn the drawing in place
            else:
                self.camera.rotate_drag(dx, dy)
        elif event.buttons() & (Qt.MouseButton.RightButton
                                | Qt.MouseButton.MiddleButton):
            self.camera.pan_drag(dx, dy, self.height())
        self._last_pos = event.position()
        self.update()

    def mouseReleaseEvent(self, event):
        self._last_pos = None

    def wheelEvent(self, event):
        self.camera.zoom(event.angleDelta().y() / 120.0)
        self.update()
