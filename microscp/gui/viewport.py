"""Interactive 3D viewport: rotate/pan/zoom, picking, measurements, labels."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtOpenGLWidgets import QOpenGLWidget

from ..core import editing, geometry
from ..core.molecule import Molecule
from ..render.camera import OrthoCamera, orientation_along, orientation_from_plane
from ..render.glrenderer import MoleculeRenderer
from ..render.scene import build_scene
from ..render.styles import Style

_MARKER_COLOR = QColor(235, 130, 20)
_LINE_COLOR = QColor(60, 60, 60)
_PIN_COLOR = QColor(95, 95, 95)
_TEXT_COLOR = QColor(25, 25, 25)
_HALO_COLOR = QColor(255, 255, 255)

LABEL_MODES = ("none", "element", "element+number", "number")


class MoleculeViewport(QOpenGLWidget):
    selectionChanged = Signal(str)
    structureEdited = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.camera = OrthoCamera()
        self.style = Style()
        self.molecule: Molecule | None = None
        self.selection: list[int] = []
        self.pinned: list[list[int]] = []
        self.label_mode = "none"
        self._renderer = MoleculeRenderer()
        self._scene = None
        self._scene_dirty = False
        self._last_pos = None
        self._press_pos = None
        self._overlay_ok = True
        self._anim_timer: QTimer | None = None
        self._anim_base = None
        self._anim_disp = None
        self._anim_phase = 0.0
        self._undo_stack: list[tuple] = []
        self._redo_stack: list[tuple] = []
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(400, 300)

    # ------------------------------------------------------------------ data

    def set_molecule(self, molecule: Molecule | None, keep_camera: bool = False) -> None:
        self.stop_animation()
        old_natoms = self.molecule.natoms if self.molecule is not None else -1
        self.molecule = molecule
        self.selection.clear()
        if molecule is not None:
            if molecule.bonds is None:
                molecule.perceive_bonds()
            self._scene = build_scene(molecule, self.style)
            if molecule.natoms != old_natoms:
                self.pinned.clear()
            if not keep_camera:
                center, radius = molecule.bounding_sphere()
                self.camera.fit(center, radius)
        else:
            self._scene = None
            self.pinned.clear()
        self._scene_dirty = True
        self.selectionChanged.emit("")
        self.update()

    def reset_view(self) -> None:
        if self.molecule is not None:
            self.camera.reset_orientation()
            center, radius = self.molecule.bounding_sphere()
            self.camera.fit(center, radius)
            self.update()

    # ------------------------------------------------------------------ labels

    def set_label_mode(self, mode: str) -> None:
        if mode in LABEL_MODES:
            self.label_mode = mode
            self.update()

    def cycle_label_mode(self) -> str:
        idx = LABEL_MODES.index(self.label_mode)
        self.label_mode = LABEL_MODES[(idx + 1) % len(LABEL_MODES)]
        self.update()
        return self.label_mode

    def _atom_label(self, i: int) -> str:
        sym = self.molecule.symbols[i]
        if self.label_mode == "element":
            return sym
        if self.label_mode == "number":
            return str(i + 1)
        return f"{sym}{i + 1}"

    # ------------------------------------------------------------------ measurements

    @staticmethod
    def _measure_value(mol: Molecule, idxs: list[int]) -> str:
        pts = [mol.coords[i] for i in idxs]
        if len(idxs) == 2:
            return f"{geometry.distance(*pts):.3f} Å"
        if len(idxs) == 3:
            return f"{geometry.angle(*pts):.1f}°"
        return f"{geometry.dihedral(*pts):.1f}°"

    def measurement_text(self) -> str:
        if self.molecule is None or not self.selection:
            return ""
        if len(self.selection) == 1:
            i = self.selection[0]
            return f"{self.molecule.symbols[i]}{i + 1} selected"
        tags = "–".join(f"{self.molecule.symbols[i]}{i + 1}" for i in self.selection)
        kind = {2: "d", 3: "∠", 4: "φ"}[len(self.selection)]
        return f"{kind}({tags}) = {self._measure_value(self.molecule, self.selection)}"

    def pin_selection(self) -> bool:
        """Keep the current 2–4 atom measurement permanently displayed."""
        if self.molecule is None or not 2 <= len(self.selection) <= 4:
            return False
        idxs = list(self.selection)
        if idxs not in self.pinned:
            self.pinned.append(idxs)
        self.selection.clear()
        self.selectionChanged.emit("")
        self.update()
        return True

    def clear_pinned(self) -> None:
        self.pinned.clear()
        self.update()

    # ------------------------------------------------------------------ vibration animation

    @property
    def is_animating(self) -> bool:
        return self._anim_timer is not None and self._anim_timer.isActive()

    def animate_mode(self, displacements, amplitude: float = 0.35) -> bool:
        """Animate a normal mode: coords oscillate along *displacements*."""
        self.stop_animation()
        if self.molecule is None or displacements is None:
            return False
        d = np.asarray(displacements, dtype=float)
        if d.shape != self.molecule.coords.shape:
            return False
        peak = np.abs(d).max()
        if peak < 1e-9:
            return False
        self._anim_base = self.molecule.coords.copy()
        self._anim_disp = d / peak * amplitude
        self._anim_phase = 0.0
        if self._anim_timer is None:
            self._anim_timer = QTimer(self)
            self._anim_timer.timeout.connect(self._anim_tick)
        self._anim_timer.start(33)
        return True

    def stop_animation(self) -> None:
        if self.is_animating:
            self._anim_timer.stop()
            if self.molecule is not None and self._anim_base is not None:
                self.molecule.coords = self._anim_base
                self._rebuild_scene()
        self._anim_base = None
        self._anim_disp = None

    def _anim_tick(self):
        self._anim_phase += 0.3
        self.molecule.coords = self._anim_base + np.sin(self._anim_phase) * self._anim_disp
        self._rebuild_scene()

    def _rebuild_scene(self):
        self._scene = build_scene(self.molecule, self.style)
        self._scene_dirty = True
        self.update()

    # ------------------------------------------------------------------ editing

    def snapshot(self) -> tuple:
        mol = self.molecule
        return (mol.coords.copy(), list(mol.symbols),
                None if mol.bonds is None else mol.bonds.copy())

    def apply_snapshot(self, snap: tuple) -> None:
        coords, symbols, bonds = snap
        mol = self.molecule
        if mol is not None and len(symbols) == mol.natoms:
            mol.coords = coords.copy()
            mol.symbols = list(symbols)
            mol.bonds = None if bonds is None else bonds.copy()
            self._rebuild_scene()
            self.selectionChanged.emit(self.measurement_text())
        else:
            new = Molecule(list(symbols), coords.copy(),
                           charge=mol.charge if mol else 0,
                           multiplicity=mol.multiplicity if mol else 1,
                           title=mol.title if mol else "")
            new.bonds = None if bonds is None else bonds.copy()
            self.set_molecule(new, keep_camera=True)

    def clear_history(self) -> None:
        self._undo_stack.clear()
        self._redo_stack.clear()

    def _push_undo(self, snap: tuple | None = None) -> None:
        self._undo_stack.append(snap if snap is not None else self.snapshot())
        del self._undo_stack[:-100]
        self._redo_stack.clear()
        self.structureEdited.emit()

    def undo(self) -> bool:
        if not self._undo_stack or self.molecule is None:
            return False
        self.stop_animation()
        self._redo_stack.append(self.snapshot())
        self.apply_snapshot(self._undo_stack.pop())
        return True

    def redo(self) -> bool:
        if not self._redo_stack or self.molecule is None:
            return False
        self.stop_animation()
        self._undo_stack.append(self.snapshot())
        self.apply_snapshot(self._redo_stack.pop())
        return True

    def preview_adjust(self, base_coords, value: float,
                       move_fragment: bool = True) -> bool:
        """Recompute the selected internal coordinate from *base_coords*."""
        mol = self.molecule
        sel = list(self.selection)
        if mol is None or not 2 <= len(sel) <= 4:
            return False
        mol.coords = base_coords.copy()
        try:
            if len(sel) == 2:
                mol.coords = editing.set_distance(mol, *sel, value, move_fragment)
            elif len(sel) == 3:
                mol.coords = editing.set_angle(mol, *sel, value, move_fragment)
            else:
                mol.coords = editing.set_dihedral(mol, *sel, value, move_fragment)
        except editing.EditError:
            self._rebuild_scene()
            return False
        self._rebuild_scene()
        self.selectionChanged.emit(self.measurement_text())
        return True

    def commit_adjust(self, base_snapshot: tuple) -> None:
        """Record the pre-edit state so the applied preview becomes undoable."""
        self._push_undo(base_snapshot)

    def cancel_adjust(self, base_coords) -> None:
        if self.molecule is not None:
            self.molecule.coords = base_coords.copy()
            self._rebuild_scene()
            self.selectionChanged.emit(self.measurement_text())

    def delete_selected(self) -> bool:
        if self.molecule is None or not self.selection:
            return False
        self.stop_animation()
        snap = self.snapshot()
        try:
            new = editing.delete_atoms(self.molecule, self.selection)
        except editing.EditError:
            return False
        self._push_undo(snap)
        self.set_molecule(new, keep_camera=True)
        return True

    def recompute_bonds(self) -> None:
        if self.molecule is not None:
            self.molecule.perceive_bonds()
            self._rebuild_scene()

    # ------------------------------------------------------------------ view helpers

    def center_on_selection(self) -> bool:
        """Make the last-selected atom the rotation/pan center (GaussView-style)."""
        if self.molecule is None or not self.selection:
            return False
        self.camera.center = self.molecule.coords[self.selection[-1]].copy()
        self.update()
        return True

    def center_on_molecule(self) -> None:
        if self.molecule is None:
            return
        center, _ = self.molecule.bounding_sphere()
        self.camera.center = center
        self.update()

    def align_to_selection(self) -> str | None:
        """2 atoms: look along the bond. 3 atoms: put their plane in the screen."""
        if self.molecule is None:
            return None
        coords = self.molecule.coords
        sel = self.selection
        if len(sel) == 2:
            rot = orientation_along(coords[sel[1]] - coords[sel[0]])
            if rot is None:
                return None
            self.camera.rotation = rot
            self.update()
            return "Looking along the selected bond"
        if len(sel) == 3:
            rot = orientation_from_plane(*(coords[i] for i in sel))
            if rot is None:
                return None
            self.camera.rotation = rot
            self.update()
            return "Selected plane aligned to the screen"
        return None

    # ------------------------------------------------------------------ GL

    def initializeGL(self):
        self._renderer.initialize()
        self._scene_dirty = True

    def paintGL(self):
        dpr = self.devicePixelRatioF()
        w = max(int(self.width() * dpr), 1)
        h = max(int(self.height() * dpr), 1)
        if self._scene_dirty:
            if self._scene is not None:
                self._renderer.set_scene(self._scene)
            self._scene_dirty = False
        if self._scene is None:
            from OpenGL import GL
            GL.glClearColor(*self.style.background, 1.0)
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
            return
        view = self.camera.view_matrix()
        proj = self.camera.proj_matrix(w / h)
        self._renderer.draw(view, proj, w, h, background=(*self.style.background, 1.0))
        overlay_needed = bool(self.selection or self.pinned or self.label_mode != "none")
        if overlay_needed and self._overlay_ok:
            try:
                self._draw_overlay()
            except Exception:
                self._overlay_ok = False

    # ------------------------------------------------------------------ overlay

    @staticmethod
    def _draw_halo_text(painter: QPainter, x: float, y: float, text: str) -> None:
        x, y = int(x), int(y)
        painter.setPen(QPen(_HALO_COLOR))
        for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            painter.drawText(x + dx, y + dy, text)
        painter.setPen(QPen(_TEXT_COLOR))
        painter.drawText(x, y, text)

    def _overlay_metrics(self) -> dict:
        """Sizes that track the zoom level so text stays legible at any scale."""
        px_per_ang = self.height() / (2.0 * self.camera.half_height)
        return {
            "label_px": int(np.clip(0.17 * px_per_ang, 9, 34)),
            "value_px": int(np.clip(0.20 * px_per_ang, 10, 40)),
            "marker_r": float(np.clip(0.14 * px_per_ang, 8, 60)),
            "offset": int(np.clip(0.09 * px_per_ang, 5, 24)),
            "line_w": float(np.clip(0.018 * px_per_ang, 1.3, 3.5)),
        }

    def _draw_overlay(self):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        mol = self.molecule
        pts = self.camera.project(mol.coords, self.width(), self.height())
        m = self._overlay_metrics()

        if self.label_mode != "none":
            font = QFont()
            font.setPixelSize(m["label_px"])
            painter.setFont(font)
            for i in range(mol.natoms):
                self._draw_halo_text(painter, pts[i, 0] + m["offset"],
                                     pts[i, 1] - m["offset"], self._atom_label(i))

        if self.pinned:
            font = QFont()
            font.setPixelSize(m["value_px"])
            painter.setFont(font)
            for idxs in self.pinned:
                p = pts[idxs]
                pen = QPen(_PIN_COLOR)
                pen.setWidthF(m["line_w"])
                pen.setStyle(Qt.PenStyle.DashLine)
                painter.setPen(pen)
                for k in range(len(p) - 1):
                    painter.drawLine(int(p[k, 0]), int(p[k, 1]),
                                     int(p[k + 1, 0]), int(p[k + 1, 1]))
                if len(idxs) == 2:
                    ax, ay = (p[0] + p[1]) / 2.0
                elif len(idxs) == 3:
                    ax, ay = p[1]
                else:
                    ax, ay = (p[1] + p[2]) / 2.0
                self._draw_halo_text(painter, ax + m["offset"], ay - m["offset"],
                                     self._measure_value(mol, idxs))

        if self.selection:
            p = pts[self.selection]
            pen = QPen(_LINE_COLOR)
            pen.setWidthF(m["line_w"])
            pen.setStyle(Qt.PenStyle.DashLine)
            painter.setPen(pen)
            for k in range(len(p) - 1):
                painter.drawLine(int(p[k, 0]), int(p[k, 1]),
                                 int(p[k + 1, 0]), int(p[k + 1, 1]))
            marker_pen = QPen(_MARKER_COLOR)
            marker_pen.setWidthF(max(2.0, m["line_w"]))
            painter.setPen(marker_pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            r = m["marker_r"]
            for x, y in p:
                painter.drawEllipse(int(x - r), int(y - r), int(2 * r), int(2 * r))
            text = self.measurement_text()
            if text:
                font = QFont()
                font.setPointSize(11)   # HUD text: fixed, screen-anchored
                painter.setFont(font)
                self._draw_halo_text(painter, 12, self.height() - 14, text)
        painter.end()

    # ------------------------------------------------------------------ input

    def mousePressEvent(self, event):
        self._press_pos = event.position()
        self._last_pos = event.position()

    def mouseMoveEvent(self, event):
        if self._last_pos is None:
            return
        pos = event.position()
        dx = pos.x() - self._last_pos.x()
        dy = pos.y() - self._last_pos.y()
        if event.buttons() & Qt.MouseButton.LeftButton:
            self.camera.rotate_drag(dx, dy)
        elif event.buttons() & (Qt.MouseButton.RightButton | Qt.MouseButton.MiddleButton):
            self.camera.pan_drag(dx, dy, self.height())
        self._last_pos = pos
        self.update()

    def mouseReleaseEvent(self, event):
        if (self._press_pos is not None
                and event.button() == Qt.MouseButton.LeftButton
                and (event.position() - self._press_pos).manhattanLength() < 4
                and self.molecule is not None):
            self._pick(event.position())
        self._press_pos = None
        self._last_pos = None

    def wheelEvent(self, event):
        self.camera.zoom(event.angleDelta().y() / 120.0)
        self.update()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Escape:
            self.selection.clear()
            self.selectionChanged.emit("")
            self.update()
        else:
            super().keyPressEvent(event)

    def _pick(self, pos):
        dpr = self.devicePixelRatioF()
        w = max(int(self.width() * dpr), 1)
        h = max(int(self.height() * dpr), 1)
        x = pos.x() * dpr
        y = h - 1 - pos.y() * dpr
        self.makeCurrent()
        try:
            idx = self._renderer.pick_atom(
                x, y, self.camera.view_matrix(), self.camera.proj_matrix(w / h), w, h)
        finally:
            self.doneCurrent()
        if idx is not None and 0 <= idx < self.molecule.natoms:
            if idx in self.selection:
                self.selection.remove(idx)
            else:
                self.selection.append(idx)
                if len(self.selection) > 4:
                    self.selection.pop(0)
        else:
            self.selection.clear()
        self.selectionChanged.emit(self.measurement_text())
        self.update()
