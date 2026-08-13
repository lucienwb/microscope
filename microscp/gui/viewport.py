"""Interactive 3D viewport: rotate/pan/zoom, picking, measurements, labels."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtOpenGLWidgets import QOpenGLWidget

from ..core import editing
from ..core.molecule import Molecule
from ..core.volume import VolumeData
from ..render.camera import OrthoCamera, orientation_along, orientation_from_plane
from ..render.glrenderer import MoleculeRenderer
from ..render.scene import build_scene, build_surface_meshes
from ..render.styles import Style, make_style
from . import annotations

_MARKER_COLOR = QColor(235, 130, 20)
_LINE_COLOR = QColor(60, 60, 60)
_PIN_COLOR = QColor(95, 95, 95)

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
        self.atom_reps: np.ndarray | None = None   # per-atom REP_* codes
        self.label_mode = "none"
        self._renderer = MoleculeRenderer()
        self._scene = None
        self._scene_dirty = False
        self.volume: VolumeData | None = None
        self.isovalue: float = 0.02
        self.surface_visible = True
        self._surface_meshes: list = []
        self._mesh_dirty = False
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
        if self.volume is not None:      # a new structure invalidates the grid
            self.volume = None
            self._surface_meshes = []
            self._mesh_dirty = True
        old_natoms = self.molecule.natoms if self.molecule is not None else -1
        self.molecule = molecule
        self.selection.clear()
        if molecule is not None:
            if molecule.bonds is None:
                molecule.perceive_bonds()
            if molecule.natoms != old_natoms:
                self.pinned.clear()
                self.atom_reps = None    # mixed representations reset too
            self._scene = build_scene(molecule, self.style, self.atom_reps)
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

    # ------------------------------------------------------------------ volume

    @property
    def has_volume(self) -> bool:
        return self.volume is not None

    def set_volume(self, volume: VolumeData | None, isovalue: float | None = None) -> None:
        """Attach volumetric data (cube grid) and show its isosurface."""
        self.volume = volume
        if volume is not None:
            self.isovalue = float(isovalue) if isovalue else volume.suggest_isovalue()
            self.surface_visible = True
        self._update_surface()

    def set_isosurface(self, isovalue: float | None = None,
                       visible: bool | None = None,
                       opacity: float | None = None,
                       positive_color: tuple | None = None,
                       negative_color: tuple | None = None) -> None:
        if isovalue is not None:
            self.isovalue = float(isovalue)
        if visible is not None:
            self.surface_visible = bool(visible)
        if opacity is not None:
            self.style.surface_opacity = float(opacity)
        if positive_color is not None:
            self.style.surface_positive = tuple(positive_color)
        if negative_color is not None:
            self.style.surface_negative = tuple(negative_color)
        self._update_surface()

    def _update_surface(self) -> None:
        if self.volume is not None and self.surface_visible and self.isovalue:
            self._surface_meshes = build_surface_meshes(
                self.volume, self.isovalue, self.style)
        else:
            self._surface_meshes = []
        self._mesh_dirty = True
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

    # ------------------------------------------------------------------ measurements

    def measurement_text(self) -> str:
        if self.molecule is None or not self.selection:
            return ""
        if len(self.selection) == 1:
            i = self.selection[0]
            return f"{self.molecule.symbols[i]}{i + 1} selected"
        if len(self.selection) > 4:
            return f"{len(self.selection)} atoms selected"
        tags = "–".join(f"{self.molecule.symbols[i]}{i + 1}" for i in self.selection)
        kind = {2: "d", 3: "∠", 4: "φ"}[len(self.selection)]
        value = annotations.measurement_value(self.molecule, self.selection)
        return f"{kind}({tags}) = {value}"

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
        if self.atom_reps is not None and len(self.atom_reps) != self.molecule.natoms:
            self.atom_reps = None        # atom count changed (delete/undo)
        self._scene = build_scene(self.molecule, self.style, self.atom_reps)
        self._scene_dirty = True
        self.update()

    def set_atom_representation(self, rep: int, atoms: list[int] | None = None) -> None:
        """Apply a REP_* code to *atoms* (None/empty = every atom)."""
        if self.molecule is None:
            return
        if self.atom_reps is None or len(self.atom_reps) != self.molecule.natoms:
            self.atom_reps = np.zeros(self.molecule.natoms, dtype=int)
        if atoms:
            self.atom_reps[list(atoms)] = rep
        else:
            self.atom_reps[:] = rep
        self._rebuild_scene()

    def set_representation(self, name: str) -> None:
        """Switch style preset (cylview / houk), keeping user adjustments."""
        if name == self.style.name:
            return
        new = make_style(name)
        new.show_hbonds = self.style.show_hbonds
        new.surface_positive = self.style.surface_positive
        new.surface_negative = self.style.surface_negative
        new.surface_opacity = self.style.surface_opacity
        self.style = new
        if self.molecule is not None:
            self._rebuild_scene()
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
        if self._mesh_dirty:
            self._renderer.set_meshes(self._surface_meshes)
            self._mesh_dirty = False
        if self._scene is None:
            from OpenGL import GL
            GL.glClearColor(*self.style.background, 1.0)
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
            return
        view = self.camera.view_matrix()
        proj = self.camera.proj_matrix(w / h)
        self._renderer.set_style_params(self.style.quadrant_color,
                                        self.style.quadrant_width)
        self._renderer.draw(view, proj, w, h, background=(*self.style.background, 1.0))
        overlay_needed = bool(self.selection or self.pinned or self.label_mode != "none")
        if overlay_needed and self._overlay_ok:
            try:
                self._draw_overlay()
            except Exception:
                self._overlay_ok = False

    # ------------------------------------------------------------------ overlay

    def _draw_overlay(self):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        mol = self.molecule
        w, h = self.width(), self.height()
        pts = self.camera.project(mol.coords, w, h)
        m = annotations.overlay_metrics(h, self.camera.half_height)

        if self.label_mode != "none":
            annotations.draw_atom_labels(painter, mol, pts, m, self.label_mode)

        for idxs in self.pinned:
            annotations.draw_measurement(painter, mol, self.camera, w, h,
                                         idxs, m, color=_PIN_COLOR)

        if self.selection:
            if 2 <= len(self.selection) <= 4:   # larger selections: markers only
                annotations.draw_measurement(painter, mol, self.camera, w, h,
                                             self.selection, m, color=_LINE_COLOR)
            marker_pen = QPen(_MARKER_COLOR)
            marker_pen.setWidthF(max(2.0, m["line_w"]))
            painter.setPen(marker_pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            r = m["marker_r"]
            for x, y in pts[self.selection]:
                painter.drawEllipse(int(x - r), int(y - r), int(2 * r), int(2 * r))
            text = self.measurement_text()
            if text:
                font = QFont()
                font.setPointSize(11)   # HUD text: fixed, screen-anchored
                painter.setFont(font)
                annotations.draw_halo_text(painter, 12, h - 14, text)
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
                self.selection.append(idx)   # unlimited; info shown for 2-4
        else:
            self.selection.clear()
        self.selectionChanged.emit(self.measurement_text())
        self.update()
