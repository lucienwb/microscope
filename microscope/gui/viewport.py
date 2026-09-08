"""Interactive 3D viewport: rotate/pan/zoom, picking, measurements, labels."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter
from PySide6.QtOpenGLWidgets import QOpenGLWidget

from ..core import editing, measure
from ..core.history import EditHistory, Snapshot
from ..core.molecule import Molecule
from ..core.volume import VolumeData
from ..render.camera import OrthoCamera, orientation_along, orientation_from_plane
from ..render.glrenderer import MoleculeRenderer
from ..render.scene import build_scene, build_surface_meshes
from ..render.styles import Style, make_style
from . import annotations, gizmo

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
        self.show_gizmo = True          # manipulator on the current selection
        self.show_axes = False          # world orientation triad in the corner
        self._gizmo_hover = None
        self._gizmo_drag = None         # handle being dragged
        self._drag_base = None          # coords when the drag started
        self._drag_snapshot = None
        self._drag_origin = None        # gizmo centre at drag start (world)
        self._drag_start_px = None
        self._drag_angle = 0.0          # accumulated rotation, degrees
        self._drag_last_angle = 0.0     # last screen angle, radians
        self._drag_moved = False
        self._anim_timer: QTimer | None = None
        self._anim_base = None
        self._anim_disp = None
        self._anim_phase = 0.0
        self._history = EditHistory()
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)     # so gizmo handles light up on hover
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
        return measure.describe(self.molecule, self.selection)

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

    def snapshot(self) -> Snapshot:
        return Snapshot.of(self.molecule)

    def apply_snapshot(self, snap: Snapshot) -> None:
        """Put a recorded state back, redrawing whatever it changed."""
        mol = self.molecule
        if snap.fits(mol):
            snap.restore_into(mol)
            self._rebuild_scene()
            self.selectionChanged.emit(self.measurement_text())
            return
        # a different number of atoms: it has to become a new molecule
        new = Molecule(list(snap.symbols), snap.coords.copy(),
                       charge=mol.charge if mol else 0,
                       multiplicity=mol.multiplicity if mol else 1,
                       title=mol.title if mol else "")
        new.bonds = None if snap.bonds is None else snap.bonds.copy()
        self.set_molecule(new, keep_camera=True)

    def clear_history(self) -> None:
        self._history.clear()

    def _push_undo(self, snap: Snapshot | None = None) -> None:
        self._history.push(snap if snap is not None else self.snapshot())
        self.structureEdited.emit()

    def undo(self) -> bool:
        return self._step(self._history.undo)

    def redo(self) -> bool:
        return self._step(self._history.redo)

    def _step(self, move) -> bool:
        """Undo and redo differ only in which way they walk the history."""
        if self.molecule is None:
            return False
        snap = move(self.snapshot())
        if snap is None:
            return False
        self.stop_animation()
        self.apply_snapshot(snap)
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

    def select_fragment(self) -> int:
        """Grow the selection to every atom connected to it (a whole group).

        The manipulator moves exactly what is selected, so this is how a
        substituent or a ligand gets picked up in one click.
        """
        if self.molecule is None or not self.selection:
            return 0
        if self.molecule.bonds is None:
            self.molecule.perceive_bonds()
        self.selection = editing.connected_fragment(self.molecule.bonds,
                                                    self.selection)
        self.selectionChanged.emit(self.measurement_text())
        self.update()
        return len(self.selection)

    # ------------------------------------------------------------------ manipulator

    @property
    def gizmo_visible(self) -> bool:
        return bool(self.show_gizmo and self.molecule is not None and self.selection)

    def set_show_gizmo(self, on: bool) -> None:
        self.show_gizmo = bool(on)
        self._gizmo_hover = None
        self.update()

    def set_show_axes(self, on: bool) -> None:
        self.show_axes = bool(on)
        self.update()

    def _gizmo_center(self):
        if not self.gizmo_visible:
            return None
        return gizmo.selection_center(self.molecule.coords, self.selection)

    def _gizmo_hit(self, pos):
        center = self._gizmo_center()
        if center is None:
            return None
        return gizmo.hit_test(self.camera, self.width(), self.height(), center,
                              (pos.x(), pos.y()))

    def _begin_gizmo_drag(self, handle, pos) -> None:
        self.stop_animation()
        self._gizmo_drag = handle
        self._drag_snapshot = self.snapshot()
        self._drag_base = self.molecule.coords.copy()
        self._drag_origin = self._gizmo_center()
        self._drag_start_px = np.array([pos.x(), pos.y()], dtype=float)
        self._drag_angle = 0.0
        self._drag_moved = False
        if handle[0] == "ring":
            self._drag_last_angle = gizmo.screen_angle(
                self.camera, self.width(), self.height(), self._drag_origin,
                self._drag_start_px)

    def _update_gizmo_drag(self, pos, snap: bool) -> None:
        kind, k = self._gizmo_drag
        px = np.array([pos.x(), pos.y()], dtype=float)
        w, h = self.width(), self.height()
        mol = self.molecule
        mol.coords = self._drag_base            # transforms read from the base
        try:
            if kind == "ring":
                ang = gizmo.screen_angle(self.camera, w, h, self._drag_origin, px)
                step = gizmo.wrap_angle(ang - self._drag_last_angle)
                self._drag_last_angle = ang
                self._drag_angle += np.degrees(step) * gizmo.rotation_sign(
                    self.camera, k)
                value = self._drag_angle
                if snap:
                    value = round(value / 15.0) * 15.0
                coords = editing.rotate_atoms(mol, self.selection,
                                              gizmo.AXIS_VECTORS[k], value,
                                              pivot=self._drag_origin)
                self._drag_moved = abs(value) > 1e-9
                message = f"rotate {gizmo.AXIS_NAMES[k]} {value:+.1f}°"
            else:
                if kind == "center":
                    delta = gizmo.plane_translation(self.camera, h,
                                                    px - self._drag_start_px)
                else:
                    t = gizmo.axis_translation(self.camera, w, h,
                                               self._drag_origin, k,
                                               px - self._drag_start_px)
                    if snap:
                        t = round(t / 0.1) * 0.1
                    delta = gizmo.AXIS_VECTORS[k] * t
                coords = editing.translate_atoms(mol, self.selection, delta)
                self._drag_moved = bool(np.linalg.norm(delta) > 1e-9)
                message = ("move " + ("in view plane" if kind == "center"
                                      else gizmo.AXIS_NAMES[k])
                           + f" {np.linalg.norm(delta):.3f} Å")
        except editing.EditError:
            mol.coords = self._drag_base.copy()
            return
        mol.coords = coords
        self._rebuild_scene()
        self.selectionChanged.emit(
            f"{len(self.selection)} atom(s): {message}"
            if len(self.selection) > 4 else
            f"{self.measurement_text()}  ·  {message}")

    def _end_gizmo_drag(self) -> None:
        if self._drag_moved and self._drag_snapshot is not None:
            self._push_undo(self._drag_snapshot)
        self._gizmo_drag = None
        self._drag_base = None
        self._drag_snapshot = None
        self._drag_origin = None
        self._drag_start_px = None
        self._drag_moved = False
        self.selectionChanged.emit(self.measurement_text())
        self.update()

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
        overlay_needed = bool(self.selection or self.pinned or self.show_axes
                              or self.label_mode != "none")
        if overlay_needed and self._overlay_ok:
            try:
                self._draw_overlay()
            except Exception:
                self._overlay_ok = False

    # ------------------------------------------------------------------ overlay

    def _draw_overlay(self):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        annotations.draw_overlay(
            painter, self.molecule, self.camera, self.width(), self.height(),
            label_mode=self.label_mode, pinned=self.pinned,
            selection=self.selection, selection_text=self.measurement_text(),
            show_axes=self.show_axes,
            gizmo_center=self._gizmo_center() if self.gizmo_visible else None,
            gizmo_hover=self._gizmo_hover, gizmo_active=self._gizmo_drag)
        painter.end()

    # ------------------------------------------------------------------ input

    def mousePressEvent(self, event):
        self._press_pos = event.position()
        self._last_pos = event.position()
        if event.button() == Qt.MouseButton.LeftButton and self.gizmo_visible:
            handle = self._gizmo_hit(event.position())
            if handle is not None:
                self._begin_gizmo_drag(handle, event.position())
                self.update()

    def mouseMoveEvent(self, event):
        pos = event.position()
        if self._gizmo_drag is not None:
            self._update_gizmo_drag(
                pos, bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier))
            self._last_pos = pos
            self.update()
            return
        if not event.buttons():                     # hover: light up a handle
            hover = self._gizmo_hit(pos) if self.gizmo_visible else None
            if hover != self._gizmo_hover:
                self._gizmo_hover = hover
                self.update()
            return
        if self._last_pos is None:
            return
        dx = pos.x() - self._last_pos.x()
        dy = pos.y() - self._last_pos.y()
        if event.buttons() & Qt.MouseButton.LeftButton:
            self.camera.rotate_drag(dx, dy)
        elif event.buttons() & (Qt.MouseButton.RightButton | Qt.MouseButton.MiddleButton):
            self.camera.pan_drag(dx, dy, self.height())
        self._last_pos = pos
        self.update()

    def mouseReleaseEvent(self, event):
        if self._gizmo_drag is not None:
            self._end_gizmo_drag()
        elif (self._press_pos is not None
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
