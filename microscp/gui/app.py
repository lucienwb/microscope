"""Main application window."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import (
    QAction, QActionGroup, QColor, QImage, QImageWriter, QKeySequence, QPainter,
    QSurfaceFormat,
)
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QColorDialog, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QMainWindow,
    QMessageBox, QPushButton, QSlider, QSpinBox, QToolBar, QWidget,
)

from .. import __version__
from .. import io as mio
from ..core import geometry
from ..core.results import ParseResult
from ..render.offscreen import render_molecule_image
from ..render.scene import REP_BALL, REP_LINE, REP_STICK
from .annotations import draw_annotations
from .spectra import SpectraDock
from .viewport import MoleculeViewport

OPEN_FILTER = (
    "Molecular files (*.xyz *.log *.out *.fchk *.fck *.fch *.gjf *.com *.gau "
    "*.pdb *.molden *.cube *.cub);;All files (*)"
)
SAVE_FILTER = ("XYZ (*.xyz);;Gaussian input (*.gjf *.com);;"
               "ORCA input (*.inp);;Q-Chem input (*.in *.qcin);;PDB (*.pdb)")
_FILTER_DEFAULT_EXT = {"XYZ": ".xyz", "Gaussian": ".gjf", "ORCA": ".inp",
                       "Q-Chem": ".in", "PDB": ".pdb"}

EXPORT_FILTER = "PNG image (*.png);;TIFF image, uncompressed (*.tif *.tiff)"
_EXPORT_DEFAULT_EXT = {"PNG": ".png", "TIFF": ".tif"}

# Curated isosurface color pairs (+ lobe, - lobe); custom colors via the picker.
SURFACE_PALETTES = (
    ("Blue / Red", (0.29, 0.44, 0.86), (0.88, 0.38, 0.22)),
    ("Red / Blue", (0.88, 0.38, 0.22), (0.29, 0.44, 0.86)),
    ("Teal / Orange", (0.13, 0.59, 0.62), (0.95, 0.52, 0.16)),
    ("Purple / Gold", (0.55, 0.36, 0.76), (0.93, 0.69, 0.13)),
    ("Green / Magenta", (0.30, 0.63, 0.36), (0.79, 0.29, 0.62)),
    ("Slate / Silver", (0.36, 0.42, 0.52), (0.72, 0.75, 0.80)),
)


def _write_image(image: QImage, path: str) -> None:
    """Save as PNG or TIFF — both keep the transparent background intact."""
    writer = QImageWriter(path)
    if Path(path).suffix.lower() in (".tif", ".tiff"):
        writer.setCompression(0)                # lossless master copy
    if not writer.write(image):
        raise RuntimeError(writer.errorString() or f"could not write {path}")


class ExportImageDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export Image")
        layout = QFormLayout(self)
        self.scale = QSpinBox()
        self.scale.setRange(1, 8)
        self.scale.setValue(3)
        self.transparent = QCheckBox()
        self.transparent.setChecked(True)
        self.annotations = QCheckBox()
        self.annotations.setChecked(True)
        self.annotations.setToolTip(
            "Draw atom labels and pinned measurements into the exported image")
        layout.addRow("Resolution scale:", self.scale)
        layout.addRow("Transparent background:", self.transparent)
        layout.addRow("Labels and measurements:", self.annotations)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)


class AdjustDialog(QDialog):
    """Adjust the selected distance/angle/dihedral with live preview."""

    def __init__(self, viewport: MoleculeViewport, parent=None):
        super().__init__(parent)
        self.viewport = viewport
        viewport.stop_animation()
        mol = viewport.molecule
        sel = list(viewport.selection)
        self._base_snapshot = viewport.snapshot()
        self._base_coords = mol.coords.copy()
        self._accepted = False

        pts = [mol.coords[x] for x in sel]
        tags = "–".join(f"{mol.symbols[x]}{x + 1}" for x in sel)
        self.value = QDoubleSpinBox()
        if len(sel) == 2:
            kind = "distance"
            current = geometry.distance(*pts)
            self.value.setRange(0.3, 99.0)
            self.value.setDecimals(3)
            self.value.setSingleStep(0.01)
            self.value.setSuffix(" Å")
        elif len(sel) == 3:
            kind = "angle"
            current = geometry.angle(*pts)
            self.value.setRange(1.0, 179.9)
            self.value.setDecimals(2)
            self.value.setSingleStep(1.0)
            self.value.setSuffix("°")
        else:
            kind = "dihedral"
            current = geometry.dihedral(*pts)
            self.value.setRange(-180.0, 180.0)
            self.value.setDecimals(2)
            self.value.setSingleStep(5.0)
            self.value.setWrapping(True)
            self.value.setSuffix("°")
        self.value.setValue(current)
        self.setWindowTitle(f"Adjust {kind}")

        self.fragment = QCheckBox("Move attached fragment")
        self.fragment.setChecked(True)
        self.fragment.setToolTip(
            "Move every atom bonded beyond the adjusted coordinate.\n"
            "For ring bonds only the end atom moves regardless.")

        layout = QFormLayout(self)
        layout.addRow(QLabel(f"{kind} {tags}"))
        layout.addRow("New value:", self.value)
        layout.addRow(self.fragment)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                   | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

        self.value.valueChanged.connect(self._preview)
        self.fragment.toggled.connect(self._preview)

    def _preview(self):
        self.viewport.preview_adjust(self._base_coords, self.value.value(),
                                     self.fragment.isChecked())

    def accept(self):
        self._preview()
        self.viewport.commit_adjust(self._base_snapshot)
        self._accepted = True
        super().accept()

    def reject(self):
        self.viewport.cancel_adjust(self._base_coords)
        super().reject()


class SurfaceDialog(QDialog):
    """Modeless isosurface controls for the loaded cube file (live update)."""

    def __init__(self, viewport: MoleculeViewport, volumes: list, parent=None):
        super().__init__(parent)
        self.viewport = viewport
        self.volumes = volumes
        self.setWindowTitle("Isosurface")
        layout = QFormLayout(self)

        if len(volumes) > 1:
            self.which = QComboBox()
            self.which.addItems([v.label or f"grid {i + 1}"
                                 for i, v in enumerate(volumes)])
            self.which.setCurrentIndex(
                next((i for i, v in enumerate(volumes) if v is viewport.volume), 0))
            self.which.currentIndexChanged.connect(self._switch_volume)
            layout.addRow("Data:", self.which)
        elif volumes and volumes[0].label:
            layout.addRow(QLabel(volumes[0].label))

        self.iso = QDoubleSpinBox()
        self.iso.setRange(0.0001, 10.0)
        self.iso.setDecimals(4)
        self.iso.setSingleStep(0.005)
        self.iso.setValue(viewport.isovalue)
        self.iso.valueChanged.connect(lambda v: viewport.set_isosurface(isovalue=v))
        layout.addRow("Isovalue (±):", self.iso)

        self.opacity = QSlider(Qt.Orientation.Horizontal)
        self.opacity.setRange(10, 100)
        self.opacity.setValue(int(viewport.style.surface_opacity * 100))
        self.opacity.valueChanged.connect(
            lambda v: viewport.set_isosurface(opacity=v / 100.0))
        layout.addRow("Opacity:", self.opacity)

        self.palette_box = QComboBox()
        for name, _pos, _neg in SURFACE_PALETTES:
            self.palette_box.addItem(name)
        current = (viewport.style.surface_positive, viewport.style.surface_negative)
        match = next((i for i, (_, pos, neg) in enumerate(SURFACE_PALETTES)
                      if self._rgb_close(pos, current[0])
                      and self._rgb_close(neg, current[1])), None)
        if match is None:
            self._select_custom_entry()
        else:
            self.palette_box.setCurrentIndex(match)
        self.palette_box.currentIndexChanged.connect(self._palette_chosen)
        layout.addRow("Colors:", self.palette_box)

        self.pos_btn = QPushButton()
        self.neg_btn = QPushButton()
        for btn, tip in ((self.pos_btn, "+ lobe"), (self.neg_btn, "− lobe")):
            btn.setFixedSize(56, 22)
            btn.setToolTip(f"{tip} color — click to pick any color (RGB or HEX)")
        self.pos_btn.clicked.connect(lambda: self._pick_color(positive=True))
        self.neg_btn.clicked.connect(lambda: self._pick_color(positive=False))
        swatches = QWidget()
        row = QHBoxLayout(swatches)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.pos_btn)
        row.addWidget(self.neg_btn)
        row.addStretch()
        layout.addRow("Lobes (+ / −):", swatches)
        self._refresh_swatches()

        self.visible = QCheckBox("Show surface")
        self.visible.setChecked(viewport.surface_visible)
        self.visible.toggled.connect(lambda on: viewport.set_isosurface(visible=on))
        layout.addRow(self.visible)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        buttons.clicked.connect(self.close)
        layout.addRow(buttons)

    def _switch_volume(self, index: int):
        self.viewport.set_volume(self.volumes[index], isovalue=self.iso.value())

    # ------------------------------------------------------------------ colors

    @staticmethod
    def _rgb_close(a: tuple, b: tuple) -> bool:
        return max(abs(x - y) for x, y in zip(a, b)) < 1e-3

    def _refresh_swatches(self):
        for btn, rgb in ((self.pos_btn, self.viewport.style.surface_positive),
                         (self.neg_btn, self.viewport.style.surface_negative)):
            btn.setStyleSheet(
                f"background-color: {QColor.fromRgbF(*rgb).name()}; "
                "border: 1px solid #888; border-radius: 3px;")

    def _palette_chosen(self, index: int):
        if index >= len(SURFACE_PALETTES):      # the "Custom" row
            return
        _, pos, neg = SURFACE_PALETTES[index]
        self.viewport.set_isosurface(positive_color=pos, negative_color=neg)
        self._refresh_swatches()

    def _pick_color(self, positive: bool):
        style = self.viewport.style
        rgb = style.surface_positive if positive else style.surface_negative
        color = QColorDialog.getColor(QColor.fromRgbF(*rgb), self,
                                      "Isosurface color")
        if not color.isValid():
            return
        chosen = (color.redF(), color.greenF(), color.blueF())
        if positive:
            self.viewport.set_isosurface(positive_color=chosen)
        else:
            self.viewport.set_isosurface(negative_color=chosen)
        self.palette_box.blockSignals(True)
        self._select_custom_entry()
        self.palette_box.blockSignals(False)
        self._refresh_swatches()

    def _select_custom_entry(self):
        if self.palette_box.findText("Custom") < 0:
            self.palette_box.addItem("Custom")
        self.palette_box.setCurrentIndex(self.palette_box.findText("Custom"))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("microscp")
        self.result: ParseResult | None = None
        self.current_frame = 0
        self._surface_dialog: SurfaceDialog | None = None

        self.viewport = MoleculeViewport(self)
        self.setCentralWidget(self.viewport)
        self.setAcceptDrops(True)

        self._edited = False
        self._file_label = QLabel("Open a file to begin  (File → Open…)")
        self.statusBar().addPermanentWidget(self._file_label)
        self.viewport.selectionChanged.connect(self.statusBar().showMessage)
        self.viewport.structureEdited.connect(self._mark_edited)

        self.spectra_dock = SpectraDock(self.viewport, self)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.spectra_dock)
        self.spectra_dock.hide()

        self._build_menus()
        self._build_frame_bar()

    # ------------------------------------------------------------------ UI

    def _build_menus(self):
        file_menu = self.menuBar().addMenu("&File")
        self._add_action(file_menu, "&Open…", QKeySequence.StandardKey.Open, self.open_dialog)
        self._add_action(file_menu, "&Save As…", QKeySequence.StandardKey.Save, self.save_dialog)
        self._add_action(file_menu, "&Export Image…", "Ctrl+E", self.export_image)
        file_menu.addSeparator()
        self._add_action(file_menu, "&Quit", QKeySequence.StandardKey.Quit, self.close)

        edit_menu = self.menuBar().addMenu("&Edit")
        self._add_action(edit_menu, "&Undo", QKeySequence.StandardKey.Undo, self._undo)
        self._add_action(edit_menu, "&Redo", QKeySequence.StandardKey.Redo, self._redo)
        edit_menu.addSeparator()
        self._add_action(edit_menu, "&Adjust Selection…", "E", self._adjust_selection)
        self._add_action(edit_menu, "Delete Selected Atoms", "X", self._delete_atoms)
        self._add_action(edit_menu, "Recompute &Bonds", "B", self._recompute_bonds)

        view_menu = self.menuBar().addMenu("&View")
        self._add_action(view_menu, "&Reset View", "Ctrl+R", self.viewport.reset_view)
        view_menu.addSeparator()

        repr_menu = view_menu.addMenu("Re&presentation")
        repr_group = QActionGroup(self)
        self._repr_actions = {}
        for name, text in (("cylview", "&CYLview"), ("houk", "&Houk (Houkmol)")):
            action = QAction(text, self, checkable=True)
            action.triggered.connect(lambda _=False, n=name: self._set_representation(n))
            repr_group.addAction(action)
            repr_menu.addAction(action)
            self._repr_actions[name] = action
        self._repr_actions["cylview"].setChecked(True)
        self._add_action(repr_menu, "&Toggle Representation", "V",
                         self._toggle_representation)
        repr_menu.addSeparator()
        for rep, text, key in ((REP_BALL, "Selection → &Ball && Stick", "1"),
                               (REP_STICK, "Selection → &Stick", "2"),
                               (REP_LINE, "Selection → &Line", "3")):
            # swallow QAction.triggered's checked argument (it would land in r)
            self._add_action(repr_menu, text, key,
                             lambda _=False, r=rep: self._apply_atom_rep(r))

        labels_menu = view_menu.addMenu("Atom &Labels")
        group = QActionGroup(self)
        self._label_actions = {}
        for mode, text in (("none", "&None"), ("element", "&Element"),
                           ("element+number", "Element + Num&ber"),
                           ("number", "N&umber")):
            action = QAction(text, self, checkable=True)
            action.triggered.connect(lambda _=False, m=mode: self.viewport.set_label_mode(m))
            group.addAction(action)
            labels_menu.addAction(action)
            self._label_actions[mode] = action
        self._label_actions["none"].setChecked(True)
        self._add_action(labels_menu, "&Cycle Labels", "L", self._cycle_labels)

        self._hbond_action = QAction("Show &Hydrogen Bonds", self, checkable=True)
        self._hbond_action.setChecked(True)
        self._hbond_action.setShortcut(QKeySequence("H"))
        self._hbond_action.toggled.connect(self._toggle_hbonds)
        view_menu.addAction(self._hbond_action)
        self._add_action(view_menu, "&Isosurface…", "I", self._show_surface_dialog)

        view_menu.addSeparator()
        self._add_action(view_menu, "&Align View to Selection", "A", self._align_view)
        self._add_action(view_menu, "&Center on Selected Atom", "C", self._center_atom)
        self._add_action(view_menu, "Center on &Molecule", "Home", self._center_molecule)

        measure_menu = self.menuBar().addMenu("&Measure")
        self._add_action(measure_menu, "&Pin Measurement", "M", self._pin_measurement)
        self._add_action(measure_menu, "&Clear Pinned Measurements", "Shift+M",
                         self._clear_pinned)

        spectra_menu = self.menuBar().addMenu("&Spectra")
        self._add_action(spectra_menu, "&Show/Hide Spectra", "S", self._toggle_spectra)
        self._add_action(spectra_menu, "Stop &Animation", "Space",
                         self.viewport.stop_animation)

        help_menu = self.menuBar().addMenu("&Help")
        self._add_action(help_menu, "&About", None, self.about)

    def _add_action(self, menu, text, shortcut, slot):
        action = QAction(text, self)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
        action.triggered.connect(slot)
        menu.addAction(action)
        return action

    def _build_frame_bar(self):
        self._frame_bar = QToolBar("Frames")
        self._frame_bar.setMovable(False)
        self.addToolBar(Qt.ToolBarArea.BottomToolBarArea, self._frame_bar)
        self._play_btn = QPushButton("▶")
        self._play_btn.setFixedWidth(34)
        self._play_btn.setCheckable(True)
        self._play_btn.toggled.connect(self._toggle_play)
        self._play_timer = QTimer(self)
        self._play_timer.setInterval(150)
        self._play_timer.timeout.connect(self._play_step)
        self._frame_slider = QSlider(Qt.Orientation.Horizontal)
        self._frame_slider.setMinimumWidth(300)
        self._frame_slider.valueChanged.connect(self._frame_changed)
        self._frame_label = QLabel("")
        self._frame_bar.addWidget(self._play_btn)
        self._frame_bar.addWidget(self._frame_slider)
        self._frame_bar.addWidget(self._frame_label)
        self._frame_bar.setVisible(False)

    def _toggle_play(self, playing: bool):
        self._play_btn.setText("⏸" if playing else "▶")
        if playing:
            self._play_timer.start()
        else:
            self._play_timer.stop()

    def _play_step(self):
        if self.result is None or self.result.nframes < 2:
            return
        self._frame_slider.setValue((self._frame_slider.value() + 1)
                                    % self.result.nframes)

    # ------------------------------------------------------------------ files

    def open_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open molecular file", "", OPEN_FILTER)
        if path:
            self.open_file(path)

    def open_file(self, path):
        try:
            result = mio.load(path)
        except Exception as exc:
            QMessageBox.critical(self, "Could not open file", f"{path}\n\n{exc}")
            return
        self.result = result
        self.current_frame = result.nframes - 1
        self._edited = False
        if self._surface_dialog is not None:
            self._surface_dialog.close()
            self._surface_dialog = None
        mol = result.frames[self.current_frame]
        self.viewport.set_molecule(mol)
        self.viewport.clear_history()
        self.spectra_dock.set_result(result)
        if result.volumes:
            self.viewport.set_volume(result.volumes[0])
            self.statusBar().showMessage(
                f"Isosurface shown at ±{self.viewport.isovalue:g} — press I to adjust",
                6000)
        name = Path(path).name
        self.setWindowTitle(f"microscp — {name}")
        info = f"{name}  ·  {mol.formula()}  ·  {mol.natoms} atoms  ·  {result.program}"
        if result.normal_termination is False:
            info += "  ·  ⚠ abnormal termination"
        self._file_label.setText(info)

        multi = result.nframes > 1
        self._frame_bar.setVisible(multi)
        self._play_btn.setChecked(False)
        if multi:
            self._frame_slider.blockSignals(True)
            self._frame_slider.setRange(0, result.nframes - 1)
            self._frame_slider.setValue(self.current_frame)
            self._frame_slider.blockSignals(False)
            self._update_frame_label()

    def _frame_changed(self, value: int):
        if self.result is None:
            return
        self.current_frame = int(value)
        self.viewport.set_molecule(self.result.frames[self.current_frame], keep_camera=True)
        self.viewport.clear_history()
        self._update_frame_label()

    def _update_frame_label(self):
        n = self.result.nframes if self.result else 0
        text = f"  frame {self.current_frame + 1}/{n}"
        energies = self.result.scf_energies if self.result else []
        if self.current_frame < len(energies):
            text += f"  ·  E = {energies[self.current_frame]:.6f} Ha"
        self._frame_label.setText(text)

    def save_dialog(self):
        if self.viewport.molecule is None:
            QMessageBox.information(self, "Nothing to save", "Open a file first.")
            return
        path, chosen = QFileDialog.getSaveFileName(self, "Save structure", "",
                                                   SAVE_FILTER)
        if not path:
            return
        if not Path(path).suffix:
            for name, ext in _FILTER_DEFAULT_EXT.items():
                if chosen.startswith(name):
                    path += ext
                    break
        try:
            mio.save_molecule(path, self.viewport.molecule)
            self.statusBar().showMessage(f"Saved {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Could not save", str(exc))

    def export_image(self):
        if self.viewport.molecule is None:
            QMessageBox.information(self, "Nothing to export", "Open a file first.")
            return
        dialog = ExportImageDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        path, chosen = QFileDialog.getSaveFileName(self, "Export image",
                                                   "molecule.png", EXPORT_FILTER)
        if not path:
            return
        if not Path(path).suffix:
            for name, ext in _EXPORT_DEFAULT_EXT.items():
                if chosen.startswith(name):
                    path += ext
                    break
        scale = dialog.scale.value()
        vp = self.viewport
        try:
            volume = vp.volume if (vp.has_volume and vp.surface_visible) else None
            image = render_molecule_image(
                vp.molecule, vp.style, vp.camera,
                vp.width() * scale, vp.height() * scale,
                supersample=2, transparent=dialog.transparent.isChecked(),
                volume=volume, isovalue=vp.isovalue, reps=vp.atom_reps)
            if dialog.annotations.isChecked() and (vp.pinned or vp.label_mode != "none"):
                painter = QPainter(image)
                draw_annotations(painter, vp.molecule, vp.camera,
                                 image.width(), image.height(),
                                 vp.pinned, vp.label_mode, scale=scale)
                painter.end()
            _write_image(image, path)
            self.statusBar().showMessage(f"Exported {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    # ------------------------------------------------------------------ view/measure

    def _cycle_labels(self):
        mode = self.viewport.cycle_label_mode()
        self._label_actions[mode].setChecked(True)
        self.statusBar().showMessage(f"Atom labels: {mode}", 3000)

    def _align_view(self):
        msg = self.viewport.align_to_selection()
        self.statusBar().showMessage(
            msg or "Select 2 atoms (view along bond) or 3 atoms (plane to screen) first",
            4000)

    def _center_atom(self):
        ok = self.viewport.center_on_selection()
        self.statusBar().showMessage(
            "Rotation center set to the selected atom" if ok else "Select an atom first",
            4000)

    def _center_molecule(self):
        self.viewport.center_on_molecule()
        self.statusBar().showMessage("Rotation center reset to the molecule", 3000)

    def _pin_measurement(self):
        ok = self.viewport.pin_selection()
        self.statusBar().showMessage(
            "Measurement pinned" if ok else "Select 2–4 atoms first", 4000)

    def _clear_pinned(self):
        self.viewport.clear_pinned()
        self.statusBar().showMessage("Pinned measurements cleared", 3000)

    def _adjust_selection(self):
        n = len(self.viewport.selection)
        if self.viewport.molecule is None or not 2 <= n <= 4:
            self.statusBar().showMessage(
                "Select 2 atoms (distance), 3 (angle), or 4 (dihedral) first", 4000)
            return
        AdjustDialog(self.viewport, self).exec()

    def _delete_atoms(self):
        count = len(self.viewport.selection)
        if self.viewport.delete_selected():
            self.statusBar().showMessage(f"Deleted {count} atom(s)", 4000)
        else:
            self.statusBar().showMessage("Select atoms to delete first", 4000)

    def _recompute_bonds(self):
        self.viewport.recompute_bonds()
        self.statusBar().showMessage("Bonds recomputed from covalent radii", 3000)

    def _undo(self):
        self.statusBar().showMessage(
            "Undone" if self.viewport.undo() else "Nothing to undo", 2500)

    def _redo(self):
        self.statusBar().showMessage(
            "Redone" if self.viewport.redo() else "Nothing to redo", 2500)

    def _mark_edited(self):
        if not self._edited:
            self._edited = True
            self._file_label.setText(self._file_label.text() + "  ·  ✎ edited")

    def _set_representation(self, name: str):
        self.viewport.set_representation(name)
        self._repr_actions[name].setChecked(True)
        label = "CYLview" if name == "cylview" else "Houk (Houkmol)"
        self.statusBar().showMessage(f"Representation: {label}", 3000)

    def _toggle_representation(self):
        current = self.viewport.style.name
        self._set_representation("houk" if current == "cylview" else "cylview")

    def _apply_atom_rep(self, rep: int):
        vp = self.viewport
        if vp.molecule is None:
            return
        target = list(vp.selection)
        vp.set_atom_representation(rep, target)
        name = {REP_BALL: "ball & stick", REP_STICK: "stick", REP_LINE: "line"}[rep]
        scope = f"{len(target)} selected atom(s)" if target else "all atoms"
        self.statusBar().showMessage(f"Representation of {scope}: {name}", 3000)

    def _toggle_hbonds(self, checked: bool):
        self.viewport.style.show_hbonds = checked
        if self.viewport.molecule is not None:
            self.viewport._rebuild_scene()
        self.statusBar().showMessage(
            f"Hydrogen bonds {'shown' if checked else 'hidden'}", 3000)

    def _show_surface_dialog(self):
        if not self.viewport.has_volume:
            self.statusBar().showMessage(
                "No volumetric data — open a cube file (.cube/.cub) first", 4000)
            return
        if self._surface_dialog is None:
            volumes = self.result.volumes if self.result else [self.viewport.volume]
            self._surface_dialog = SurfaceDialog(self.viewport, volumes, self)
            self._surface_dialog.finished.connect(self._surface_dialog_closed)
        self._surface_dialog.show()
        self._surface_dialog.raise_()
        self._surface_dialog.activateWindow()

    def _surface_dialog_closed(self):
        self._surface_dialog = None

    def _toggle_spectra(self):
        if not self.spectra_dock.has_data:
            self.statusBar().showMessage(
                "No spectral data in this file (need a freq, TD-DFT, or NMR job)", 4000)
            return
        self.spectra_dock.setVisible(not self.spectra_dock.isVisible())

    # ------------------------------------------------------------------ misc

    def about(self):
        QMessageBox.about(
            self, "About microscp",
            f"<b>microscp {__version__}</b> — a microscope for molecules<br>"
            "Structure and spectroscopy viewer for quantum chemistry.<br><br>"
            "Rotate: left-drag · Pan: right-drag · Zoom: scroll<br>"
            "Measure: click 2–4 atoms · Esc clears · M pins the measurement<br>"
            "L cycles atom labels · A aligns view to selection (2=bond, 3=plane)<br>"
            "C centers rotation on the selected atom · Home re-centers the molecule<br>"
            "E adjusts the selected distance/angle/dihedral · X deletes atoms<br>"
            "Ctrl+Z / Ctrl+Shift+Z undo/redo · S spectra · Space stops animation<br>"
            "I opens isosurface controls for cube files (MOs, densities)")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            self.open_file(url.toLocalFile())
            break


def main():
    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    fmt.setSamples(8)
    QSurfaceFormat.setDefaultFormat(fmt)

    app = QApplication(sys.argv)
    app.setApplicationName("microscp")
    window = MainWindow()
    window.resize(1000, 720)
    window.show()
    if len(sys.argv) > 1:
        window.open_file(sys.argv[1])
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
