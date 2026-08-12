"""Main application window."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QKeySequence, QSurfaceFormat
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
    QLabel, QMainWindow, QMessageBox, QPushButton, QSlider, QSpinBox, QToolBar,
)

from .. import __version__
from .. import io as mio
from ..core import geometry
from ..core.results import ParseResult
from ..render.offscreen import render_molecule_image
from .spectra import SpectraDock
from .viewport import MoleculeViewport

OPEN_FILTER = (
    "Molecular files (*.xyz *.log *.out *.fchk *.fck *.fch *.gjf *.com *.gau "
    "*.pdb *.molden);;All files (*)"
)
SAVE_FILTER = "XYZ (*.xyz);;Gaussian input (*.gjf *.com);;PDB (*.pdb)"


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
        layout.addRow("Resolution scale:", self.scale)
        layout.addRow("Transparent background:", self.transparent)
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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("microscp")
        self.result: ParseResult | None = None
        self.current_frame = 0

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
        mol = result.frames[self.current_frame]
        self.viewport.set_molecule(mol)
        self.viewport.clear_history()
        self.spectra_dock.set_result(result)
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
        path, _ = QFileDialog.getSaveFileName(self, "Save structure", "", SAVE_FILTER)
        if not path:
            return
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
        path, _ = QFileDialog.getSaveFileName(self, "Export image", "molecule.png",
                                              "PNG image (*.png)")
        if not path:
            return
        scale = dialog.scale.value()
        try:
            image = render_molecule_image(
                self.viewport.molecule, self.viewport.style, self.viewport.camera,
                self.viewport.width() * scale, self.viewport.height() * scale,
                supersample=2, transparent=dialog.transparent.isChecked())
            image.save(path)
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
            "Ctrl+Z / Ctrl+Shift+Z undo/redo · S spectra · Space stops animation")

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
