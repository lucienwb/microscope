"""The viewer window: menus, shortcuts, and what each command does."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QPainter
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSlider,
    QStackedWidget,
    QToolBar,
)

from .. import __version__
from .. import io as mio
from ..core.results import ParseResult
from ..render.formats import IMAGE_SUFFIXES, VECTOR_SUFFIXES
from ..render.imagefile import write_image as _write_image
from ..render.lewisdraw import render_lewis_image, write_lewis_vector
from ..render.offscreen import render_molecule_image
from ..render.scene import REP_BALL, REP_LINE, REP_STICK
from .annotations import draw_annotations
from .dialogs import AdjustDialog, ExportImageDialog, StyleDialog, SurfaceDialog
from .filetypes import (
    _CHEMDRAW_DEFAULT_EXT,
    _EXPORT_DEFAULT_EXT,
    _FILTER_DEFAULT_EXT,
    CHEMDRAW_FILTER,
    EXPORT_FILTER,
    LEWIS_EXPORT_FILTER,
    OPEN_FILTER,
    SAVE_FILTER,
    with_extension,
)
from .lewisview import LewisView
from .menus import build_menus
from .spectra import SpectraDock
from .viewport import MoleculeViewport


def _charge_word(structure) -> str:
    """How to name what is hanging off the brackets, for a status message."""
    parts = []
    if structure.net_charge:
        parts.append(f"{structure.net_charge:+d} charge")
    if structure.net_radicals:
        parts.append("radical" if structure.net_radicals == 1 else "radicals")
    return " and ".join(parts)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("microscope")
        self.result: ParseResult | None = None
        self.current_frame = 0
        self._surface_dialog: SurfaceDialog | None = None
        self._style_dialog: StyleDialog | None = None

        self.viewport = MoleculeViewport(self)
        # the flat view shares the camera, so both stay on the same orientation
        self.lewis = LewisView(self)
        self.lewis.set_camera(self.viewport.camera)
        self._pages = QStackedWidget(self)
        self._pages.addWidget(self.viewport)
        self._pages.addWidget(self.lewis)
        self.setCentralWidget(self._pages)
        self.setAcceptDrops(True)

        self._edited = False
        self._file_label = QLabel("Open a file to begin  (File → Open…)")
        self.statusBar().addPermanentWidget(self._file_label)
        self.viewport.selectionChanged.connect(self.statusBar().showMessage)
        self.viewport.structureEdited.connect(self._mark_edited)

        self.spectra_dock = SpectraDock(self.viewport, self)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.spectra_dock)
        self.spectra_dock.hide()

        build_menus(self)
        self._build_frame_bar()

    # ------------------------------------------------------------------ UI


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
        self.lewis.set_molecule(mol)
        self.viewport.clear_history()
        self.spectra_dock.set_result(result)
        if result.volumes:
            self.viewport.set_volume(result.volumes[0])
            self.statusBar().showMessage(
                f"Isosurface shown at ±{self.viewport.isovalue:g} — press I to adjust",
                6000)
        name = Path(path).name
        self.setWindowTitle(f"microscope — {name}")
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
        frame = self.result.frames[self.current_frame]
        self.viewport.set_molecule(frame, keep_camera=True)
        self.lewis.set_molecule(frame)
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
        path = with_extension(path, chosen, _FILTER_DEFAULT_EXT)
        try:
            mio.save_molecule(path, self.viewport.molecule)
            self.statusBar().showMessage(f"Saved {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Could not save", str(exc))

    @property
    def lewis_active(self) -> bool:
        return self._pages.currentWidget() is self.lewis

    def export_image(self):
        if self.viewport.molecule is None:
            QMessageBox.information(self, "Nothing to export", "Open a file first.")
            return
        if self.lewis_active:
            self._export_lewis()
            return
        dialog = ExportImageDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        path, chosen = QFileDialog.getSaveFileName(self, "Export image",
                                                   "molecule.png", EXPORT_FILTER)
        if not path:
            return
        path = with_extension(path, chosen, _EXPORT_DEFAULT_EXT)
        scale = dialog.scale.value()
        vp = self.viewport
        try:
            volume = vp.volume if (vp.has_volume and vp.surface_visible) else None
            image = render_molecule_image(
                vp.molecule, vp.style, vp.camera,
                vp.width() * scale, vp.height() * scale,
                supersample=2, transparent=dialog.transparent.isChecked(),
                volume=volume, isovalue=vp.isovalue, reps=vp.atom_reps)
            if dialog.annotations.isChecked() and (vp.pinned or vp.show_axes
                                                   or vp.label_mode != "none"):
                painter = QPainter(image)
                draw_annotations(painter, vp.molecule, vp.camera,
                                 image.width(), image.height(),
                                 vp.pinned, vp.label_mode, scale=scale,
                                 axes=vp.show_axes)
                painter.end()
            _write_image(image, path)
            self.statusBar().showMessage(f"Exported {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _export_lewis(self):
        """Export the flat drawing — as pixels, or as the vector art it is."""
        dialog = ExportImageDialog(self)
        dialog.annotations.setEnabled(False)      # a drawing has none to add
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        path, chosen = QFileDialog.getSaveFileName(self, "Export drawing",
                                                   "structure.png",
                                                   LEWIS_EXPORT_FILTER)
        if not path:
            return
        path = with_extension(path, chosen, _EXPORT_DEFAULT_EXT)
        view = self.lewis
        transparent = dialog.transparent.isChecked()
        suffix = Path(path).suffix.lower()
        try:
            if suffix in VECTOR_SUFFIXES:
                # vector art has no resolution to scale, so the page is the view
                write_lewis_vector(path, view.structure, view.camera,
                                   view.width(), view.height(), view.options,
                                   None if transparent else (1.0, 1.0, 1.0))
            elif suffix in IMAGE_SUFFIXES:
                scale = dialog.scale.value()
                image = render_lewis_image(view.structure, view.camera,
                                           view.width() * scale,
                                           view.height() * scale, view.options,
                                           transparent=transparent)
                _write_image(image, path)
            else:
                what = suffix or "a file with no extension"
                raise ValueError(f"cannot write {what}; "
                                 "use .png, .tif, .svg or .pdf")
            self.statusBar().showMessage(f"Exported {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def save_chemdraw(self):
        """Write the drawing as a structure file to carry on editing elsewhere."""
        if self.viewport.molecule is None:
            QMessageBox.information(self, "Nothing to save", "Open a file first.")
            return
        view = self.lewis
        # re-read the geometry: it may have been edited since the last look
        view.set_molecule(self.viewport.molecule)
        path, chosen = QFileDialog.getSaveFileName(self, "Save for ChemDraw",
                                                   "structure.cdxml",
                                                   CHEMDRAW_FILTER)
        if not path:
            return
        path = with_extension(path, chosen, _CHEMDRAW_DEFAULT_EXT)
        try:
            mio.save_lewis(path, view.structure,
                           view.projected(view.width(), view.height()),
                           view.drawn_atoms())
        except Exception as exc:
            QMessageBox.critical(self, "Could not save", str(exc))
            return
        message = f"Saved {path} — opens at this angle"
        if view.structure.bracketed:
            # neither format has a slot for a charge that sits on no atom
            message += (f"  ·  ⚠ the {_charge_word(view.structure)} on the "
                        f"brackets is not carried into a structure file")
        if not view.structure.matches_file:
            message += (f"  ·  ⚠ perceived charge "
                        f"{view.structure.total_charge:+d}, the file says "
                        f"{int(view.molecule.charge):+d}")
        self.statusBar().showMessage(message, 8000)

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

    def _select_fragment(self):
        n = self.viewport.select_fragment()
        self.statusBar().showMessage(
            f"Selected the connected fragment ({n} atoms)" if n
            else "Select an atom of the fragment first", 4000)

    def _toggle_gizmo(self, on: bool):
        self.viewport.set_show_gizmo(on)
        self.statusBar().showMessage(
            "Move/rotate handles on — drag an arrow, a ring, or the centre dot "
            "(Shift snaps to 0.1 Å / 15°)" if on
            else "Move/rotate handles off", 4000)

    def _reset_view(self):
        self.viewport.reset_view()
        self.lewis.update()

    def _toggle_lewis(self, on: bool):
        if on and self.viewport.molecule is None:
            self._lewis_action.blockSignals(True)     # do not come back here
            self._lewis_action.setChecked(False)
            self._lewis_action.blockSignals(False)
            self.statusBar().showMessage("Open a file first", 4000)
            return
        if on:
            self.viewport.stop_animation()
            self.lewis.set_molecule(self.viewport.molecule)
        self._pages.setCurrentWidget(self.lewis if on else self.viewport)
        for action in self._edit_actions:
            action.setEnabled(not on)                 # this mode is read-only
        widget = self.lewis if on else self.viewport
        widget.setFocus()
        self.statusBar().showMessage(
            "Lewis structure — drag to turn it, Shift+drag to spin it in the "
            "page, scroll to zoom. Read-only: pick the angle, then export or "
            "save for ChemDraw." if on else "Back to the 3D view", 8000)

    def _toggle_axes(self, on: bool):
        self.viewport.set_show_axes(on)
        self.statusBar().showMessage(
            f"XYZ axes {'shown' if on else 'hidden'}", 3000)

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

    def edit_style(self):
        """Open the style editor — live, and savable as a file to share."""
        if self._style_dialog is None:
            self._style_dialog = StyleDialog(self.viewport, self)
        self._style_dialog.show()
        self._style_dialog.raise_()

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
            self, "About microscope",
            f"<b>microscope {__version__}</b> — a microscope for "
            "molecules (command: <code>scope</code>)<br>"
            "Structure and spectroscopy viewer for quantum chemistry.<br><br>"
            "Rotate: left-drag · Pan: right-drag · Zoom: scroll<br>"
            "Measure: click 2–4 atoms · Esc clears · M pins the measurement<br>"
            "L cycles atom labels · A aligns view to selection (2=bond, 3=plane)<br>"
            "C centers rotation on the selected atom · Home re-centers the molecule<br>"
            "E adjusts the selected distance/angle/dihedral · X deletes atoms<br>"
            "Ctrl+Z / Ctrl+Shift+Z undo/redo · S spectra · Space stops animation<br>"
            "I opens isosurface controls for cube files (MOs, densities)<br>"
            "Shift+L draws the flat Lewis structure (ChemDraw style)")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            self.open_file(url.toLocalFile())
            break
