"""The modal dialogs: image export, geometry adjustment, isosurface setup.

Each one edits a viewport it is handed and nothing else, so they can be read
and changed without the main window in view.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QFontDatabase
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPushButton,
    QSlider,
    QSpinBox,
    QWidget,
)

from .. import io as mio
from ..core import elements, geometry
from ..core.orbitals import HARTREE_TO_EV
from ..render.styles import StyleError, load_style, save_style
from .filetypes import CUBE_FILTER, STYLE_FILTER
from .viewport import MoleculeViewport

# Curated isosurface color pairs (+ lobe, - lobe); custom colors via the picker.
SURFACE_PALETTES = (
    ("Blue / Red", (0.29, 0.44, 0.86), (0.88, 0.38, 0.22)),
    ("Red / Blue", (0.88, 0.38, 0.22), (0.29, 0.44, 0.86)),
    ("Teal / Orange", (0.13, 0.59, 0.62), (0.95, 0.52, 0.16)),
    ("Purple / Gold", (0.55, 0.36, 0.76), (0.93, 0.69, 0.13)),
    ("Green / Magenta", (0.30, 0.63, 0.36), (0.79, 0.29, 0.62)),
    ("Slate / Silver", (0.36, 0.42, 0.52), (0.72, 0.75, 0.80)),
)


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
    """Modeless isosurface controls (live update): the grids a cube file holds,
    or the orbitals of a wavefunction file, put on a grid as they are chosen.

    ``orbital_volume(spin, index)`` makes the grid for one orbital; the window
    passes its own, which keeps the last few so stepping back is instant.
    """

    DEBOUNCE_MS = 120      # holding an arrow key moves the list, not the CPU

    def __init__(self, viewport: MoleculeViewport, volumes: list, parent=None,
                 orbitals=None, orbital_volume=None):
        super().__init__(parent)
        self.viewport = viewport
        self.volumes = volumes
        self.orbitals = orbitals
        self.orbital_volume = orbital_volume
        self.setWindowTitle("Orbitals & Isosurface" if orbitals else "Isosurface")
        layout = QFormLayout(self)

        if orbitals is not None:
            self._build_orbital_list(layout)
        elif len(volumes) > 1:
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
        save = buttons.addButton("Save Cube…", QDialogButtonBox.ButtonRole.ActionRole)
        save.setToolTip("Write the grid on show as a Gaussian cube file")
        save.clicked.connect(self._save_cube)
        buttons.rejected.connect(self.close)
        layout.addRow(buttons)

        if orbitals is not None:
            self._select_start()

    def _switch_volume(self, index: int):
        self.viewport.set_volume(self.volumes[index], isovalue=self.iso.value())

    # ---------------------------------------------------------------- orbitals

    def _build_orbital_list(self, layout: QFormLayout):
        if self.orbitals.convention == "unrecognized":
            warning = QLabel("⚠ This file's basis-set conventions were not recognized:\n"
                             "its orbitals may be drawn wrong.")
            warning.setStyleSheet("color: #b35c00;")
            layout.addRow(warning)
        if not self.orbitals.restricted:
            self.spin_box = QComboBox()
            self.spin_box.addItems(["α (alpha)", "β (beta)"])
            self.spin_box.currentIndexChanged.connect(self._spin_changed)
            layout.addRow("Spin:", self.spin_box)
        self.orbital_list = QListWidget()
        self.orbital_list.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.orbital_list.setMinimumSize(340, 260)
        self.orbital_list.currentRowChanged.connect(self._orbital_row_changed)
        layout.addRow(self.orbital_list)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._show_chosen_orbital)
        self._spin = "alpha"
        self._fill_orbitals()

    def _fill_orbitals(self):
        orbitals = self.orbitals.spin(self._spin)
        self.orbital_list.blockSignals(True)
        self.orbital_list.clear()
        for index in range(orbitals.nmo):
            energy = float(orbitals.energies[index]) * HARTREE_TO_EV
            shown = f"{energy:>10.2f} eV" if np.isfinite(energy) else " " * 13
            text = (f"{self.orbitals.name(self._spin, index):<8} {index + 1:>5}"
                    f"{shown}   occ {float(orbitals.occupations[index]):g}")
            if index < len(orbitals.symmetries) and orbitals.symmetries[index]:
                text += f"   {orbitals.symmetries[index]}"
            self.orbital_list.addItem(text)
        self.orbital_list.blockSignals(False)

    def _select_start(self):
        """The orbital on show if it is one, else the HOMO."""
        meta = self.viewport.volume.meta if self.viewport.has_volume else {}
        if "orbital" in meta:
            if meta.get("spin") == "beta" and hasattr(self, "spin_box"):
                self.spin_box.setCurrentIndex(1)
            row = meta["orbital"]
        else:
            row = max(self.orbitals.spin(self._spin).homo, 0)
        self.orbital_list.setCurrentRow(row)
        self.orbital_list.scrollToItem(self.orbital_list.item(row),
                                       QListWidget.ScrollHint.PositionAtCenter)
        self._timer.stop()
        if "orbital" not in meta:
            self._show_chosen_orbital()

    def _spin_changed(self, index: int):
        self._spin = "beta" if index == 1 else "alpha"
        row = self.orbital_list.currentRow()
        self._fill_orbitals()
        self.orbital_list.setCurrentRow(min(max(row, 0), self.orbital_list.count() - 1))

    def _orbital_row_changed(self, _row: int):
        self._timer.start(self.DEBOUNCE_MS)

    def _show_chosen_orbital(self):
        row = self.orbital_list.currentRow()
        if row < 0 or self.orbital_volume is None:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            volume = self.orbital_volume(self._spin, row)
        finally:
            QApplication.restoreOverrideCursor()
        keep = self.viewport.isovalue if self.viewport.has_volume else None
        self.viewport.set_volume(volume, isovalue=keep)
        for widget, set_value in ((self.iso, lambda: self.iso.setValue(self.viewport.isovalue)),
                                  (self.visible, lambda: self.visible.setChecked(True))):
            widget.blockSignals(True)       # the viewport already has both
            set_value()
            widget.blockSignals(False)

    def _save_cube(self):
        if not self.viewport.has_volume or self.viewport.molecule is None:
            return
        volume = self.viewport.volume
        name = volume.label.split(" · ")[0].replace(" ", "_") or "grid"
        path, _ = QFileDialog.getSaveFileName(self, "Save Cube File", f"{name}.cube",
                                              CUBE_FILTER)
        if not path:
            return
        if not Path(path).suffix:
            path += ".cube"
        try:
            mio.save_cube(path, volume, self.viewport.molecule)
        except OSError as exc:
            QMessageBox.warning(self, "Save Cube File", f"Could not write {path}:\n{exc}")

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

def _swatch(button: QPushButton, rgb) -> None:
    """Paint a colour button with the colour it stands for."""
    colour = QColor.fromRgbF(*rgb) if rgb is not None else QColor(255, 255, 255)
    edge = "#888" if rgb is not None else "#ccc"
    button.setStyleSheet(
        f"background-color: {colour.name()}; border: 1px solid {edge};")


class StyleDialog(QDialog):
    """Edit the rendering style, live, and keep it as a file.

    The point of the file is a group standard: save it once, then both the
    viewer and `scope -s --style ours.json` draw every figure the same way.
    Only the elements actually in the open structure get a colour button —
    a periodic table of them would be a worse way to find carbon.
    """

    def __init__(self, viewport: MoleculeViewport, parent=None):
        super().__init__(parent)
        self.viewport = viewport
        self.setWindowTitle("Style")
        layout = QFormLayout(self)

        self.preset = QComboBox()
        self.preset.addItems(["cylview", "houk"])
        if viewport.style.name in ("cylview", "houk"):
            self.preset.setCurrentText(viewport.style.name)
        else:
            self.preset.addItem(viewport.style.name)
            self.preset.setCurrentText(viewport.style.name)
        self.preset.activated.connect(self._preset_chosen)
        layout.addRow("Start from:", self.preset)

        self.atom_scale = self._slider(layout, "Atom size:", 10, 100,
                                       viewport.style.atom_scale * 100,
                                       lambda v: self._set(atom_scale=v / 100.0))
        self.bond_radius = self._slider(layout, "Bond width:", 2, 40,
                                        viewport.style.bond_radius * 100,
                                        lambda v: self._set(bond_radius=v / 100.0))

        self.background = QPushButton()
        self.background.setFixedSize(56, 22)
        self.background.clicked.connect(self._pick_background)
        layout.addRow("Background:", self.background)

        self.uniform_bonds = QCheckBox("one colour for every bond")
        self.uniform_bonds.setChecked(viewport.style.bond_color is not None)
        self.uniform_bonds.toggled.connect(self._uniform_toggled)
        self.bond_colour = QPushButton()
        self.bond_colour.setFixedSize(56, 22)
        self.bond_colour.clicked.connect(self._pick_bond_colour)
        bonds = QWidget()
        row = QHBoxLayout(bonds)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.uniform_bonds)
        row.addWidget(self.bond_colour)
        row.addStretch(1)
        layout.addRow("Bonds:", bonds)

        self.hbonds = QCheckBox("show hydrogen bonds")
        self.hbonds.setChecked(viewport.style.show_hbonds)
        self.hbonds.toggled.connect(lambda on: self._set(show_hbonds=bool(on)))
        layout.addRow("", self.hbonds)

        self.elements = QWidget()
        self.element_row = QHBoxLayout(self.elements)
        self.element_row.setContentsMargins(0, 0, 0, 0)
        layout.addRow("Elements:", self.elements)
        self._build_element_buttons()

        buttons = QDialogButtonBox()
        buttons.addButton("&Load…", QDialogButtonBox.ButtonRole.ResetRole
                          ).clicked.connect(self._load)
        buttons.addButton("&Save As…", QDialogButtonBox.ButtonRole.ApplyRole
                          ).clicked.connect(self._save)
        buttons.addButton(QDialogButtonBox.StandardButton.Close
                          ).clicked.connect(self.close)
        layout.addRow(buttons)
        self._refresh()

    # ------------------------------------------------------------------ parts

    def _slider(self, layout, label, low, high, value, on_change):
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(low, high)
        slider.setValue(int(round(value)))
        slider.valueChanged.connect(on_change)
        layout.addRow(label, slider)
        return slider

    def _build_element_buttons(self) -> None:
        while self.element_row.count():
            self.element_row.takeAt(0).widget().deleteLater()
        self.element_buttons = {}
        molecule = self.viewport.molecule
        present = sorted({int(z) for z in molecule.atomic_numbers}) if molecule else []
        for z in present[:12]:                      # a long row helps nobody
            button = QPushButton(elements.SYMBOLS[z])
            button.setFixedSize(38, 22)
            button.clicked.connect(lambda _=False, z=z: self._pick_element(z))
            self.element_row.addWidget(button)
            self.element_buttons[z] = button
        self.element_row.addStretch(1)

    # ----------------------------------------------------------------- edits

    def _set(self, **changes) -> None:
        for key, value in changes.items():
            setattr(self.viewport.style, key, value)
        self.viewport.refresh_style()
        self._refresh()

    def _preset_chosen(self, _index: int) -> None:
        name = self.preset.currentText()
        if name in ("cylview", "houk"):
            self.viewport.set_representation(name)
            self._sync_controls()

    def _uniform_toggled(self, on: bool) -> None:
        self._set(bond_color=(0.07, 0.07, 0.07) if on else None)

    def _pick_background(self) -> None:
        self._ask_colour(self.viewport.style.background,
                         lambda rgb: self._set(background=rgb))

    def _pick_bond_colour(self) -> None:
        if self.viewport.style.bond_color is None:
            self.uniform_bonds.setChecked(True)
        self._ask_colour(self.viewport.style.bond_color,
                         lambda rgb: self._set(bond_color=rgb))

    def _pick_element(self, z: int) -> None:
        def apply(rgb):
            palette = dict(self.viewport.style.palette)
            palette[z] = rgb
            self._set(palette=palette)
        self._ask_colour(self.viewport.style.atom_color(z), apply)

    def _ask_colour(self, current, apply) -> None:
        start = QColor.fromRgbF(*current) if current else QColor(0, 0, 0)
        chosen = QColorDialog.getColor(start, self, "Choose a colour")
        if chosen.isValid():
            apply((chosen.redF(), chosen.greenF(), chosen.blueF()))

    # ------------------------------------------------------------------ files

    def _save(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save style", "style.json",
                                              STYLE_FILTER)
        if not path:
            return
        style = self.viewport.style
        if style.name in ("cylview", "houk"):
            style.name = Path(path).stem
        save_style(path, style)
        self.preset.setCurrentText(style.name) if self.preset.findText(
            style.name) >= 0 else self.preset.addItem(style.name)
        self.preset.setCurrentText(style.name)

    def _load(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open style", "", STYLE_FILTER)
        if not path:
            return
        try:
            style = load_style(path)
        except StyleError as exc:
            QMessageBox.warning(self, "Could not read that style", str(exc))
            return
        self.viewport.apply_style(style)
        if self.preset.findText(style.name) < 0:
            self.preset.addItem(style.name)
        self.preset.setCurrentText(style.name)
        self._sync_controls()

    # ---------------------------------------------------------------- display

    def _sync_controls(self) -> None:
        style = self.viewport.style
        for slider, value in ((self.atom_scale, style.atom_scale * 100),
                              (self.bond_radius, style.bond_radius * 100)):
            slider.blockSignals(True)
            slider.setValue(int(round(value)))
            slider.blockSignals(False)
        self.uniform_bonds.blockSignals(True)
        self.uniform_bonds.setChecked(style.bond_color is not None)
        self.uniform_bonds.blockSignals(False)
        self.hbonds.blockSignals(True)
        self.hbonds.setChecked(style.show_hbonds)
        self.hbonds.blockSignals(False)
        self._refresh()

    def _refresh(self) -> None:
        style = self.viewport.style
        _swatch(self.background, style.background)
        _swatch(self.bond_colour, style.bond_color)
        self.bond_colour.setEnabled(style.bond_color is not None)
        for z, button in self.element_buttons.items():
            _swatch(button, style.atom_color(z))
