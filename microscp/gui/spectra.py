"""Spectra dock: interactive IR / UV-Vis / NMR plots.

The IR tab is linked to the 3D viewport: clicking a band animates that
vibrational mode; clicking it again (or clicking elsewhere) stops it.
"""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox, QDockWidget, QDoubleSpinBox, QFileDialog, QHBoxLayout, QLabel,
    QPushButton, QTabWidget, QVBoxLayout, QWidget,
)

from ..core.results import ParseResult
from ..spectra import ir_spectrum, nmr_spectrum, uvvis_spectrum

_CURVE = "#303030"
_STICK = "#9a9a9a"
_ACTIVE = "#e08214"


class _SpectrumTab(QWidget):
    """Shared plumbing: control row, matplotlib canvas, CSV/PNG export."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(5.0, 2.6), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.controls = QHBoxLayout()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.addLayout(self.controls)
        layout.addWidget(self.canvas)
        self._spec = None

    def _add_control(self, label: str, widget) -> None:
        self.controls.addWidget(QLabel(label))
        self.controls.addWidget(widget)

    def _finish_controls(self, hint: str = "") -> None:
        if hint:
            lbl = QLabel(hint)
            lbl.setStyleSheet("color: gray;")
            self.controls.addWidget(lbl)
        self.controls.addStretch(1)
        for text, slot in (("CSV", self.export_csv), ("PNG", self.export_png)):
            btn = QPushButton(text)
            btn.setFixedWidth(52)
            btn.clicked.connect(slot)
            self.controls.addWidget(btn)

    def export_csv(self):
        if self._spec is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export spectrum data",
                                              "spectrum.csv", "CSV (*.csv)")
        if not path:
            return
        import pandas as pd
        pd.DataFrame({self._spec.xlabel: self._spec.x,
                      self._spec.ylabel: self._spec.y}).to_csv(path, index=False)

    def export_png(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figure", "spectrum.png",
            "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)")
        if path:
            self.figure.savefig(path, dpi=300)

    def _style_axes(self, spec):
        self.ax.set_xlabel(spec.xlabel)
        self.ax.set_ylabel(spec.ylabel)
        if spec.invert_x:
            self.ax.set_xlim(float(spec.x.max()), float(spec.x.min()))
        for side in ("top", "right"):
            self.ax.spines[side].set_visible(False)


class IRTab(_SpectrumTab):
    def __init__(self, vibrations, viewport, parent=None):
        super().__init__(parent)
        self.vibrations = vibrations
        self.viewport = viewport
        self._active_index: int | None = None

        self.fwhm = QDoubleSpinBox()
        self.fwhm.setRange(1.0, 200.0)
        self.fwhm.setValue(8.0)
        self.fwhm.setSuffix(" cm⁻¹")
        self.scale = QDoubleSpinBox()
        self.scale.setRange(0.5, 1.5)
        self.scale.setDecimals(4)
        self.scale.setSingleStep(0.001)
        self.scale.setValue(1.0)
        self._add_control("FWHM", self.fwhm)
        self._add_control("scale", self.scale)
        self._finish_controls("click a band to animate the mode")
        self.fwhm.valueChanged.connect(self.redraw)
        self.scale.valueChanged.connect(self.redraw)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.redraw()

    def redraw(self):
        try:
            self._spec = ir_spectrum(self.vibrations, scale=self.scale.value(),
                                     fwhm=self.fwhm.value())
        except ValueError:
            return
        spec = self._spec
        self.ax.clear()
        top = float(spec.y.max()) or 1.0
        self.ax.vlines(spec.stick_x, 0.0, spec.stick_y, color=_STICK, lw=1.0)
        self.ax.plot(spec.x, spec.y, color=_CURVE, lw=1.3)
        if self._active_index is not None:
            k = self._active_index
            self.ax.vlines([spec.stick_x[k]], 0.0, [max(spec.stick_y[k], 0.05 * top)],
                           color=_ACTIVE, lw=2.2)
            self.ax.annotate(f"{spec.stick_x[k]:.0f} cm⁻¹",
                             (spec.stick_x[k], max(spec.stick_y[k], 0.05 * top)),
                             textcoords="offset points", xytext=(4, 4),
                             color=_ACTIVE, fontsize=9)
        self._style_axes(spec)
        self.canvas.draw_idle()

    def _on_click(self, event):
        if event.inaxes != self.ax or self._spec is None or event.xdata is None:
            return
        sx = self._spec.stick_x
        if len(sx) == 0:
            return
        k = int(np.argmin(np.abs(sx - event.xdata)))
        span = abs(float(self._spec.x.max()) - float(self._spec.x.min()))
        if abs(sx[k] - event.xdata) > 0.03 * span:
            if self.viewport.is_animating:
                self.viewport.stop_animation()
                self._active_index = None
                self.redraw()
            return
        if k == self._active_index and self.viewport.is_animating:
            self.viewport.stop_animation()
            self._active_index = None
        else:
            mode = self._spec.meta["modes"][k]
            if self.viewport.animate_mode(mode.displacements):
                self._active_index = k
            else:
                self._active_index = None
        self.redraw()


class UVVisTab(_SpectrumTab):
    def __init__(self, states, parent=None):
        super().__init__(parent)
        self.states = states
        self.ax2 = self.ax.twinx()
        self.fwhm = QDoubleSpinBox()
        self.fwhm.setRange(0.02, 2.0)
        self.fwhm.setSingleStep(0.05)
        self.fwhm.setValue(0.30)
        self.fwhm.setSuffix(" eV")
        self.unit = QComboBox()
        self.unit.addItems(["nm", "eV"])
        self._add_control("FWHM", self.fwhm)
        self._add_control("axis", self.unit)
        self._finish_controls()
        self.fwhm.valueChanged.connect(self.redraw)
        self.unit.currentTextChanged.connect(self.redraw)
        self.redraw()

    def redraw(self):
        try:
            self._spec = uvvis_spectrum(self.states, fwhm_ev=self.fwhm.value(),
                                        unit=self.unit.currentText())
        except ValueError:
            return
        spec = self._spec
        self.ax.clear()
        self.ax2.clear()
        self.ax2.vlines(spec.stick_x, 0.0, spec.stick_y, color=_STICK, lw=1.0)
        self.ax2.set_ylabel("oscillator strength f", color=_STICK)
        self.ax2.tick_params(axis="y", labelcolor=_STICK)
        self.ax2.set_ylim(bottom=0.0)
        self.ax.plot(spec.x, spec.y, color=_CURVE, lw=1.3)
        self.ax.set_ylim(bottom=0.0)
        self._style_axes(spec)
        self.canvas.draw_idle()


class NMRTab(_SpectrumTab):
    def __init__(self, shieldings, parent=None):
        super().__init__(parent)
        self.shieldings = shieldings
        elems = sorted({s.symbol for s in shieldings},
                       key=lambda e: {"H": 0, "C": 1}.get(e, 2))
        self.element = QComboBox()
        self.element.addItems(elems)
        self.reference = QDoubleSpinBox()
        self.reference.setRange(0.0, 10000.0)
        self.reference.setDecimals(2)
        self.reference.setSuffix(" ppm")
        self.reference.setToolTip(
            "Reference shielding σ₀ (e.g. TMS computed at the same level).\n"
            "0 = plot raw isotropic shieldings.")
        self.fwhm = QDoubleSpinBox()
        self.fwhm.setRange(0.01, 20.0)
        self.fwhm.setSingleStep(0.1)
        self.fwhm.setValue(0.5)
        self.fwhm.setSuffix(" ppm")
        self._add_control("nucleus", self.element)
        self._add_control("σ₀ ref", self.reference)
        self._add_control("FWHM", self.fwhm)
        self._finish_controls("σ₀ = 0 plots raw shielding")
        for w in (self.reference, self.fwhm):
            w.valueChanged.connect(self.redraw)
        self.element.currentTextChanged.connect(self.redraw)
        self.redraw()

    def redraw(self):
        try:
            self._spec = nmr_spectrum(self.shieldings,
                                      element=self.element.currentText(),
                                      reference=self.reference.value(),
                                      fwhm=self.fwhm.value())
        except ValueError:
            return
        spec = self._spec
        self.ax.clear()
        self.ax.vlines(spec.stick_x, 0.0, spec.stick_y, color=_STICK, lw=1.0)
        self.ax.plot(spec.x, spec.y, color=_CURVE, lw=1.3)
        self.ax.set_ylim(bottom=0.0)
        self._style_axes(spec)
        self.canvas.draw_idle()


class SpectraDock(QDockWidget):
    def __init__(self, viewport, parent=None):
        super().__init__("Spectra", parent)
        self.viewport = viewport
        self.setObjectName("spectraDock")
        self.tabs = QTabWidget()
        self.setWidget(self.tabs)
        self.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea
                             | Qt.DockWidgetArea.RightDockWidgetArea)

    def set_result(self, result: ParseResult | None) -> None:
        self.viewport.stop_animation()
        while self.tabs.count():
            w = self.tabs.widget(0)
            self.tabs.removeTab(0)
            w.deleteLater()
        if result is None:
            self.hide()
            return
        if result.vibrations:
            self.tabs.addTab(IRTab(result.vibrations, self.viewport), "IR")
        if result.excited_states:
            self.tabs.addTab(UVVisTab(result.excited_states), "UV-Vis")
        if result.nmr_shieldings:
            self.tabs.addTab(NMRTab(result.nmr_shieldings), "NMR")
        self.setVisible(self.tabs.count() > 0)

    @property
    def has_data(self) -> bool:
        return self.tabs.count() > 0
