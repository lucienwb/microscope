"""Spectrum construction from parsed calculation results."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.results import ExcitedState, NMRShielding, Vibration
from .broadening import broaden

EV_NM = 1239.841984


@dataclass
class Spectrum:
    x: np.ndarray
    y: np.ndarray
    stick_x: np.ndarray
    stick_y: np.ndarray
    xlabel: str
    ylabel: str
    invert_x: bool = False
    meta: dict = field(default_factory=dict)


def ir_spectrum(vibrations: list[Vibration], scale: float = 1.0,
                fwhm: float = 8.0, npoints: int = 4000) -> Spectrum:
    """Broadened IR spectrum; imaginary modes are excluded from the plot."""
    real = [v for v in vibrations if v.frequency > 0]
    if not real:
        raise ValueError("no real vibrational frequencies")
    freqs = np.array([v.frequency * scale for v in real])
    heights = np.array([v.ir_intensity if v.ir_intensity is not None else 0.0
                        for v in real])
    xmax = max(float(freqs.max()) * 1.08, 500.0)
    grid = np.linspace(0.0, xmax, npoints)
    y = broaden(freqs, heights, grid, fwhm, "lorentzian")
    return Spectrum(grid, y, freqs, heights,
                    "wavenumber (cm⁻¹)", "IR intensity (km/mol)",
                    invert_x=True, meta={"modes": real})


def uvvis_spectrum(states: list[ExcitedState], fwhm_ev: float = 0.30,
                   unit: str = "nm", npoints: int = 2000) -> Spectrum:
    """Broadened UV-Vis spectrum (Gaussian in the energy domain)."""
    if not states:
        raise ValueError("no excited states")
    ev = np.array([s.energy_ev for s in states])
    f = np.array([s.osc_strength for s in states])
    emin = max(float(ev.min()) - 3.0 * fwhm_ev, 0.1)
    emax = float(ev.max()) + 3.0 * fwhm_ev
    grid_ev = np.linspace(emin, emax, npoints)
    y = broaden(ev, f, grid_ev, fwhm_ev, "gaussian")
    if unit == "nm":
        x = EV_NM / grid_ev
        order = np.argsort(x)
        return Spectrum(x[order], y[order], EV_NM / ev, f,
                        "wavelength (nm)", "intensity (arb. units)",
                        meta={"states": states})
    return Spectrum(grid_ev, y, ev, f, "energy (eV)", "intensity (arb. units)",
                    meta={"states": states})


def nmr_spectrum(shieldings: list[NMRShielding], element: str = "H",
                 reference: float = 0.0, fwhm: float = 0.5,
                 npoints: int = 4000) -> Spectrum:
    """NMR spectrum for one element.

    With reference == 0 the raw isotropic shieldings are plotted; a non-zero
    reference sigma_0 switches to the chemical-shift scale delta = sigma_0 - sigma.
    """
    sel = [s for s in shieldings if s.symbol == element]
    if not sel:
        raise ValueError(f"no {element} shieldings available")
    sigma = np.array([s.isotropic for s in sel])
    if reference:
        sticks = reference - sigma
        xlabel = "chemical shift δ (ppm)"
        invert = True
    else:
        sticks = sigma
        xlabel = "isotropic shielding σ (ppm)"
        invert = False
    heights = np.ones(len(sel))
    span = float(np.ptp(sticks))
    pad = max(4.0 * fwhm, 0.08 * span, 1.0)
    grid = np.linspace(float(sticks.min()) - pad, float(sticks.max()) + pad, npoints)
    y = broaden(sticks, heights, grid, fwhm, "lorentzian")
    return Spectrum(grid, y, sticks, heights, xlabel, "intensity (arb. units)",
                    invert_x=invert, meta={"nuclei": sel, "element": element})
