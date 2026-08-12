"""Line-shape broadening of stick spectra."""

from __future__ import annotations

import numpy as np


def broaden(centers, heights, grid, fwhm: float, shape: str = "lorentzian") -> np.ndarray:
    """Sum of unit-peak-height line shapes at *centers* scaled by *heights*.

    Peaks are evaluated on *grid*; each contributes `height` at its center.
    """
    centers = np.asarray(centers, dtype=float).ravel()
    heights = np.asarray(heights, dtype=float).ravel()
    grid = np.asarray(grid, dtype=float)
    if centers.size == 0:
        return np.zeros_like(grid)
    x = grid[:, None]
    if shape == "gaussian":
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        peaks = np.exp(-((x - centers) ** 2) / (2.0 * sigma * sigma))
    elif shape == "lorentzian":
        hw = 0.5 * fwhm
        peaks = hw * hw / ((x - centers) ** 2 + hw * hw)
    else:
        raise ValueError(f"unknown line shape {shape!r}")
    return peaks @ heights
