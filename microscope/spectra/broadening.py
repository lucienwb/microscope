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
    if shape not in ("gaussian", "lorentzian"):
        raise ValueError(f"unknown line shape {shape!r}")
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    hw = 0.5 * fwhm

    # The answer is a weighted sum over peaks, so the (grid x peaks) matrix is
    # never needed whole: a few thousand normal modes on a few thousand grid
    # points would be hundreds of megabytes of it. Accumulate in blocks.
    out = np.zeros_like(grid)
    block = max(1, 2_000_000 // max(grid.size, 1))
    x = grid[:, None]
    for start in range(0, centers.size, block):
        c = centers[start:start + block]
        if shape == "gaussian":
            peaks = np.exp(-((x - c) ** 2) / (2.0 * sigma * sigma))
        else:
            peaks = hw * hw / ((x - c) ** 2 + hw * hw)
        out += peaks @ heights[start:start + block]
    return out
