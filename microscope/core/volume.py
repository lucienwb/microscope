"""Volumetric scalar data on a regular grid (cube files, future MO grids)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class VolumeData:
    """Scalar field sampled on a (possibly non-orthogonal) regular grid.

    ``origin`` and the three ``axes`` step vectors are in Angstrom; grid point
    (i, j, k) sits at ``origin + i*axes[0] + j*axes[1] + k*axes[2]``.
    """

    origin: np.ndarray                     # (3,)
    axes: np.ndarray                       # (3, 3), rows = step vectors
    values: np.ndarray                     # (n1, n2, n3)
    label: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.values.shape

    @property
    def is_signed(self) -> bool:
        """True for orbital-like data with +/- lobes (vs. densities >= 0)."""
        return float(self.values.min()) < 0.0

    def grid_to_world(self, points: np.ndarray) -> np.ndarray:
        """Map fractional grid indices (N, 3) to world coordinates (N, 3)."""
        return np.asarray(points, dtype=float) @ self.axes + self.origin

    def suggest_isovalue(self) -> float:
        """A sensible default isovalue for display.

        Orbitals are conventionally shown at |psi| = 0.02 au and densities at
        rho = 0.002 au; fall back to a value inside the data range if the
        conventional one would produce an empty surface.
        """
        amax = float(np.abs(self.values).max())
        if amax == 0.0:
            return 0.02
        default = 0.02 if self.is_signed else 0.002
        if default < amax:
            return default
        return 0.1 * amax
