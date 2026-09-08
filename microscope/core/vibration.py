"""A normal mode oscillating about a geometry.

Clicking an IR band sets the molecule walking through that vibration. The
walking is arithmetic — a displacement scaled by a sine — and lives here so
it can be checked without a timer or a widget; the viewport keeps the timer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

STEP = 0.3            # radians of phase per frame, at the viewport's 33 ms
AMPLITUDE = 0.35      # Angstrom at the atom that moves most


@dataclass
class ModeAnimation:
    """Coordinates swinging along one normal mode."""

    base: np.ndarray                   # the geometry it returns to
    displacement: np.ndarray           # already scaled to the amplitude
    phase: float = 0.0

    @classmethod
    def of(cls, coords: np.ndarray, displacements,
           amplitude: float = AMPLITUDE) -> ModeAnimation | None:
        """None when this mode cannot be animated: wrong shape, or no motion."""
        if displacements is None:
            return None
        d = np.asarray(displacements, dtype=float)
        if d.shape != coords.shape:
            return None
        peak = float(np.abs(d).max())
        if peak < 1e-9:
            return None
        return cls(base=coords.copy(), displacement=d / peak * amplitude)

    def at(self, phase: float) -> np.ndarray:
        return self.base + np.sin(phase) * self.displacement

    def step(self, delta: float = STEP) -> np.ndarray:
        """Advance one frame and return where the atoms are now."""
        self.phase += delta
        return self.at(self.phase)
