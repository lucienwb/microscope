"""The orbital energy-level diagram, laid out before anything is drawn.

One short line per orbital at its energy, degenerate orbitals side by side,
electrons as arrows, alpha and beta in two columns when the wavefunction is
unrestricted. Pure numpy, like lewis2d: levelsdraw turns it into lines on any
QPainter device, the orbital dialog and an exported SVG alike.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.orbitals import HARTREE_TO_EV, Orbitals

DEGENERATE_EV = 0.02        # closer than this, two orbitals are degenerate
WINDOW = 5                  # orbitals shown beyond the frontier or the selection


@dataclass
class Level:
    spin: str               # "alpha" (or restricted) / "beta"
    index: int              # 0-based, as Orbitals counts them
    energy: float           # eV
    column: int             # 0 restricted or alpha, 1 beta
    degenerate: int         # how many orbitals share this energy, itself included
    electrons: int          # arrows to draw: 0, 1 or 2
    name: str               # HOMO, LUMO+1, ...


@dataclass
class LevelDiagram:
    levels: list[Level] = field(default_factory=list)
    gaps: dict[str, float] = field(default_factory=dict)   # spin -> LUMO - HOMO, eV
    columns: int = 1

    @property
    def span(self) -> tuple[float, float]:
        energies = [lv.energy for lv in self.levels]
        return (min(energies), max(energies)) if energies else (0.0, 1.0)

    def find(self, spin: str, index: int) -> Level | None:
        return next((lv for lv in self.levels if lv.spin == spin and lv.index == index),
                    None)


def level_diagram(orbitals: Orbitals, around: tuple[str, int] | None = None,
                  window: int = WINDOW) -> LevelDiagram:
    """The levels from *window* below the HOMO to *window* above the LUMO, and
    around the orbital *around* too when it lies further out."""
    spins = ("alpha",) if orbitals.restricted else ("alpha", "beta")
    diagram = LevelDiagram(columns=len(spins))
    for column, spin in enumerate(spins):
        mos = orbitals.spin(spin)
        energies = np.asarray(mos.energies, dtype=float) * HARTREE_TO_EV
        homo = mos.homo
        low, high = homo - window + 1, homo + window          # as many each side
        if around is not None and around[0] == spin:
            low, high = min(low, around[1] - 2), max(high, around[1] + 2)
        shown = [i for i in range(max(low, 0), min(high, mos.nmo - 1) + 1)
                 if np.isfinite(energies[i])]
        full = 2 if orbitals.restricted else 1
        for i in shown:
            degenerate = int(np.sum(np.abs(energies[shown] - energies[i]) < DEGENERATE_EV))
            electrons = int(round(min(max(float(mos.occupations[i]), 0.0), full)))
            diagram.levels.append(Level(spin, i, float(energies[i]), column, degenerate,
                                        electrons, orbitals.name(spin, i)))
        if 0 <= homo < mos.nmo - 1 and homo in shown and homo + 1 in shown:
            diagram.gaps[spin] = float(energies[homo + 1] - energies[homo])
    return diagram


def side_by_side(heights: list[float], apart: float) -> list[tuple[int, int]]:
    """(slot, partners) for levels drawn at *heights* (pixels): a level closer
    than *apart* to the one below it joins it, side by side. Degenerate
    orbitals do, and so do ones too close in energy to be told apart at the
    scale drawn - which only the drawing knows."""
    order = sorted(range(len(heights)), key=lambda k: heights[k])
    groups: list[list[int]] = []
    for k in order:
        if groups and heights[k] - heights[groups[-1][-1]] < apart:
            groups[-1].append(k)
        else:
            groups.append([k])
    placed = [(0, 1)] * len(heights)
    for group in groups:
        # left to right in rising energy: the lowest level on screen is last
        for slot, k in enumerate(reversed(group)):
            placed[k] = (slot, len(group))
    return placed
