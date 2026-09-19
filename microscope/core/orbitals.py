"""Molecular orbitals read from a wavefunction file, sampled on a grid to draw.

The orbitals are the program's own: its basis set and its coefficients, read
from an fchk or Molden file. Putting one on a grid is what cubegen or
orca_plot do, done here so that the file draws its orbitals without a cube
file being made first. The grid comes out as a :class:`VolumeData`, so it goes
through the same isosurface path a cube file does.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from .basis import (
    BOHR_TO_ANGSTROM,
    BasisSet,
    cartesian_components,
    orthonormality_error,
)
from .volume import VolumeData

HARTREE_TO_EV = 27.211386245988
SPACING = 0.12                 # Angstrom between grid points
PADDING = 4.0                  # Angstrom of space around the outermost atoms, to start with
MAX_PADDING = 12.0             # a diffuse orbital gets more room, up to this
EDGE = 0.005                   # the box grows while the orbital is above this at its walls
MAX_POINTS = 8_000_000         # a bigger molecule gets a coarser grid, not a slower one
CUTOFF = 1e-6                  # a shell is left out where it adds less than this

SPINS = ("alpha", "beta")
ORTHONORMAL = 1e-3             # a file's orbitals are orthonormal to this, or it was misread
SAMPLE = 600                   # orbitals checked at most; a spread of them tells as well as all


@dataclass
class OrbitalSet:
    """The orbitals of one spin (or of both, for a restricted wavefunction)."""

    energies: np.ndarray           # (nmo,) hartree
    occupations: np.ndarray        # (nmo,) electrons
    coefficients: np.ndarray       # (nmo, nbf), float32: one row per orbital
    symmetries: list[str] = field(default_factory=list)

    @property
    def nmo(self) -> int:
        return len(self.coefficients)

    @property
    def homo(self) -> int:
        """Index of the highest occupied orbital, or -1 if nothing is."""
        occupied = np.flatnonzero(np.asarray(self.occupations) >= 0.5)
        return int(occupied[-1]) if len(occupied) else -1

    def sample(self, limit: int = SAMPLE) -> np.ndarray:
        """Rows to check a reading with: all of them for a small file, else an
        even spread from the lowest core orbital to the highest virtual, the
        HOMO always among them."""
        if self.nmo <= limit:
            return self.coefficients
        rows = np.unique(np.concatenate([np.linspace(0, self.nmo - 1, limit - 1).astype(int),
                                         [max(self.homo, 0)]]))
        return self.coefficients[rows]


@dataclass
class Orbitals:
    """A wavefunction: the basis set, and the orbitals expanded in it.

    Attributes:
        basis: The basis set, with every convention of the file folded in.
        alpha: The orbitals of a restricted wavefunction, or the alpha ones.
        beta: The beta orbitals of an unrestricted wavefunction, else None.
        convention: How the file had to be read: "standard", a named
            variant such as ORCA's, or "unrecognized" when no reading made
            its orbitals orthonormal - drawn anyway, but not to be trusted.
    """

    basis: BasisSet
    alpha: OrbitalSet
    beta: OrbitalSet | None = None
    convention: str = "standard"      # how the file's basis had to be read

    @property
    def restricted(self) -> bool:
        return self.beta is None

    def spin(self, name: str = "alpha") -> OrbitalSet:
        if name == "beta":
            if self.beta is None:
                raise ValueError("a restricted wavefunction has no separate beta orbitals")
            return self.beta
        return self.alpha

    def find(self, spec: str) -> tuple[str, int]:
        """Which orbital a name means, as (spin, 0-based index).

        ``homo``, ``lumo``, ``homo-2``, ``lumo+1``, or a 1-based number as the
        program counts them, optionally after ``alpha:`` or ``beta:``.
        """
        text = str(spec).strip().lower().replace(" ", "")
        spin = "alpha"
        if ":" in text:
            spin, text = text.split(":", 1)
            spin = {"a": "alpha", "b": "beta"}.get(spin, spin)
            if spin not in SPINS:
                raise ValueError(f"{spec!r}: the spin is alpha or beta, not {spin!r}")
        orbitals = self.spin(spin)
        match = re.fullmatch(r"(homo|lumo)([+-]\d+)?", text)
        if match:
            if orbitals.homo < 0:
                raise ValueError(f"{spec!r}: the file gives no occupations, so there is no "
                                 "HOMO to count from; give the orbital's number instead")
            index = orbitals.homo + (match.group(1) == "lumo") + int(match.group(2) or 0)
        elif text.isdigit():
            index = int(text) - 1
        else:
            raise ValueError(f"{spec!r} is not an orbital: use homo, lumo, homo-1, "
                             "lumo+2 or a number, optionally after alpha: or beta:")
        if not 0 <= index < orbitals.nmo:
            raise ValueError(f"{spec!r} is outside the {orbitals.nmo} {spin} orbitals"
                             if not self.restricted else
                             f"{spec!r} is outside the {orbitals.nmo} orbitals")
        return spin, index

    def name(self, spin: str, index: int) -> str:
        """HOMO, HOMO-3, LUMO+1: where an orbital sits relative to the frontier."""
        homo = self.spin(spin).homo
        if index <= homo:
            return "HOMO" if index == homo else f"HOMO-{homo - index}"
        return "LUMO" if index == homo + 1 else f"LUMO+{index - homo - 1}"

    def describe(self, spin: str, index: int) -> str:
        """A one-line label: 'HOMO-1 · MO 20 · -9.52 eV', with the spin if it matters."""
        orbitals = self.spin(spin)
        parts = [self.name(spin, index), f"MO {index + 1}"]
        energy = float(orbitals.energies[index])
        if np.isfinite(energy):
            parts.append(f"{energy * HARTREE_TO_EV:.2f} eV")
        text = " · ".join(parts)
        return text if self.restricted else f"{'α' if spin == 'alpha' else 'β'} {text}"

    def volume(self, which: str | tuple[str, int] = "homo", *, spacing: float = SPACING,
               padding: float = PADDING, max_points: int = MAX_POINTS) -> VolumeData:
        """One orbital sampled on a grid around the molecule, ready to draw.

        The box starts *padding* past the outermost atoms and grows while the
        orbital is still above EDGE at its walls, so a diffuse virtual is not
        drawn with its lobes sliced flat - but a bigger box is a coarser grid,
        never more than *max_points*.
        """
        spin, index = self.find(which) if isinstance(which, str) else which
        coefficients = self.spin(spin).coefficients[index]
        while True:
            origin, step, shape = grid_around(self.centers(), spacing, padding, max_points)
            values = evaluate(self.basis, coefficients, origin, step, shape)
            if padding >= MAX_PADDING or _at_walls(values) < EDGE:
                break
            padding = min(padding * 1.5, MAX_PADDING)
        return VolumeData(origin=origin * BOHR_TO_ANGSTROM,
                          axes=np.diag(step) * BOHR_TO_ANGSTROM, values=values,
                          label=self.describe(spin, index),
                          meta={"spin": spin, "orbital": index})

    def reading_error(self) -> float:
        """How far a sample of the orbitals is from orthonormal, as read."""
        return max(orthonormality_error(self.basis.shells, s.sample())
                   for s in (self.alpha, self.beta) if s is not None)

    def centers(self) -> np.ndarray:
        """Where the basis functions sit (bohr): the atoms, ghosts included."""
        return np.array([s.center for s in self.basis.shells]).reshape(-1, 3)

    def orthonormality_error(self, occupied_only: bool = True) -> float:
        """Largest deviation of C^T S C from the identity over the orbitals.

        Near zero when the file's basis was read the way the program meant it;
        the occupied orbitals are the ones a program always writes carefully.
        """
        worst = 0.0
        for spin in SPINS if not self.restricted else ("alpha",):
            orbitals = self.spin(spin)
            C = orbitals.coefficients
            if occupied_only:
                C = C[:max(orbitals.homo + 1, 1)]
            worst = max(worst, orthonormality_error(self.basis.shells, C))
        return worst


# --------------------------------------------------------------------- grid

def grid_around(centers: np.ndarray, spacing: float = SPACING, padding: float = PADDING,
                max_points: int = MAX_POINTS):
    """(origin, step, shape) of a grid boxing the centres, in bohr."""
    pad = padding / BOHR_TO_ANGSTROM
    lo = centers.min(axis=0) - pad
    hi = centers.max(axis=0) + pad
    extent = hi - lo
    h = spacing / BOHR_TO_ANGSTROM
    h = max(h, float(np.prod(extent) / max_points) ** (1 / 3))
    shape = tuple(int(n) for n in np.ceil(extent / h).astype(int) + 1)
    step = np.full(3, h)
    # centre the grid on the box, so rounding the count up costs both sides equally
    origin = lo - (step * (np.array(shape) - 1) - extent) / 2
    return origin, step, shape


def _at_walls(values: np.ndarray) -> float:
    """The largest magnitude on the six faces of a grid."""
    return float(max(np.abs(face).max() for face in (
        values[0], values[-1], values[:, 0], values[:, -1], values[:, :, 0], values[:, :, -1])))


def _reach(shell, amplitude: float, cutoff: float) -> float:
    """How far (bohr) the shell stays above the cutoff for coefficients this big."""
    weights = np.abs(shell.weights) * amplitude * len(shell.weights) / cutoff
    alive = weights > 1.0
    if not alive.any():
        return 0.0
    k = np.log(weights[alive])
    alpha = shell.exponents[alive]
    r = np.sqrt(k / alpha) + np.sqrt(shell.l / (2 * alpha))
    for _ in range(4):          # r^l exp(-alpha r^2) = 1/weight, outer root
        r = np.sqrt(np.maximum(k + shell.l * np.log(np.maximum(r, 1e-3)), 0.0) / alpha)
    return float(r.max())


def evaluate(basis: BasisSet, coefficients: np.ndarray, origin: np.ndarray,
             step: np.ndarray, shape: tuple[int, int, int],
             cutoff: float = CUTOFF) -> np.ndarray:
    """An orbital's values on an axis-aligned grid (bohr), float32.

    Shell by shell: a Gaussian on such a grid factorizes into one factor per
    axis, so each shell is built from 1-D arrays and joined with one matrix
    product, over only the sub-box it reaches for this orbital's coefficients.
    """
    values = np.zeros(shape, dtype=np.float32)
    coefficients = np.asarray(coefficients, dtype=float)
    axes = [origin[k] + step[k] * np.arange(shape[k]) for k in range(3)]
    for shell, start in zip(basis.shells, basis.offsets()):
        c = coefficients[start:start + shell.nfunc] @ shell.transform    # (ncart,)
        amplitude = float(np.abs(c).sum())
        if amplitude == 0.0:
            continue
        reach = _reach(shell, amplitude, cutoff)
        if reach == 0.0:
            continue
        box = []
        for k in range(3):
            lo = int(np.ceil((shell.center[k] - reach - origin[k]) / step[k]))
            hi = int(np.floor((shell.center[k] + reach - origin[k]) / step[k])) + 1
            box.append((max(lo, 0), min(hi, shape[k])))
        if any(lo >= hi for lo, hi in box):
            continue
        l = shell.l
        powers = np.arange(l + 1)
        factors = []
        for k in range(3):
            d = axes[k][box[k][0]:box[k][1]] - shell.center[k]              # (n,)
            gauss = np.exp(-shell.exponents[:, None] * d[None, :] ** 2)     # (p, n)
            factors.append(gauss[:, None, :] * d[None, None, :] ** powers[None, :, None])
        fx, fy, fz = factors                                                # (p, l+1, n)
        fx = fx * shell.weights[:, None, None]
        tensor = np.zeros((l + 1, l + 1, l + 1))
        comps = cartesian_components(l)
        tensor[comps[:, 0], comps[:, 1], comps[:, 2]] = c
        # sum_abc T[a,b,c] X_pa Y_pb Z_pc, contracted one axis at a time
        yz = np.einsum("abc,pbj,pck->pajk", tensor, fy, fz, optimize=True)
        p, _, ny, nz = yz.shape
        block = fx.reshape(p * (l + 1), -1).T.astype(np.float32) @ \
            yz.reshape(p * (l + 1), ny * nz).astype(np.float32)
        (x0, x1), (y0, y1), (z0, z1) = box
        values[x0:x1, y0:y1, z0:z1] += block.reshape(x1 - x0, ny, nz)
    return values
