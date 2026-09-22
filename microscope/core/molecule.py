"""Molecule data model."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import elements
from .neighbors import close_pairs


@dataclass
class Molecule:
    symbols: list[str]
    coords: np.ndarray                 # (N, 3) float64, Angstrom
    charge: int = 0
    multiplicity: int = 1
    title: str = ""
    bonds: np.ndarray | None = None    # (M, 2) int, 0-based atom indices
    charge_known: bool = True          # False when the format never states one
    # atomic numbers, and the symbols list they were worked out from
    _numbers: tuple = field(default=(None, None), init=False, repr=False, compare=False)

    def __post_init__(self):
        self.coords = np.asarray(self.coords, dtype=np.float64).reshape(-1, 3)
        self.symbols = [elements.normalize_symbol(s) for s in self.symbols]
        if len(self.symbols) != len(self.coords):
            raise ValueError(
                f"{len(self.symbols)} symbols but {len(self.coords)} coordinates"
            )

    @property
    def natoms(self) -> int:
        return len(self.symbols)

    @property
    def atomic_numbers(self) -> np.ndarray:
        """int16 atomic numbers - the periodic table stops well short of 32767.
        Kept until ``symbols`` is replaced (nothing edits it in place), since
        rendering and perception read it several times over; on a protein
        each rebuild was 13 ms."""
        source, numbers = self._numbers
        if source is not self.symbols or len(numbers) != len(self.symbols):
            table = elements._SYMBOL_TO_Z
            numbers = np.array([table.get(s) or elements.symbol_to_z(s)
                                for s in self.symbols], dtype=np.int16)
            numbers.flags.writeable = False           # shared: nobody may edit it
            self._numbers = (self.symbols, numbers)
        return numbers

    def perceive_bonds(self, tolerance: float = 0.45) -> np.ndarray:
        """Detect bonds: d(i,j) <= r_cov(i) + r_cov(j) + tolerance (Angstrom).

        Every pair closer than the longest possible bond is tested, but they
        are found through a grid of cells one bond wide (core.neighbors) rather
        than by comparing all N^2 pairs: a protein has tens of thousands of
        atoms, and the pairwise difference array alone would be hundreds of
        gigabytes.
        """
        radii = elements.covalent_radii(self.atomic_numbers)
        reach = 2.0 * float(radii.max(initial=0.0)) + tolerance
        i, j, d = close_pairs(self.coords, reach)
        keep = (d <= radii[i] + radii[j] + tolerance) & (d > 0.4)
        self.bonds = np.column_stack([i[keep], j[keep]]).astype(np.int32)
        return self.bonds

    def formula(self) -> str:
        """Molecular formula in Hill order (C, H, then alphabetical)."""
        counts = Counter(self.symbols)
        parts = []
        for sym in ("C", "H"):
            if sym in counts:
                c = counts.pop(sym)
                parts.append(sym + (str(c) if c > 1 else ""))
        for sym in sorted(counts):
            c = counts[sym]
            parts.append(sym + (str(c) if c > 1 else ""))
        return "".join(parts)

    def bounding_sphere(self) -> tuple[np.ndarray, float]:
        if self.natoms == 0:
            return np.zeros(3), 1.0
        lo, hi = self.coords.min(axis=0), self.coords.max(axis=0)
        center = (lo + hi) / 2.0
        radius = float(np.linalg.norm(self.coords - center, axis=1).max())
        return center, max(radius, 1.0)

    def copy(self) -> Molecule:
        return Molecule(
            symbols=list(self.symbols),
            coords=self.coords.copy(),
            charge=self.charge,
            multiplicity=self.multiplicity,
            title=self.title,
            bonds=None if self.bonds is None else self.bonds.copy(),
        )
