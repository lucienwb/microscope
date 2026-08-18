"""Molecule data model."""

from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter

import numpy as np

from . import elements


@dataclass
class Molecule:
    symbols: list[str]
    coords: np.ndarray                 # (N, 3) float64, Angstrom
    charge: int = 0
    multiplicity: int = 1
    title: str = ""
    bonds: np.ndarray | None = None    # (M, 2) int, 0-based atom indices

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
        return np.array([elements.symbol_to_z(s) for s in self.symbols], dtype=int)

    def perceive_bonds(self, tolerance: float = 0.45) -> np.ndarray:
        """Detect bonds: d(i,j) <= r_cov(i) + r_cov(j) + tolerance (Angstrom)."""
        n = self.natoms
        if n < 2:
            self.bonds = np.zeros((0, 2), dtype=int)
            return self.bonds
        radii = np.array([elements.covalent_radius(z) for z in self.atomic_numbers])
        diff = self.coords[:, None, :] - self.coords[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        cutoff = radii[:, None] + radii[None, :] + tolerance
        mask = (dist <= cutoff) & (dist > 0.4)
        ii, jj = np.where(np.triu(mask, k=1))
        self.bonds = np.column_stack([ii, jj]).astype(int)
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

    def copy(self) -> "Molecule":
        return Molecule(
            symbols=list(self.symbols),
            coords=self.coords.copy(),
            charge=self.charge,
            multiplicity=self.multiplicity,
            title=self.title,
            bonds=None if self.bonds is None else self.bonds.copy(),
        )
