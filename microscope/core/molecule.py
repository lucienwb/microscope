"""Molecule data model."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

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
    charge_known: bool = True          # False when the format never states one

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
        # int16: the periodic table stops well short of 32767
        return np.array([elements.symbol_to_z(s) for s in self.symbols],
                        dtype=np.int16)

    def perceive_bonds(self, tolerance: float = 0.45) -> np.ndarray:
        """Detect bonds: d(i,j) <= r_cov(i) + r_cov(j) + tolerance (Angstrom).

        Every pair closer than the longest possible bond is tested, but they
        are found through a grid of cells one bond wide rather than by
        comparing all N^2 pairs: a protein has tens of thousands of atoms, and
        the pairwise difference array alone would be hundreds of gigabytes.
        """
        n = self.natoms
        if n < 2:
            self.bonds = np.zeros((0, 2), dtype=np.int32)
            return self.bonds
        radii = np.array([elements.covalent_radius(z) for z in self.atomic_numbers])
        reach = 2.0 * float(radii.max()) + tolerance

        coords = self.coords
        cells = np.floor((coords - coords.min(axis=0)) / reach).astype(np.int64)
        buckets: dict[tuple[int, int, int], list[int]] = {}
        for index, cell in enumerate(map(tuple, cells)):
            buckets.setdefault(cell, []).append(index)

        offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                   for dz in (-1, 0, 1)]
        pairs: list[np.ndarray] = []
        for (cx, cy, cz), members in buckets.items():
            near = [a for dx, dy, dz in offsets
                    for a in buckets.get((cx + dx, cy + dy, cz + dz), ())]
            here = np.array(members)
            there = np.array(near)
            dist = np.linalg.norm(coords[here][:, None, :] - coords[there][None, :, :],
                                  axis=2)
            cutoff = radii[here][:, None] + radii[there][None, :] + tolerance
            # i < j keeps each pair once, whichever cell it is reached from
            keep = (dist <= cutoff) & (dist > 0.4) & (here[:, None] < there[None, :])
            rows, cols = np.where(keep)
            if len(rows):
                pairs.append(np.column_stack([here[rows], there[cols]]))

        if pairs:
            bonds = np.vstack(pairs)
            order = np.lexsort((bonds[:, 1], bonds[:, 0]))
            self.bonds = bonds[order].astype(np.int32)
        else:
            self.bonds = np.zeros((0, 2), dtype=np.int32)
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
