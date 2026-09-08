"""Undo and redo for structure edits, as plain data.

The viewer is where an edit is triggered and where the result is drawn, but
what may be undone is only a question of which states have been visited. That
part is kept here: no Qt, no widgets, so it can be reasoned about and tested
on its own.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DEPTH = 100          # states kept; a molecule is small but not free


@dataclass(frozen=True)
class Snapshot:
    """Everything an edit can change about a molecule."""

    coords: np.ndarray
    symbols: tuple[str, ...]
    bonds: np.ndarray | None

    @classmethod
    def of(cls, molecule) -> Snapshot:
        return cls(coords=molecule.coords.copy(),
                   symbols=tuple(molecule.symbols),
                   bonds=None if molecule.bonds is None else molecule.bonds.copy())

    def fits(self, molecule) -> bool:
        """Can this be restored in place, or does the atom count differ?"""
        return molecule is not None and len(self.symbols) == molecule.natoms

    def restore_into(self, molecule) -> None:
        """Put this state back on *molecule*, which must be the right size."""
        molecule.coords = self.coords.copy()
        molecule.symbols = list(self.symbols)
        molecule.bonds = None if self.bonds is None else self.bonds.copy()


class EditHistory:
    """The states an edit can be taken back to, and forward to again.

    Push before each edit; ``undo`` and ``redo`` are handed the current state
    so they can put it on the other stack, and return the one to go to (or
    None when there is nowhere to go).
    """

    def __init__(self, depth: int = DEPTH):
        self.depth = depth
        self._undo: list[Snapshot] = []
        self._redo: list[Snapshot] = []

    def __len__(self) -> int:
        return len(self._undo)

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    def clear(self) -> None:
        self._undo.clear()
        self._redo.clear()

    def push(self, state: Snapshot) -> None:
        """Record where we were; anything undone from here is no longer reachable."""
        self._undo.append(state)
        del self._undo[:-self.depth]
        self._redo.clear()

    def undo(self, current: Snapshot) -> Snapshot | None:
        if not self._undo:
            return None
        self._redo.append(current)
        return self._undo.pop()

    def redo(self, current: Snapshot) -> Snapshot | None:
        if not self._redo:
            return None
        self._undo.append(current)
        return self._redo.pop()
