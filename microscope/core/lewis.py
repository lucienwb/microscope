"""Lewis-structure perception: bond orders, formal charges, lone pairs.

The viewer already infers *connectivity* from the geometry (see
``Molecule.perceive_bonds``); this module goes one step further and works out
how to *draw* that connectivity the way a chemist would — double and triple
lines, formal charges, lone-pair dots. It is drawing perception, not
chemistry: nothing is optimized or evaluated here, the geometry is read
exactly as the quantum-chemistry program wrote it.

Bond orders come from the bond length. Every element pair has well known
single/double/triple distances, so the closest reference wins; that guess is
then filtered through the valence each atom has left. Handing out the most
contracted bonds first is what makes a benzene ring come out as a Kekule
structure instead of six double bonds.

Perception cannot be right for everything a quantum chemist opens — a
carbocation and a carbanion have the same connectivity, and metal complexes
have no Lewis structure at all. So the result carries the perceived total
charge, and :attr:`LewisStructure.matches_file` says whether it agrees with
the charge the file itself declares; the viewer shows a warning when it does
not, rather than quietly presenting a wrong drawing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import elements
from .molecule import Molecule

# Reference bond lengths in Angstrom: (single, double, triple), None where the
# order does not occur. Keys are element symbol pairs in alphabetical order.
_REFERENCE: dict[tuple[str, str], tuple[float, float | None, float | None]] = {
    ("B", "C"): (1.56, 1.42, None),
    ("B", "N"): (1.55, 1.40, None),
    # boron-oxygen is single at every length that occurs: 1.36 in a trigonal
    # boronate ester or boroxine, 1.50 in a four-coordinate borate. A B=O
    # double bond is not a thing outside the gas phase, and listing one made
    # every pinacol boronate come out as a borate anion.
    ("B", "O"): (1.40, None, None),
    ("Br", "C"): (1.94, 1.83, None),
    ("C", "C"): (1.54, 1.34, 1.20),
    ("C", "Cl"): (1.77, 1.68, None),
    ("C", "F"): (1.35, 1.28, None),
    ("C", "I"): (2.14, 2.03, None),
    ("C", "N"): (1.47, 1.28, 1.16),
    # the triple is 1.128 in free CO, but back-donation stretches a metal
    # carbonyl to 1.13-1.18, and those are all drawn M-C(triple)O; 1.15 puts
    # the split with a real C=O (1.19 and up) between the two families
    ("C", "O"): (1.43, 1.21, 1.15),
    ("C", "P"): (1.85, 1.67, 1.54),
    ("C", "S"): (1.82, 1.60, 1.54),
    ("C", "Se"): (1.95, 1.73, None),
    ("C", "Si"): (1.87, 1.70, None),
    ("Cl", "O"): (1.70, 1.44, None),
    ("N", "N"): (1.45, 1.25, 1.10),
    ("N", "O"): (1.40, 1.21, 1.06),
    ("N", "P"): (1.72, 1.58, None),
    ("N", "S"): (1.71, 1.54, None),
    ("N", "Si"): (1.74, 1.57, None),
    ("O", "O"): (1.48, 1.21, None),
    ("O", "P"): (1.62, 1.48, None),
    ("O", "S"): (1.66, 1.44, None),
    ("O", "Se"): (1.78, 1.61, None),
    ("O", "Si"): (1.63, 1.51, None),
    ("P", "P"): (2.21, 2.03, 1.89),
    ("P", "S"): (2.10, 1.95, None),
    ("S", "S"): (2.05, 1.89, None),
}

# The most bonds an atom will be given. Period 3 and below may expand their
# octet (P 5, S 6, halogens 7), which is what draws sulfate and perchlorate.
_MAX_BONDS = {
    1: 1, 2: 0, 3: 1, 4: 2, 5: 4, 6: 4, 7: 4, 8: 3, 9: 1, 10: 0,
    11: 1, 12: 2, 13: 4, 14: 4, 15: 5, 16: 6, 17: 7, 18: 0,
    19: 1, 20: 2, 31: 4, 32: 4, 33: 5, 34: 6, 35: 7, 36: 0,
    37: 1, 38: 2, 49: 4, 50: 4, 51: 5, 52: 6, 53: 7, 54: 8,
    55: 1, 56: 2, 81: 4, 82: 4, 83: 5, 84: 6, 85: 7, 86: 8,
}

# Valence electrons by main group; everything else is left out of the Lewis
# bookkeeping altogether (transition metals have no formal charge here).
_GROUP_MEMBERS = {
    1: (1, 3, 11, 19, 37, 55, 87), 2: (4, 12, 20, 38, 56, 88),
    3: (5, 13, 31, 49, 81), 4: (6, 14, 32, 50, 82),
    5: (7, 15, 33, 51, 83), 6: (8, 16, 34, 52, 84),
    7: (9, 17, 35, 53, 85), 8: (2, 10, 18, 36, 54, 86),
}
VALENCE_ELECTRONS = {z: n for n, members in _GROUP_MEMBERS.items() for z in members}
VALENCE_ELECTRONS[2] = 2                      # helium fills its shell with two

# Electrons that count as a filled shell: a duet for H/He/Li, a sextet for
# boron and its electron-deficient neighbours, an octet for everything else.
_SHELL = {1: 2, 2: 2, 3: 2, 4: 4, 5: 6, 13: 6}

# Bond lengths within this fraction of each other count as the same length.
_TIE_WIDTH = 0.02


def reference_lengths(a: str, b: str):
    """Single/double/triple reference lengths for an element pair, or None."""
    return _REFERENCE.get(tuple(sorted((a, b))))


def is_main_group(z: int) -> bool:
    return int(z) in VALENCE_ELECTRONS


def max_bonds(z: int) -> int:
    return _MAX_BONDS.get(int(z), 0)


class Adjacency:
    """Who is bonded to whom, held as two arrays instead of a list per atom.

    Indexing gives a view of one atom's neighbours, so it reads like the list
    of lists it replaces. That list cost more than every other array in a
    structure put together: a protein has tens of thousands of atoms, and each
    one carried a Python list object holding boxed integers — 8 MB for
    photosystem II against 0.7 MB here.
    """

    __slots__ = ("offsets", "flat")

    def __init__(self, natoms: int, bonds: np.ndarray):
        ends = np.asarray(bonds, dtype=np.int32).reshape(-1)       # i0 j0 i1 j1
        others = np.asarray(bonds, dtype=np.int32)[:, ::-1].reshape(-1)
        order = np.argsort(ends, kind="stable")     # keeps bond order per atom
        self.flat = np.ascontiguousarray(others[order])
        self.offsets = np.zeros(natoms + 1, dtype=np.int32)
        self.offsets[1:] = np.cumsum(np.bincount(ends, minlength=natoms))

    def __getitem__(self, atom: int) -> np.ndarray:
        return self.flat[self.offsets[atom]:self.offsets[atom + 1]]

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def __iter__(self):
        return (self[a] for a in range(len(self)))

    @property
    def nbytes(self) -> int:
        return self.offsets.nbytes + self.flat.nbytes


@dataclass
class LewisStructure:
    """A molecule plus everything needed to draw it as a Lewis structure."""

    molecule: Molecule
    bonds: np.ndarray                 # (M, 2) 0-based atom indices
    orders: np.ndarray                # (M,) 1, 2 or 3
    charges: np.ndarray               # (N,) formal charges
    lone_pairs: np.ndarray            # (N,) non-bonding electron pairs
    radicals: np.ndarray              # (N,) unpaired electrons
    hydrogens: np.ndarray             # (N,) attached hydrogens
    neighbors: Adjacency | None = None
    hydrogens_missing: bool = False    # the file left them out (X-ray, say)
    net_charge: int = 0                # charge on the species, not on an atom
    net_radicals: int = 0              # unpaired electrons ditto

    @property
    def natoms(self) -> int:
        return self.molecule.natoms

    @property
    def total_charge(self) -> int:
        """Everything the drawing accounts for, on atoms and on the species."""
        return int(self.charges.sum()) + int(self.net_charge)

    @property
    def bracketed(self) -> bool:
        """Is there a charge or a radical that belongs to no single atom?"""
        return bool(self.net_charge or self.net_radicals)

    @property
    def matches_file(self) -> bool:
        """Does the perceived charge agree with the one in the file?

        A format that never states a charge cannot be disagreed with.
        """
        if not self.molecule.charge_known:
            return True
        return self.total_charge == int(self.molecule.charge)

    def heavy_neighbors(self, atom: int) -> list[int]:
        symbols = self.molecule.symbols
        return [j for j in self.neighbors[atom] if symbols[j] != "H"]

    def order_of(self, i: int, j: int) -> int:
        for b, (a, c) in enumerate(self.bonds):
            if (a, c) in ((i, j), (j, i)):
                return int(self.orders[b])
        return 0


def perceive(molecule: Molecule) -> LewisStructure:
    """Work out bond orders, formal charges and lone pairs for *molecule*."""
    bonds = molecule.bonds
    if bonds is None:
        bonds = molecule.perceive_bonds()
    bonds = np.asarray(bonds, dtype=np.int32).reshape(-1, 2)
    n = molecule.natoms
    z = molecule.atomic_numbers
    symbols = molecule.symbols

    neighbors = Adjacency(n, bonds)

    orders = _assign_orders(molecule, bonds, z, symbols, neighbors)
    bonded = _bonded_totals(bonds, orders, n)

    # int8 throughout: a formal charge is a single digit, a lone-pair count
    # smaller still, and at 100k atoms the default int is four megabytes
    charges, lone_pairs, radicals = _electron_bookkeeping(z, bonded)

    absent = _hydrogens_are_missing(z, bonded)
    if absent:
        # An X-ray structure has no hydrogens, so every carbon looks short of
        # a bond. Counting those as formal charges buries a protein under tens
        # of thousands of carbanions; the honest answer is that the file does
        # not say, so nothing is claimed.
        charges[:] = 0
        lone_pairs[:] = 0
        radicals[:] = 0
    else:
        _apply_radicals(molecule, charges, lone_pairs, radicals)

    net_charge, net_radicals = _reconcile(molecule, z, charges, lone_pairs,
                                          radicals, neighbors)

    hydrogens = _hydrogen_counts(z, bonds, n)
    return LewisStructure(molecule=molecule, bonds=bonds, orders=orders,
                          charges=charges, lone_pairs=lone_pairs,
                          radicals=radicals, hydrogens=hydrogens,
                          neighbors=neighbors, hydrogens_missing=absent,
                          net_charge=net_charge, net_radicals=net_radicals)


def _bonded_totals(bonds, orders, natoms: int) -> np.ndarray:
    """Sum of bond orders at each atom."""
    counts = orders.astype(np.int32)
    return (np.bincount(bonds[:, 0], counts, minlength=natoms)
            + np.bincount(bonds[:, 1], counts, minlength=natoms)).astype(np.int16)


def _hydrogen_counts(z, bonds, natoms: int) -> np.ndarray:
    """Hydrogens attached to each atom."""
    counts = np.zeros(natoms, dtype=np.int8)
    for near, far in ((0, 1), (1, 0)):
        hydrogen = z[bonds[:, far]] == 1
        counts += np.bincount(bonds[hydrogen, near],
                              minlength=natoms).astype(np.int8)
    return counts


# Per-element constants as arrays indexed by atomic number, so the electron
# count below is arithmetic over every atom at once rather than a dict lookup
# for each of a protein's hundred thousand.
_VALENCE_BY_Z = np.array([VALENCE_ELECTRONS.get(zi, 0) for zi in range(119)],
                         dtype=np.int16)
_SHELL_BY_Z = np.array([_SHELL.get(zi, 8) for zi in range(119)], dtype=np.int16)
_MAIN_GROUP_BY_Z = np.array([is_main_group(zi) for zi in range(119)], dtype=bool)


def _electron_bookkeeping(z, bonded):
    """Formal charges, lone pairs and unpaired electrons, from the bond orders.

    A free atom keeps all its electrons and is neutral; a bonded one fills its
    shell around the bonds it has, and whatever the count comes to short is
    its formal charge. Metals are left out of it entirely.
    """
    idx = np.clip(z, 0, 118)
    electrons = _VALENCE_BY_Z[idx]
    shell = _SHELL_BY_Z[idx]
    main = _MAIN_GROUP_BY_Z[idx]
    free = main & (bonded == 0)
    held = main & (bonded > 0)

    lone_pairs = np.zeros(len(z), dtype=np.int8)
    radicals = np.zeros(len(z), dtype=np.int8)
    charges = np.zeros(len(z), dtype=np.int8)
    lone_pairs[free] = electrons[free] // 2
    radicals[free] = electrons[free] % 2
    lone_pairs[held] = np.maximum(0, (shell[held] - 2 * bonded[held]) // 2)
    charges[held] = electrons[held] - 2 * lone_pairs[held] - bonded[held]
    return charges, lone_pairs, radicals


def _reconcile(molecule, z, charges, lone_pairs, radicals, neighbors):
    """Account for a charge the file declares but the drawing has not placed.

    An ion has the same connectivity as the neutral molecule, so geometry
    alone cannot find it — but the file says what it is, and that number was
    only ever being compared against, never used. Where the missing charge
    can go on one atom it goes there; where it cannot, as on a delocalized
    radical cation, it belongs to the species and is drawn the way a chemist
    draws it, in brackets with the charge outside.

    Returns the (charge, unpaired electrons) left over for those brackets.
    Metal complexes are left alone: their mismatch comes from drawing dative
    bonds as plain lines, which is a different argument.
    """
    if not molecule.charge_known:
        return 0, 0
    if any(not is_main_group(int(zi)) for zi in z):
        return 0, 0
    missing = int(molecule.charge) - int(charges.sum())
    if missing == 0:
        return 0, 0

    loose = [a for a in range(len(z)) if not len(neighbors[a])]
    if len(loose) == 1 and abs(missing) <= 8:     # a bare ion, or one in a
                                                  # fragment; anything wilder
                                                  # is not a formal charge
        a = loose[0]
        charges[a] += missing
        electrons = VALENCE_ELECTRONS[int(z[a])] - int(charges[a])
        if electrons >= 0:
            lone_pairs[a] = electrons // 2
            radicals[a] = electrons % 2
            return 0, 0
        charges[a] -= missing                 # more than the atom can give up

    unpaired = max(0, int(molecule.multiplicity) - 1 - int(radicals.sum()))
    return missing, unpaired


def _hydrogens_are_missing(z, bonded) -> bool:
    """Did the file simply leave the hydrogens out?

    Crystal structures rarely resolve them, and a structure solved with some
    of them still leaves most of its carbons looking under-coordinated.
    Carbon is the tell: it is four-valent with no room to argue, so a file
    where most carbons are short of four is one whose hydrogens are absent
    rather than one made of carbanions. A molecule that genuinely has no
    hydrogens — carbon dioxide, nitrate, a metal halide — has its carbons
    satisfied by multiple bonds and keeps its charges, which is where the
    interesting Lewis structures are.
    """
    carbons = [a for a in range(len(z)) if z[a] == 6]
    if len(carbons) < 4:
        return False
    short = sum(1 for a in carbons if bonded[a] < 4)
    return short > len(carbons) // 2


def _bond_lengths(coords, bonds) -> np.ndarray:
    """Every bond length at once — one call, not one per bond."""
    delta = coords[bonds[:, 0]] - coords[bonds[:, 1]]
    return np.sqrt(np.einsum("ij,ij->i", delta, delta))


def _wanted_orders(z, bonds, lengths):
    """The order each bond's length is closest to, and how contracted it is.

    Bonds are grouped by which pair of elements they join, so the reference
    lengths are looked up once per pair rather than once per bond: a protein
    has tens of thousands of bonds and perhaps twenty distinct pairs.
    """
    wanted = np.ones(len(bonds), dtype=np.int8)
    contraction = np.ones(len(bonds))
    za, zb = z[bonds[:, 0]].astype(np.int32), z[bonds[:, 1]].astype(np.int32)
    key = np.minimum(za, zb) * 256 + np.maximum(za, zb)
    for k in np.unique(key):
        refs = reference_lengths(elements.SYMBOLS[k // 256],
                                 elements.SYMBOLS[k % 256])
        if refs is None:
            continue
        here = key == k
        d = lengths[here]
        available = [(order, ref) for order, ref in enumerate(refs, start=1)
                     if ref is not None]
        errors = np.stack([np.abs(d - ref) for _, ref in available])
        # argmin takes the first on a tie, which is the lower order, exactly
        # as the tuple comparison it replaces did
        wanted[here] = np.array([o for o, _ in available])[errors.argmin(axis=0)]
        contraction[here] = d / refs[0]
    return wanted, contraction


# Indexed by atomic number: how many bonds an atom may hold, 0 for anything
# outside the main group, so capacity is a lookup rather than a loop.
_CAPACITY_BY_Z = np.array(
    [max_bonds(zi) if is_main_group(zi) else 0 for zi in range(119)],
    dtype=np.int16)


def _capacity(z, degree) -> np.ndarray:
    """How many further bond orders each atom can take (0 for the metals)."""
    top = _CAPACITY_BY_Z[np.clip(z, 0, len(_CAPACITY_BY_Z) - 1)]
    return np.maximum(0, top - degree).astype(np.int16)


def _assign_orders(molecule, bonds, z, symbols, neighbors) -> np.ndarray:
    """Bond orders from bond lengths, capped by the valence left on each atom."""
    orders = np.ones(len(bonds), dtype=np.int8)
    if not len(bonds):
        return orders

    lengths = _bond_lengths(molecule.coords, bonds)
    wanted, contraction = _wanted_orders(z, bonds, lengths)
    capacity = _capacity(z, np.diff(neighbors.offsets).astype(np.int16))

    # A bond that plainly wants a triple claims its valence before the doubles
    # are handed out, but only where both atoms can afford every triple they
    # want. Carbon dioxide's carbon wants two and can afford one, so it falls
    # through to the double pass and comes out O=C=O; a nitrile carbon wants
    # one and can afford it, so it is not left as C=N with a nitrogen anion.
    triple = wanted >= 3
    demand = (np.bincount(bonds[triple, 0], minlength=molecule.natoms)
              + np.bincount(bonds[triple, 1], minlength=molecule.natoms))
    affordable = 2 * demand <= capacity
    orders[triple & affordable[bonds[:, 0]] & affordable[bonds[:, 1]]] = 3

    for target in (2, 3):
        candidates = [b for b in range(len(bonds))
                      if wanted[b] >= target and orders[b] == target - 1]
        spare = _spare(capacity, bonds, orders)
        for b in sorted(candidates, key=_chooser(bonds, contraction, candidates,
                                                 molecule.natoms)):
            i, j = bonds[b]
            if spare[i] > 0 and spare[j] > 0:
                orders[b] = target
                spare[i] -= 1
                spare[j] -= 1
        while _augment(bonds, orders, spare, candidates, target):
            spare = _spare(capacity, bonds, orders)

    return orders

def _spare(capacity, bonds, orders) -> np.ndarray:
    """Valence each atom has left once the multiple bonds so far are counted."""
    extra = (orders - 1).astype(np.int32)
    n = len(capacity)
    return (capacity
            - np.bincount(bonds[:, 0], extra, minlength=n)
            - np.bincount(bonds[:, 1], extra, minlength=n)).astype(np.int16)


def _augment(bonds, orders, spare, candidates, target) -> bool:
    """Rescue a stranded atom by flipping an alternating chain of candidates.

    Greed hands the multiple bonds out one at a time and never takes one
    back, so it can paint itself into a corner. Ordering the choices (see
    :func:`_chooser`) is enough for one ring, but not for fused ones: from
    tetracene on, the acenes come out two double bonds short, with a pair of
    carbanions in the middle. Freeing them needs a chain that alternates
    unassigned and assigned candidates between two atoms with valence to
    spare — flipping it moves every partner along by one and leaves the
    structure one multiple bond richer.

    Returns True when it flipped something, so the caller can keep going.
    An odd-membered ring can hide a chain this depth-first walk misses, which
    is one more reason the perceived charge is checked against the file's.
    """
    free: dict[int, list[int]] = {}
    taken: dict[int, list[int]] = {}
    for b in candidates:
        side = taken if orders[b] == target else free
        i, j = bonds[b]
        side.setdefault(int(i), []).append(b)
        side.setdefault(int(j), []).append(b)

    def other(b, atom):
        i, j = bonds[b]
        return int(j) if int(i) == atom else int(i)

    def chain(atom, seen):
        """Alternating walk from *atom*, which is waiting for a partner."""
        for b in free.get(atom, ()):
            v = other(b, atom)
            if v in seen:
                continue
            seen.add(v)
            if spare[v] > 0:                          # v can simply take it
                return [b]
            for c in taken.get(v, ()):                # free v by moving its
                w = other(c, v)                       # partner w along
                if w in seen:
                    continue
                seen.add(w)
                rest = chain(w, seen)
                if rest is not None:
                    return [b, c] + rest
        return None

    for start in sorted(free):
        if spare[start] <= 0:
            continue
        flips = chain(start, {start})
        if flips is None:
            continue
        for step, b in enumerate(flips):              # assign, unassign, ...
            orders[b] = target if step % 2 == 0 else target - 1
        return True
    return False


def _chooser(bonds, contraction, candidates, natoms):
    """Sort key deciding which bond gets a multiple bond first.

    Most contracted first, so a real double bond takes the valence before a
    merely shortish neighbour can claim it. Bonds of much the same length are
    a tie — an aromatic ring is six of them — and there the most constrained
    bond wins: the one whose atoms have the fewest other candidates. Sorting
    once and walking the order beats re-picking the minimum after every
    assignment, which costs a protein tens of millions of comparisons;
    whatever the fixed order still strands, :func:`_augment` moves around.
    """
    options = np.zeros(natoms, dtype=int)
    for b in candidates:
        options[bonds[b][0]] += 1
        options[bonds[b][1]] += 1

    def key(b):
        i, j = bonds[b]
        return (round(contraction[b] / _TIE_WIDTH), options[i] + options[j],
                contraction[b])

    return key


def _apply_radicals(molecule, charges, lone_pairs, radicals) -> None:
    """Turn spare lone pairs into unpaired electrons for open-shell species.

    A doublet radical looks like an anion to the closed-shell rules above
    (methyl radical comes out as CH3-). The correction is only applied when it
    reconciles the perceived charge with the charge the file declares, so a
    guess never overrides what the calculation actually says.
    """
    unpaired = max(0, int(molecule.multiplicity) - 1)
    if not unpaired:
        return
    if int(charges.sum()) - int(molecule.charge) != -unpaired:
        return
    order = sorted((a for a in range(len(charges)) if lone_pairs[a] > 0),
                   key=lambda a: (charges[a], -lone_pairs[a]))
    for a in order[:unpaired]:
        lone_pairs[a] -= 1
        radicals[a] += 1
        charges[a] += 1
