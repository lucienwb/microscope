"""Gaussian basis sets: what the orbitals in a wavefunction file are built from.

A shell is one contraction of Gaussians on one centre. Everything a program's
conventions decide - the order of the components, how each function is
normalized, pure (spherical) or Cartesian, the odd flipped sign - is folded
into the shell's ``transform`` when the file is read, so the rest of the code
sees one thing: each basis function as a combination of the bare Cartesian
functions ``x^a y^b z^c sum_i w_i exp(-alpha_i r^2)`` in a fixed order.

The overlap matrix is here to check that a file was understood: a file's
orbitals are orthonormal in its own basis, so C^T S C = I holds only if every
convention was read the way the program meant it.

Lengths are in bohr throughout, as the files write them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from math import comb, pi, sqrt

import numpy as np

BOHR_TO_ANGSTROM = 0.52917721092
ANGULAR_LETTERS = "spdfghi"


def cartesian_components(l: int) -> np.ndarray:
    """(ncart, 3) exponents of x, y, z in the order used here: xx, xy, xz, yy, yz, zz."""
    return np.array([(a, b, l - a - b) for a in range(l, -1, -1)
                     for b in range(l - a, -1, -1)], dtype=np.int64)


def ncart(l: int) -> int:
    return (l + 1) * (l + 2) // 2


def component_exponents(label: str) -> tuple[int, int, int]:
    """'xyy' -> (1, 2, 0); 's' or '' -> (0, 0, 0)."""
    label = label.lower()
    if label == "s":
        return (0, 0, 0)
    return (label.count("x"), label.count("y"), label.count("z"))


def _double_factorial(n: int) -> int:
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def primitive_norm(exponents: np.ndarray, l: int) -> np.ndarray:
    """What makes x^l exp(-alpha r^2) a unit-norm function."""
    a = np.asarray(exponents, dtype=float)
    return (2 * a / pi) ** 0.75 * (4 * a) ** (l / 2) / sqrt(_double_factorial(2 * l - 1))


def solid_harmonic(l: int, m: int) -> dict[tuple[int, int, int], float]:
    """A real regular solid harmonic as a polynomial in x, y, z (unnormalized).

    Helgaker, Jorgensen & Olsen, *Molecular Electronic-Structure Theory*, eq.
    6.4.47. The signs are the common ones: d(+1) = xz, f(+3) = x^3 - 3xy^2,
    f(-3) = 3x^2y - y^3.
    """
    am = abs(m)
    parity = 0 if m >= 0 else 1          # 2 * v_m
    terms: dict[tuple[int, int, int], float] = {}
    for t in range((l - am) // 2 + 1):
        for u in range(t + 1):
            for k in range(parity, am + 1, 2):      # k = 2v
                sign = -1 if (t + (k - parity) // 2) % 2 else 1
                c = (sign * 0.25 ** t * comb(l, t) * comb(l - t, am + t)
                     * comb(t, u) * comb(am, k))
                key = (2 * t + am - 2 * u - k, 2 * u + k, l - 2 * t - am)
                terms[key] = terms.get(key, 0.0) + c
    return terms


def pure_m_order(l: int) -> list[int]:
    """The order nearly every program lists pure functions in: 0, +1, -1, +2, -2, ..."""
    return [0] + [s * m for m in range(1, l + 1) for s in (1, -1)]


def pure_rows(l: int, ms: list[int]) -> np.ndarray:
    """(len(ms), ncart) solid harmonics in the Cartesian components, unnormalized."""
    index = {tuple(c): i for i, c in enumerate(cartesian_components(l))}
    rows = np.zeros((len(ms), ncart(l)))
    for r, m in enumerate(ms):
        for key, c in solid_harmonic(l, m).items():
            rows[r, index[key]] = c
    return rows


def cartesian_rows(l: int, labels: list[str]) -> np.ndarray:
    """(n, ncart) selecting the Cartesian components in a file's order."""
    index = {tuple(c): i for i, c in enumerate(cartesian_components(l))}
    rows = np.zeros((len(labels), ncart(l)))
    for r, label in enumerate(labels):
        rows[r, index[component_exponents(label)]] = 1.0
    return rows


@dataclass
class Shell:
    """One contracted shell, already in the form the evaluator wants.

    ``transform`` has one row per basis function as the file counts them, in
    the file's order, giving that function in the Cartesian components of
    :func:`cartesian_components`; ``weights`` multiply the bare Gaussians.
    """

    center: np.ndarray            # (3,) bohr
    l: int
    exponents: np.ndarray         # (nprim,) bohr^-2
    weights: np.ndarray           # (nprim,)
    transform: np.ndarray         # (nfunc, ncart)
    atom: int = -1                # which atom of the file it sits on

    @property
    def nfunc(self) -> int:
        return len(self.transform)

    @property
    def pure(self) -> bool:
        return self.l >= 2 and self.nfunc == 2 * self.l + 1


@dataclass
class BasisSet:
    shells: list[Shell] = field(default_factory=list)

    @property
    def nbf(self) -> int:
        return sum(s.nfunc for s in self.shells)

    def offsets(self) -> np.ndarray:
        """Index of each shell's first basis function."""
        sizes = [s.nfunc for s in self.shells]
        return np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(np.int64)

    @property
    def max_l(self) -> int:
        return max((s.l for s in self.shells), default=0)

    def overlap(self) -> np.ndarray:
        """The (nbf, nbf) overlap matrix of the basis functions."""
        return overlap_matrix(self.shells)


def make_shell(center, l: int, exponents, coefficients, rows: np.ndarray, *,
               normalized_primitives: bool = True, normalize: str = "function",
               atom: int = -1) -> Shell:
    """A shell from what a file says about it.

    ``rows`` are the file's functions in Cartesian components, before any
    normalization (from :func:`pure_rows` or :func:`cartesian_rows`).
    ``normalized_primitives``: the contraction coefficients multiply
    unit-norm primitives, as almost every program writes them.
    ``normalize``: "function" makes every basis function unit-norm, and
    "turbomole" makes each one sqrt((2l-1)!!) times that. "shell" and "raw"
    scale a Cartesian shell as a whole, so its components keep the ratios of
    their polynomials: "shell" makes x^l unit-norm and leaves xy at
    1/sqrt(3) of it (Psi4 1.3), "raw" leaves out the double factorial, so xy
    is the unit-norm one and xx is sqrt(3) (CFOUR).
    """
    exponents = np.asarray(exponents, dtype=float)
    weights = np.asarray(coefficients, dtype=float).copy()
    if normalized_primitives:
        weights = weights * primitive_norm(exponents, l)
    shell = Shell(center=np.asarray(center, dtype=float), l=l, exponents=exponents,
                  weights=weights, transform=np.array(rows, dtype=float), atom=atom)
    block = _one_centre(shell)                                # (ncart, ncart)
    if normalize in ("function", "turbomole"):
        norms = np.einsum("ij,jk,ik->i", shell.transform, block, shell.transform)
        # a contraction of nothing but zeros stays zero rather than turning to NaN
        shell.transform /= np.sqrt(np.where(norms > 0, norms, 1.0))[:, None]
        if normalize == "turbomole":
            shell.transform *= sqrt(_double_factorial(2 * l - 1))
    elif normalize in ("shell", "raw"):
        shell.transform /= sqrt(block[0, 0]) if block[0, 0] > 0 else 1.0   # x^l is first
        if normalize == "raw":
            shell.transform *= sqrt(_double_factorial(2 * l - 1))
    else:
        raise ValueError(f"unknown normalization {normalize!r}")
    return shell


@cache
def _angular_overlap(l: int) -> np.ndarray:
    """(ncart, ncart) of prod_k (n_k - 1)!! over the summed exponents n_k, zero
    if any is odd: the angular half of two components on one centre."""
    comps = cartesian_components(l)
    total = comps[:, None, :] + comps[None, :, :]
    even = (total % 2 == 0).all(axis=2)
    values = np.vectorize(_double_factorial)(total - 1).prod(axis=2).astype(float)
    return np.where(even, values, 0.0)


def _one_centre(shell: Shell) -> np.ndarray:
    """Overlap of a shell's bare Cartesian components with themselves:
    int x^A y^B z^C exp(-p r^2) = prod (n-1)!! / (2p)^(n/2) * (pi/p)^(3/2), so a
    radial sum over primitive pairs times a fixed matrix for each l."""
    p = shell.exponents[:, None] + shell.exponents[None, :]
    radial = (shell.weights[:, None] * shell.weights[None, :]
              * (pi / p) ** 1.5 / (2 * p) ** shell.l).sum()
    return radial * _angular_overlap(shell.l)


# ------------------------------------------------------------------ overlap

def _overlap_1d(la: int, lb: int, xpa, xpb, p):
    """Obara-Saika table s[i][j] of int x_A^i x_B^j exp(-p x_P^2) dx / sqrt(pi/p)."""
    half = 0.5 / p
    s = [[None] * (lb + 1) for _ in range(la + 1)]
    s[0][0] = np.ones_like(p)
    for i in range(la):
        s[i + 1][0] = xpa * s[i][0] + (half * i * s[i - 1][0] if i else 0.0)
    for j in range(lb):
        for i in range(la + 1):
            term = xpb * s[i][j]
            if i:
                term = term + half * i * s[i - 1][j]
            if j:
                term = term + half * j * s[i][j - 1]
            s[i][j + 1] = term
    return np.array(s)                  # (la+1, lb+1, *p.shape)


def _primitives(shells: list[Shell]):
    exps = np.concatenate([s.exponents for s in shells])
    weights = np.concatenate([s.weights for s in shells])
    centers = np.concatenate([np.repeat(s.center[None], len(s.exponents), axis=0)
                              for s in shells])
    sizes = [len(s.exponents) for s in shells]
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(np.int64)
    return exps, weights, centers, starts


def _cartesian_block(shells_a: list[Shell], shells_b: list[Shell]) -> np.ndarray:
    """(na, nb, ncart_a, ncart_b) overlaps of the bare Cartesian components.

    All shells in a list share one l; the work is vectorized over every pair of
    primitives and summed back to shells with reduceat.
    """
    la, lb = shells_a[0].l, shells_b[0].l
    ea, wa, ca, sa = _primitives(shells_a)
    eb, wb, cb, sb = _primitives(shells_b)
    p = ea[:, None] + eb[None, :]
    mu = ea[:, None] * eb[None, :] / p
    diff = ca[:, None, :] - cb[None, :, :]                   # A - B
    prefactor = (pi / p) ** 1.5 * np.exp(-mu * np.einsum("abk,abk->ab", diff, diff))
    prefactor *= wa[:, None] * wb[None, :]
    centre_p = (ea[:, None, None] * ca[:, None, :] + eb[None, :, None] * cb[None, :, :])
    centre_p /= p[..., None]
    comps_a, comps_b = cartesian_components(la), cartesian_components(lb)
    out = prefactor[:, :, None, None]
    for k in range(3):
        table = _overlap_1d(la, lb, centre_p[..., k] - ca[:, None, k],
                            centre_p[..., k] - cb[None, :, k], p)
        picked = table[comps_a[:, k][:, None], comps_b[:, k][None, :]]   # (ca, cb, Pa, Pb)
        out = out * np.moveaxis(picked, (0, 1), (2, 3))
    out = np.add.reduceat(out, sa, axis=0)
    return np.add.reduceat(out, sb, axis=1)


_CHUNK = 2_000_000       # primitive-pair components computed at once


def _overlap_blocks(shells: list[Shell]):
    """The overlap matrix in the files' own basis functions (their order, their
    normalization), one angular-momentum pair and one bounded chunk at a time:
    yields (rows, columns, block, same_l) with rows <= columns in l."""
    sizes = [s.nfunc for s in shells]
    offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
    by_l: dict[int, list[int]] = {}
    for i, s in enumerate(shells):
        by_l.setdefault(s.l, []).append(i)

    def padded(indices):
        """Transforms padded to ncart rows, with the basis-function index of each row."""
        l = shells[indices[0]].l
        T = np.zeros((len(indices), ncart(l), ncart(l)))
        rows = np.full((len(indices), ncart(l)), -1, dtype=np.int64)
        for n, i in enumerate(indices):
            f = shells[i].nfunc
            T[n, :f] = shells[i].transform
            rows[n, :f] = offsets[i] + np.arange(f)
        return T, rows

    for la in sorted(by_l):
        for lb in sorted(by_l):
            if lb < la:
                continue
            ib = by_l[lb]
            Tb, rows_b = padded(ib)
            keep_b = rows_b.ravel() >= 0
            nprim_b = sum(len(shells[i].exponents) for i in ib)
            per_prim = max(1, nprim_b * ncart(la) * ncart(lb))
            chunk: list[int] = []
            count = 0
            ia = by_l[la]
            for n, i in enumerate(ia):
                chunk.append(i)
                count += len(shells[i].exponents)
                if count * per_prim < _CHUNK and n + 1 < len(ia):
                    continue
                Ta, rows_a = padded(chunk)
                block = _cartesian_block([shells[i] for i in chunk], [shells[i] for i in ib])
                block = np.einsum("aim,abmn,bjn->aibj", Ta, block, Tb, optimize=True)
                keep_a = rows_a.ravel() >= 0
                flat = block.reshape(rows_a.size, rows_b.size)[keep_a][:, keep_b]
                yield rows_a.ravel()[keep_a], rows_b.ravel()[keep_b], flat, la == lb
                chunk, count = [], 0


def overlap_matrix(shells: list[Shell]) -> np.ndarray:
    """The (nbf, nbf) overlap matrix in the files' own basis functions."""
    nbf = sum(s.nfunc for s in shells)
    S = np.zeros((nbf, nbf))
    for rows, columns, block, _ in _overlap_blocks(shells):
        S[np.ix_(rows, columns)] = block
        S[np.ix_(columns, rows)] = block.T
    return S


def orbital_overlaps(shells: list[Shell], coefficients: np.ndarray) -> np.ndarray:
    """C S C^T for orbitals (k, nbf) in these shells' functions, without the
    overlap matrix ever being held: memory is the orbitals and one block, where
    S itself is 650 MB at 6000 basis functions."""
    C = np.asarray(coefficients, dtype=float)
    gram = np.zeros((len(C), len(C)))
    for rows, columns, block, same_l in _overlap_blocks(shells):
        part = C[:, rows] @ block @ C[:, columns].T
        # an l with itself: the chunks' rows cover the group once, so every
        # pair is here once; two different l: the mirror pairs are the transpose
        gram += part if same_l else part + part.T
    return gram


def to_cartesian(basis: BasisSet, coefficients: np.ndarray) -> np.ndarray:
    """Orbitals (nmo, nbf) re-expanded in the bare Cartesian components."""
    C = np.asarray(coefficients, dtype=float)
    out = np.empty((len(C), sum(ncart(s.l) for s in basis.shells)))
    column = 0
    for shell, start in zip(basis.shells, basis.offsets()):
        n = ncart(shell.l)
        out[:, column:column + n] = C[:, start:start + shell.nfunc] @ shell.transform
        column += n
    return out


def orthonormality_error(shells: list[Shell], coefficients: np.ndarray) -> float:
    """Largest deviation of C^T S C from the identity; coefficients are (nmo, nbf)."""
    if not len(coefficients):
        return 0.0
    gram = orbital_overlaps(shells, coefficients)
    return float(np.abs(gram - np.eye(len(gram))).max())


def angular_momentum(letter: str) -> int:
    letter = letter.lower()
    if letter not in ANGULAR_LETTERS:
        raise ValueError(f"unknown shell type {letter!r}")
    return ANGULAR_LETTERS.index(letter)
