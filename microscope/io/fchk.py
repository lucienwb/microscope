"""Gaussian formatted checkpoint (.fchk) reader: geometry, and the wavefunction
(basis set and molecular orbitals) when the file carries one."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..core import elements
from ..core.basis import BasisSet, cartesian_rows, make_shell, pure_m_order, pure_rows
from ..core.molecule import Molecule
from ..core.orbitals import ORTHONORMAL, Orbitals, OrbitalSet
from ..core.results import ParseResult
from .errors import FileFormatError

BOHR_TO_ANGSTROM = 0.52917721092

# Gaussian's order for Cartesian shells up to f; from g on it is systematic
_CARTESIAN = {
    0: ["s"],
    1: ["x", "y", "z"],
    2: ["xx", "yy", "zz", "xy", "xz", "yz"],
    3: ["xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"],
}


def cartesian_order(l: int) -> list[str]:
    """Gaussian's Cartesian components: zzzz, yzzz, yyzz, ... xxxx from g on."""
    if l in _CARTESIAN:
        return _CARTESIAN[l]
    return ["x" * a + "y" * b + "z" * (l - a - b)
            for a in range(l + 1) for b in range(l - a + 1)]


_HEADER = re.compile(r"^\S[^\n]*", re.M)     # data lines start with a space


class _Fields:
    """Where each named field is, so a large array is parsed only if asked for,
    straight from its stretch of the file: splitting millions of numbers into
    Python strings first costs several times the file in memory."""

    def __init__(self, text: str):
        self.text = text
        self.scalars: dict[str, str] = {}
        self.arrays: dict[str, tuple[str, int, int, int]] = {}   # type, start, end, count
        headers = []
        for match in _HEADER.finditer(text):
            line = match.group(0).rstrip("\r")
            if len(line) < 44:
                continue
            name, kind, rest = line[:43].strip(), line[43], line[44:]
            headers.append((name, kind, rest, match.start(), match.end()))
        for n, (name, kind, rest, _, end) in enumerate(headers):
            if "N=" in rest:
                stop = headers[n + 1][3] if n + 1 < len(headers) else len(text)
                self.arrays[name] = (kind, end, stop, int(rest.split("N=")[1]))
            elif kind in "IR":
                self.scalars[name] = rest.strip()

    def array(self, name: str, dtype=float) -> np.ndarray | None:
        if name not in self.arrays:
            return None
        _, start, stop, count = self.arrays[name]
        values = np.fromstring(self.text[start:stop], dtype=dtype, sep=" ")[:count]
        if len(values) != count:
            raise FileFormatError(f"{name!r} is cut short ({len(values)} of {count})")
        return values

    def integer(self, name: str, default: int = 0) -> int:
        try:
            return int(self.scalars[name])
        except (KeyError, ValueError):
            return default

    @property
    def title(self) -> str:
        return self.text[:self.text.find("\n")].strip() if self.text else ""


def read(path) -> ParseResult:
    fields = _Fields(Path(path).read_text(errors="replace"))
    try:
        zs = fields.array("Atomic numbers", dtype=np.int64)
        coords = fields.array("Current cartesian coordinates")
    except FileFormatError as exc:
        raise FileFormatError(f"{path}: {exc}") from None
    if zs is None or coords is None:
        raise FileFormatError(f"{path}: no geometry found in fchk file")
    coords = coords.reshape(-1, 3) * BOHR_TO_ANGSTROM
    if len(zs) != len(coords):
        raise FileFormatError(f"{path}: atom count mismatch in fchk file")

    mol = Molecule(
        symbols=[elements.SYMBOLS[z] if 0 < z < len(elements.SYMBOLS) else "X" for z in zs],
        coords=coords,
        charge=fields.integer("Charge"),
        multiplicity=fields.integer("Multiplicity", 1),
        title=fields.title,
    )
    # a wavefunction that cannot be read costs the orbitals, not the structure
    warnings = []
    try:
        orbitals = _orbitals(fields)
    except (FileFormatError, ValueError, IndexError) as exc:
        orbitals = None
        warnings.append(f"the orbitals could not be read: {exc}")
    if orbitals is not None and orbitals.reading_error() >= ORTHONORMAL:
        # Gaussian's conventions always hold for Gaussian's files; other
        # programs write fchk too, and one that means something else is
        # flagged rather than drawn wrong
        orbitals.convention = "unrecognized"
    return ParseResult(frames=[mol], program="Gaussian fchk", source=str(path),
                       orbitals=orbitals, warnings=warnings)


def _basis(fields: _Fields) -> BasisSet | None:
    types = fields.array("Shell types", dtype=np.int64)
    if types is None:
        return None
    nprims = fields.array("Number of primitives per shell", dtype=np.int64)
    atoms = fields.array("Shell to atom map", dtype=np.int64)
    exponents = fields.array("Primitive exponents")
    coefficients = fields.array("Contraction coefficients")
    sp_coefficients = fields.array("P(S=P) Contraction coefficients")
    centers = fields.array("Coordinates of each shell")
    if nprims is None or atoms is None or exponents is None or coefficients is None \
            or centers is None:
        raise FileFormatError("the basis set is incomplete")
    centers = centers.reshape(-1, 3)

    shells = []
    start = 0
    for n, (kind, count) in enumerate(zip(types, nprims)):
        kind, count = int(kind), int(count)
        exps = exponents[start:start + count]
        coefs = coefficients[start:start + count]
        where = {"center": centers[n], "atom": int(atoms[n]) - 1}
        if kind == -1:                       # SP: an s and a p sharing exponents
            if sp_coefficients is None:
                raise FileFormatError("an SP shell without its P coefficients")
            shells.append(make_shell(l=0, exponents=exps, coefficients=coefs,
                                     rows=cartesian_rows(0, ["s"]), **where))
            shells.append(make_shell(l=1, exponents=exps,
                                     coefficients=sp_coefficients[start:start + count],
                                     rows=cartesian_rows(1, cartesian_order(1)), **where))
        else:
            l = abs(kind)
            rows = (pure_rows(l, pure_m_order(l)) if kind < -1
                    else cartesian_rows(l, cartesian_order(l)))
            shells.append(make_shell(l=l, exponents=exps, coefficients=coefs,
                                     rows=rows, **where))
        start += count
    return BasisSet(shells)


def _orbitals(fields: _Fields) -> Orbitals | None:
    if "Alpha MO coefficients" not in fields.arrays:
        return None
    basis = _basis(fields)
    if basis is None:
        return None
    nbf = basis.nbf
    nalpha = fields.integer("Number of alpha electrons")
    nbeta = fields.integer("Number of beta electrons")

    def spin_set(prefix: str, occupations_for) -> OrbitalSet:
        energies = fields.array(f"{prefix} Orbital Energies")
        coefficients = fields.array(f"{prefix} MO coefficients")
        assert coefficients is not None
        if len(coefficients) % nbf:
            raise FileFormatError(f"{prefix} MO coefficients do not fit "
                                  f"{nbf} basis functions")
        coefficients = coefficients.reshape(-1, nbf)
        nmo = len(coefficients)
        if energies is None or len(energies) != nmo:
            energies = np.full(nmo, np.nan)
        return OrbitalSet(energies=energies, occupations=occupations_for(nmo),
                          coefficients=coefficients.astype(np.float32))

    def filled(nmo: int, *counts: int) -> np.ndarray:
        occupations = np.zeros(nmo, dtype=np.float32)
        for n in counts:
            occupations[:min(n, nmo)] += 1
        return occupations

    if "Beta MO coefficients" in fields.arrays:
        return Orbitals(basis=basis,
                        alpha=spin_set("Alpha", lambda n: filled(n, nalpha)),
                        beta=spin_set("Beta", lambda n: filled(n, nbeta)))
    # restricted, closed or open shell: doubly occupied, then singly
    return Orbitals(basis=basis, alpha=spin_set("Alpha", lambda n: filled(n, nalpha, nbeta)))
