"""Molden format reader: [Atoms], [GEOMETRIES], [FREQ]/[FR-NORM-COORD], and the
wavefunction in [GTO] + [MO].

Molden is written by many programs and they do not all mean the same thing by
it: how contraction coefficients and Cartesian functions are normalized, and
the signs of some pure f and g functions, differ. Rather than trust a title
line, the reader checks each reading against the file's own orbitals - they
must come out orthonormal - and keeps the one that does.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core import elements
from ..core.basis import (
    BasisSet,
    angular_momentum,
    cartesian_rows,
    make_shell,
    pure_m_order,
    pure_rows,
)
from ..core.molecule import Molecule
from ..core.orbitals import ORTHONORMAL, Orbitals, OrbitalSet
from ..core.results import ParseResult, Vibration
from .errors import FileFormatError

BOHR_TO_ANGSTROM = 0.52917721092

_SECTION_RE = re.compile(r"[ \t]*\[([A-Za-z0-9_ -]+)\][ \t]*(.*)")


def looks_like(text: str) -> bool:
    return "[molden format]" in text.lower()


def _split_sections(text: str) -> dict[str, tuple[str, str]]:
    """Section name -> (the text after the name, the section's own text).

    Headers are found by looking for '[' rather than trying a pattern at every
    line: the [MO] block of a large file is millions of lines of numbers with
    no bracket in them, and is handed on as one string, never split.
    """
    headers = []
    at = text.find("[")
    while at >= 0:
        start = text.rfind("\n", 0, at) + 1
        end = text.find("\n", at)
        end = len(text) if end < 0 else end
        m = _SECTION_RE.match(text, start, end)
        if m:
            headers.append((m.group(1).strip().lower(), m.group(2).strip(), start, end))
        at = text.find("[", end)
    sections: dict[str, tuple[str, str]] = {}
    for k, (name, arg, _, end) in enumerate(headers):
        stop = headers[k + 1][2] if k + 1 < len(headers) else len(text)
        sections[name] = (arg, text[end + 1:stop])
    return sections


def _lines(sections, name: str) -> list[str]:
    return sections[name][1].splitlines() if name in sections else []


def read(path) -> ParseResult:
    sections = _split_sections(Path(path).read_text(errors="replace"))
    if not sections:
        raise FileFormatError(f"{path}: not a Molden file")

    frames: list[Molecule] = []

    # primary geometry from [Atoms]
    natoms = 0
    centers: dict[int, np.ndarray] = {}     # atom number -> position (bohr), ghosts too
    if "atoms" in sections:
        arg, body = sections["atoms"][0], _lines(sections, "atoms")
        to_ang = BOHR_TO_ANGSTROM if _in_bohr(arg) else 1.0
        syms, xyz = [], []
        for row in body:
            parts = row.split()
            if len(parts) < 6:
                continue
            try:
                number, z = int(parts[1]), int(parts[2])
                coords3 = [float(v) * to_ang for v in parts[3:6]]
            except ValueError:
                continue
            centers[number] = np.array(coords3) / BOHR_TO_ANGSTROM
            if z == 0:                 # a ghost: basis functions, but no atom
                continue
            # The name says which element it is. The number is the nuclear
            # charge the program used, which with a pseudopotential is the
            # core charge: ORCA writes 17 for a Rh with a 28-electron core.
            named = elements.symbol_to_z(parts[0])
            syms.append(elements.SYMBOLS[named] if named else elements.normalize_symbol(z))
            xyz.append(coords3)
        if syms:
            frames.append(Molecule(syms, np.array(xyz), charge_known=False))
            natoms = len(syms)

    # optimization history from [GEOMETRIES] XYZ
    if "geometries" in sections:
        arg, body = sections["geometries"][0], _lines(sections, "geometries")
        if arg.upper().startswith("XYZ") or not arg:
            geo_frames = _parse_xyz_frames(body)
            if geo_frames:
                frames = geo_frames  # full trajectory supersedes the single frame
                natoms = frames[-1].natoms

    # frequency job geometry from [FR-COORD] (bohr) if nothing else
    if not frames and "fr-coord" in sections:
        body = _lines(sections, "fr-coord")
        syms, xyz = [], []
        for row in body:
            parts = row.split()
            if len(parts) >= 4:
                try:
                    coords3 = [float(v) * BOHR_TO_ANGSTROM for v in parts[1:4]]
                except ValueError:
                    continue
                syms.append(parts[0])
                xyz.append(coords3)
        if syms:
            frames.append(Molecule(syms, np.array(xyz), charge_known=False))
            natoms = len(syms)

    if not frames:
        raise FileFormatError(f"{path}: no geometry found in Molden file")

    vibrations: list[Vibration] = []
    if "freq" in sections:
        freqs = [float(t) for t in sections["freq"][1].split()]
        intensities = []
        if "int" in sections:
            intensities = [float(t) for t in sections["int"][1].split()]
        disps: list[np.ndarray] = []
        if "fr-norm-coord" in sections and natoms:
            current: list[list[float]] = []
            for row in _lines(sections, "fr-norm-coord"):
                if row.strip().lower().startswith("vibration"):
                    if current:
                        disps.append(np.array(current) * BOHR_TO_ANGSTROM)
                    current = []
                else:
                    parts = row.split()
                    if len(parts) >= 3:
                        try:
                            current.append([float(v) for v in parts[:3]])
                        except ValueError:
                            pass
            if current:
                disps.append(np.array(current) * BOHR_TO_ANGSTROM)
        for k, f in enumerate(freqs):
            if abs(f) < 1e-3:
                continue
            d = disps[k] if k < len(disps) and len(disps[k]) == natoms else None
            vibrations.append(Vibration(
                frequency=f,
                ir_intensity=intensities[k] if k < len(intensities) else None,
                displacements=d))

    # a wavefunction that cannot be read costs the orbitals, not the structure
    warnings = []
    try:
        orbitals = _orbitals(sections, centers)
    except (ValueError, IndexError) as exc:
        orbitals = None
        warnings.append(f"the orbitals could not be read: {exc}")
    return ParseResult(frames=frames, program="Molden", source=str(path),
                       vibrations=vibrations, orbitals=orbitals, warnings=warnings)


def _in_bohr(arg: str) -> bool:
    arg = arg.lower()
    return "au" in arg or "bohr" in arg


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _parse_xyz_frames(body: list[str]) -> list[Molecule]:
    frames: list[Molecule] = []
    i = 0
    while i < len(body):
        stripped = body[i].strip()
        if not stripped:
            i += 1
            continue
        try:
            count = int(stripped.split()[0])
        except ValueError:
            break
        rows = body[i + 2:i + 2 + count]
        syms, xyz = [], []
        for row in rows:
            parts = row.split()
            if len(parts) >= 4:
                syms.append(parts[0])
                xyz.append([float(v) for v in parts[1:4]])
        if len(syms) == count:
            frames.append(Molecule(syms, np.array(xyz), charge_known=False))
        i += 2 + count
    return frames


# ------------------------------------------------------------- wavefunction

# Molden's order for Cartesian shells (the g order is its own); h and above
# have no Cartesian order in the format, so they are always pure
_CARTESIAN = {
    0: ["s"],
    1: ["x", "y", "z"],
    2: ["xx", "yy", "zz", "xy", "xz", "yz"],
    3: ["xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"],
    4: ["xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy",
        "xxyy", "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"],
}


@dataclass
class _RawShell:
    atom: int                     # atom number in [Atoms]
    l: int
    exponents: np.ndarray
    coefficients: np.ndarray


@dataclass(frozen=True)
class Convention:
    """One way a program may have meant its [GTO] section."""

    normalized_primitives: bool = True     # coefficients multiply unit-norm primitives
    cartesian: str = "function"            # how Cartesian d, f, g are normalized
    orca_signs: bool = False               # pure functions with |m| = 3, 4 flipped

    @property
    def name(self) -> str:
        parts = []
        if not self.normalized_primitives:
            parts.append("primitive norms in the coefficients")
        if self.cartesian != "function":
            parts.append(_CARTESIAN_NAMES[self.cartesian])
        if self.orca_signs:
            parts.append("ORCA's signs for pure f, g, h")
        return ", ".join(parts) or "standard"


_CARTESIAN_NAMES = {"shell": "Cartesians normalized as x^l",
                    "raw": "Cartesians normalized as xy",
                    "turbomole": "Cartesians scaled by (2l-1)!!"}

# The Molden format as written down, then what programs actually write: older
# Psi4 puts the primitive norms into the coefficients, ORCA does that and also
# flips the sign of pure functions with |m| = 3 or 4, and Psi4 1.3, CFOUR and
# Turbomole each scale Cartesian d, f and g their own way. Every other
# combination follows, in case a program mixes them. The order only breaks
# ties, which happen when a file cannot tell two readings apart (one atom
# cannot show a sign convention; it needs a neighbour to overlap with).
_ORCA = Convention(normalized_primitives=False, orca_signs=True)
_LIKELY = (
    Convention(),
    Convention(normalized_primitives=False),
    _ORCA,
    Convention(cartesian="shell"),
    Convention(cartesian="raw"),
    Convention(cartesian="turbomole"),
)
CONVENTIONS = _LIKELY + tuple(
    c for c in (Convention(primitives, cartesian, signs)
                for primitives in (True, False) for signs in (False, True)
                for cartesian in ("function", "shell", "raw", "turbomole"))
    if c not in _LIKELY)
TIE = 4.0               # readings within this factor of the best are as good


def _read_gto(body: list[str]) -> list[_RawShell]:
    shells: list[_RawShell] = []
    atom = 0
    i = 0
    while i < len(body):
        parts = body[i].split()
        i += 1
        if not parts:
            continue
        if parts[0].isdigit():                 # "  3 0": the next atom's shells
            atom = int(parts[0])
            continue
        label, nprim = parts[0].lower(), int(parts[1])
        scale = _number(parts[2]) if len(parts) > 2 else 1.0
        rows = np.array([[_number(tok) for tok in body[i + k].split()]
                         for k in range(nprim)])
        i += nprim
        exponents = rows[:, 0] * scale ** 2
        if label == "sp":
            shells.append(_RawShell(atom, 0, exponents, rows[:, 1]))
            shells.append(_RawShell(atom, 1, exponents, rows[:, 2]))
        else:
            shells.append(_RawShell(atom, angular_momentum(label), exponents, rows[:, 1]))
    return shells


def _pure_shells(sections) -> set[int]:
    """Which angular momenta the file writes as pure functions."""
    pure = {5, 6}
    if "5d" in sections or "5d7f" in sections:
        pure |= {2, 3}
    if "5d10f" in sections:
        pure.add(2)
    if "7f" in sections:
        pure.add(3)
    if "9g" in sections:
        pure.add(4)
    return pure


def _build_basis(raw: list[_RawShell], centers: dict[int, np.ndarray], pure: set[int],
                 convention: Convention) -> BasisSet:
    shells = []
    for shell in raw:
        l = shell.l
        if l in pure:
            ms = pure_m_order(l)
            rows = pure_rows(l, ms)
            if convention.orca_signs:
                rows[[abs(m) in (3, 4) for m in ms]] *= -1
        else:
            rows = cartesian_rows(l, _CARTESIAN[l])
        if shell.atom not in centers:
            raise ValueError(f"[GTO] names atom {shell.atom}, which [Atoms] does not have")
        shells.append(make_shell(
            centers[shell.atom], l, shell.exponents, shell.coefficients, rows,
            normalized_primitives=convention.normalized_primitives,
            normalize=convention.cartesian if l >= 2 and l not in pure else "function",
            atom=shell.atom - 1))
    return BasisSet(shells)


def _key_lines(text: str):
    """(start, end, key, value) of each 'Key= value' line, found by their '='
    signs, which only the few key lines have."""
    at = text.find("=")
    while at >= 0:
        start = text.rfind("\n", 0, at) + 1
        end = text.find("\n", at)
        end = len(text) if end < 0 else end
        yield start, end, text[start:at].strip().lower(), text[at + 1:end].strip()
        at = text.find("=", end)


def _read_mo(text: str) -> list[dict]:
    """Each orbital's keys, and its coefficients as (indices, values) arrays.

    The numbers between one orbital's keys and the next's are parsed in one
    go: a large ORCA file is millions of 'index value' lines, and a Python
    loop over them was most of the time it took to open one.
    """
    orbitals: list[dict] = []
    current: dict | None = None
    position = 0

    def numbers(chunk: str) -> None:
        if current is None or not chunk.strip():
            return
        if "D" in chunk or "d" in chunk:
            chunk = chunk.replace("D", "E").replace("d", "e")
        pairs = np.array(chunk.split(), dtype=float)
        if len(pairs) % 2:
            raise ValueError("an orbital's coefficients are not 'index value' pairs")
        pairs = pairs.reshape(-1, 2)
        current["indices"].append(pairs[:, 0].astype(np.int64) - 1)
        current["values"].append(pairs[:, 1])

    for start, end, key, value in _key_lines(text):
        numbers(text[position:start])
        position = end
        if current is None or current["indices"]:
            current = {"sym": "", "energy": np.nan, "spin": "alpha", "occupation": 0.0,
                       "indices": [], "values": []}
            orbitals.append(current)
        if key == "sym":
            current["sym"] = value
        elif key == "ene":
            current["energy"] = _number(value)
        elif key == "spin":
            current["spin"] = "beta" if value.lower().startswith("b") else "alpha"
        elif key == "occup":
            current["occupation"] = _number(value)
    numbers(text[position:])
    for orbital in orbitals:
        orbital["indices"] = (np.concatenate(orbital["indices"]) if orbital["indices"]
                              else np.zeros(0, dtype=np.int64))
        orbital["values"] = (np.concatenate(orbital["values"]) if orbital["values"]
                             else np.zeros(0))
    return [o for o in orbitals if len(o["indices"])]


def _counted_pure(raw: list[_RawShell], declared: set[int], nbf: int) -> set[int]:
    """Which shells are pure, settled by counting when the flags disagree with
    the orbitals: a file that writes five d functions but no [5D] means five."""
    def count(pure):
        return sum(2 * s.l + 1 if s.l in pure else (s.l + 1) * (s.l + 2) // 2 for s in raw)
    if count(declared) == nbf:
        return declared
    ambiguous = sorted({s.l for s in raw} & {2, 3, 4})
    for size in range(1, len(ambiguous) + 1):
        for flipped in _subsets(ambiguous, size):
            pure = declared ^ set(flipped)
            if count(pure) == nbf:
                return pure
    if nbf < count(declared):      # trailing zero coefficients left out
        return declared
    raise ValueError(f"the orbitals have {nbf} coefficients, which no reading of the "
                     f"basis set gives (Cartesian: {count(set())}, pure: {count({2, 3, 4, 5, 6})})")


def _subsets(items: list[int], size: int):
    if size == 0:
        yield ()
        return
    for i, item in enumerate(items):
        for rest in _subsets(items[i + 1:], size - 1):
            yield (item, *rest)


def _orbital_set(entries: list[dict], nbf: int) -> OrbitalSet:
    coefficients = np.zeros((len(entries), nbf), dtype=np.float32)
    for row, entry in zip(coefficients, entries):
        row[entry["indices"]] = entry["values"]
    return OrbitalSet(
        energies=np.array([e["energy"] for e in entries]),
        occupations=np.array([e["occupation"] for e in entries], dtype=np.float32),
        coefficients=coefficients,
        symmetries=[e["sym"] for e in entries])


def _orbitals(sections, centers: dict[int, np.ndarray]) -> Orbitals | None:
    if "gto" not in sections or "mo" not in sections or not centers:
        return None
    raw = _read_gto(_lines(sections, "gto"))
    entries = _read_mo(sections["mo"][1])
    if not raw or not entries:
        return None
    lowest = min(int(e["indices"].min()) for e in entries)
    if lowest < 0:
        raise ValueError("a coefficient is numbered 0; Molden counts from 1")
    nbf = max(int(e["indices"].max()) for e in entries) + 1
    pure = _counted_pure(raw, _pure_shells(sections), nbf)
    basis = _build_basis(raw, centers, pure, Convention())
    alpha = [e for e in entries if e["spin"] == "alpha"]
    beta = [e for e in entries if e["spin"] == "beta"]
    if not alpha:
        return None
    orbitals = Orbitals(basis=basis, alpha=_orbital_set(alpha, basis.nbf),
                        beta=_orbital_set(beta, basis.nbf) if beta else None)
    title = sections.get("title", ("", ""))[1].lower()
    _settle_convention(orbitals, raw, centers, pure, orca="orca" in title)
    return orbitals


EXACT = 1e-6           # as orthonormal as float32 coefficients get: stop looking


def _settle_convention(orbitals: Orbitals, raw: list[_RawShell], centers, pure,
                       orca: bool = False) -> None:
    """Read the basis the way under which the file's orbitals come out most
    nearly orthonormal. Only the choices this file's shells can tell apart are
    tried, in order of how often programs make them, stopping at the first
    reading that is exact; otherwise the best is kept, earlier ones winning a
    near-tie. A file whose title says ORCA wrote it has ORCA's reading tried
    first.
    """
    matters = (any(len(s.exponents) > 1 for s in raw),
               any(s.l >= 2 and s.l not in pure for s in raw),
               any(s.l >= 3 and s.l in pure for s in raw))
    scored: list[tuple[float, Convention, BasisSet]] = []
    for convention in ((_ORCA,) if orca else ()) + CONVENTIONS:
        convention = Convention(
            convention.normalized_primitives if matters[0] else True,
            convention.cartesian if matters[1] else "function",
            convention.orca_signs if matters[2] else False)
        if any(convention == c for _, c, _ in scored):
            continue
        basis = _build_basis(raw, centers, pure, convention)
        orbitals.basis = basis
        error = orbitals.reading_error()
        scored.append((error, convention, basis))
        if error < EXACT:
            break
    best = min(error for error, _, _ in scored)
    error, convention, basis = next(s for s in scored if s[0] <= TIE * best + 1e-12)
    orbitals.basis = basis
    orbitals.convention = convention.name if error < ORTHONORMAL else "unrecognized"
