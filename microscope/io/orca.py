"""ORCA output (.out) reader."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..core.molecule import Molecule
from ..core.results import ExcitedState, NMRShielding, ParseResult, Vibration
from .errors import FileFormatError

WAVENUMBER_PER_EV = 8065.544

_RE_ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?[\d.]+)")
_RE_FREQ = re.compile(r"^\s*(\d+):\s+(-?[\d.]+)\s+cm\*\*-1")
_RE_IR_ROW = re.compile(r"^\s*(\d+):\s+(-?[\d.]+)\s+([\d.eE+-]+)\s+(-?[\d.]+)")
# ORCA >= 5: "0-1Ag -> 1-3Bu  3.129277  25239.3  396.2  0.000000000 ..."
_RE_TD_NEW = re.compile(
    r"->\s*(\S+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)")
# ORCA 4: "  1  25239.3  396.2  0.000000000 ..."
_RE_TD_OLD = re.compile(r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
_RE_NMR_ROW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+(-?[\d.]+)\s+(-?[\d.]+)")


def looks_like_output(text: str) -> bool:
    return "O   R   C   A" in text


def read_log(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    n = len(lines)

    frames_raw: list[tuple[list[str], np.ndarray]] = []
    charge = mult = None
    energies: list[float] = []
    freqs: dict[int, float] = {}
    ir_intensities: dict[int, float] = {}
    modes_matrix: np.ndarray | None = None
    excited: list[ExcitedState] = []
    nmr: list[NMRShielding] = []
    normal_term = False

    i = 0
    while i < n:
        line = lines[i]

        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            i += 2  # skip the dashes line
            syms, xyz = [], []
            while i < n:
                parts = lines[i].split()
                if len(parts) != 4:
                    break
                try:
                    coords3 = [float(p) for p in parts[1:4]]
                except ValueError:
                    break
                syms.append(parts[0])
                xyz.append(coords3)
                i += 1
            if syms:
                frames_raw.append((syms, np.array(xyz)))
            continue

        if charge is None and "Total Charge" in line and "...." in line:
            charge = int(line.split()[-1])
        elif mult is None and line.lstrip().startswith("Multiplicity") and "...." in line:
            mult = int(line.split()[-1])

        m = _RE_ENERGY.search(line)
        if m:
            energies.append(float(m.group(1)))
            i += 1
            continue

        if "VIBRATIONAL FREQUENCIES" in line:
            i += 1
            while i < n and "NORMAL MODES" not in lines[i]:
                fm = _RE_FREQ.match(lines[i])
                if fm:
                    freqs[int(fm.group(1))] = float(fm.group(2))
                i += 1
            continue

        if line.strip() == "NORMAL MODES":
            i, modes_matrix = _parse_normal_modes(lines, i + 1)
            continue

        if line.strip() == "IR SPECTRUM":
            i += 1
            started = False
            while i < n:
                row = _RE_IR_ROW.match(lines[i])
                if row:
                    started = True
                    ir_intensities[int(row.group(1))] = float(row.group(4))
                elif started:
                    break
                elif i < n and ("SPECTRUM" not in lines[i] or not started):
                    s = lines[i].strip()
                    if s and not (s.startswith(("Mode", "---", "cm**-1", "*"))):
                        if started:
                            break
                i += 1
            continue

        if ("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in line
                and "SOC" not in line):
            i += 1
            started = False
            while i < n:
                s = lines[i]
                m_new = _RE_TD_NEW.search(s)
                m_old = _RE_TD_OLD.match(s) if not m_new else None
                if m_new:
                    started = True
                    ev = float(m_new.group(2))
                    excited.append(ExcitedState(
                        index=len(excited) + 1, label=m_new.group(1),
                        energy_ev=ev, wavelength_nm=float(m_new.group(4)),
                        osc_strength=float(m_new.group(5))))
                elif m_old and started:
                    cm = float(m_old.group(2))
                    excited.append(ExcitedState(
                        index=len(excited) + 1, label="",
                        energy_ev=cm / WAVENUMBER_PER_EV,
                        wavelength_nm=float(m_old.group(3)),
                        osc_strength=float(m_old.group(4))))
                elif m_old and not started and s.strip() and s.split()[0].rstrip(":").isdigit():
                    started = True
                    cm = float(m_old.group(2))
                    excited.append(ExcitedState(
                        index=len(excited) + 1, label="",
                        energy_ev=cm / WAVENUMBER_PER_EV,
                        wavelength_nm=float(m_old.group(3)),
                        osc_strength=float(m_old.group(4))))
                elif started:
                    break
                elif not s.strip() and started:
                    break
                i += 1
            continue

        if "CHEMICAL SHIELDING SUMMARY" in line:
            nmr = []
            i += 1
            started = False
            while i < n:
                m_row = _RE_NMR_ROW.match(lines[i])
                if m_row:
                    started = True
                    nmr.append(NMRShielding(
                        atom_index=int(m_row.group(1)),
                        symbol=m_row.group(2),
                        isotropic=float(m_row.group(3))))
                elif started:
                    break
                i += 1
            continue

        if "ORCA TERMINATED NORMALLY" in line:
            normal_term = True

        i += 1

    if not frames_raw:
        raise FileFormatError(f"{path}: no geometry found in ORCA output")

    frames: list[Molecule] = []
    for syms, xyz in frames_raw:
        if frames and len(frames[-1].symbols) == len(syms) and np.allclose(frames[-1].coords, xyz):
            continue
        frames.append(Molecule(syms, xyz, charge=charge or 0, multiplicity=mult or 1,
                               charge_known=charge is not None))

    natoms = frames[-1].natoms
    vibrations: list[Vibration] = []
    for idx in sorted(freqs):
        f = freqs[idx]
        if abs(f) < 1e-3:
            continue  # translations/rotations
        disp = None
        if (modes_matrix is not None and idx < modes_matrix.shape[1]
                and modes_matrix.shape[0] == 3 * natoms):
            disp = modes_matrix[:, idx].reshape(natoms, 3)
        vibrations.append(Vibration(
            frequency=f, ir_intensity=ir_intensities.get(idx),
            displacements=disp))

    return ParseResult(
        frames=frames, program="ORCA", source=str(path),
        scf_energies=energies, vibrations=vibrations,
        excited_states=excited, nmr_shieldings=nmr,
        normal_termination=normal_term)


def write_inp(path, molecule: Molecule,
              keywords: str = "B3LYP D4 def2-SVP Opt Freq",
              nprocs: int = 8, maxcore_mb: int = 4000) -> None:
    """Write an ORCA input file with a Cartesian * xyz block."""
    kw = keywords.strip()
    if not kw.startswith("!"):
        kw = "! " + kw
    with open(path, "w") as fh:
        fh.write(kw + "\n")
        fh.write(f"%pal nprocs {nprocs} end\n")
        fh.write(f"%maxcore {maxcore_mb}\n\n")
        fh.write(f"* xyz {molecule.charge} {molecule.multiplicity}\n")
        for sym, (x, y, z) in zip(molecule.symbols, molecule.coords):
            fh.write(f" {sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        fh.write("*\n")


def _parse_normal_modes(lines: list[str], i: int) -> tuple[int, np.ndarray | None]:
    """Parse the ORCA NORMAL MODES matrix (3N x 3N, printed in column blocks)."""
    n = len(lines)
    data: dict[tuple[int, int], float] = {}
    current_cols: list[int] | None = None
    started = False
    expect_irrep = False
    while i < n:
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        tokens = s.split()
        if all(re.fullmatch(r"\d+", t) for t in tokens):
            current_cols = [int(t) for t in tokens]
            started = True
            expect_irrep = True
            i += 1
            continue
        if current_cols is not None and re.fullmatch(r"\d+", tokens[0]) \
                and len(tokens) == len(current_cols) + 1:
            try:
                vals = [float(t) for t in tokens[1:]]
            except ValueError:
                i += 1
                continue
            r = int(tokens[0])
            for c, v in zip(current_cols, vals):
                data[(r, c)] = v
            expect_irrep = False
            i += 1
            continue
        if not started:
            i += 1  # descriptive text before the matrix
            continue
        if expect_irrep and all("." not in t for t in tokens):
            # irrep labels line right below a column header (e.g. "1-Au  1-Bu ...")
            expect_irrep = False
            i += 1
            continue
        break  # anything else ends the section (e.g. the IR SPECTRUM header)
    if not data:
        return i, None
    nrows = max(r for r, _ in data) + 1
    ncols = max(c for _, c in data) + 1
    mat = np.zeros((nrows, ncols))
    for (r, c), v in data.items():
        mat[r, c] = v
    return i, mat
