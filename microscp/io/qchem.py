"""Q-Chem output (.out) reader."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..core.molecule import Molecule
from ..core.results import ExcitedState, ParseResult, Vibration
from .errors import FileFormatError

_RE_ENERGY = re.compile(r"Total energy in the final basis set\s*=\s*(-?[\d.]+)")
_RE_EXC = re.compile(r"Excited state\s+(\d+): excitation energy \(eV\)\s*=\s*(-?[\d.]+)")
_RE_STRENGTH = re.compile(r"Strength\s*:\s*([\d.eE+-]+)")
_RE_MULT = re.compile(r"Multiplicity:\s*(\S+)")


def looks_like_output(text: str) -> bool:
    return "Welcome to Q-Chem" in text or "Q-Chem, Inc" in text


def read_log(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    n = len(lines)

    frames_raw: list[tuple[list[str], np.ndarray]] = []
    charge = mult = None
    energies: list[float] = []
    vibrations: list[Vibration] = []
    excited: list[ExcitedState] = []
    normal_term = False

    i = 0
    while i < n:
        line = lines[i]

        if "Standard Nuclear Orientation" in line:
            i += 3  # column header + dashes
            syms, xyz = [], []
            while i < n and not lines[i].lstrip().startswith("----"):
                parts = lines[i].split()
                if len(parts) < 5:
                    break
                try:
                    coords3 = [float(p) for p in parts[2:5]]
                except ValueError:
                    break
                syms.append(parts[1])
                xyz.append(coords3)
                i += 1
            if syms:
                frames_raw.append((syms, np.array(xyz)))
            continue

        if charge is None and line.strip() == "$molecule" and i + 1 < n:
            parts = lines[i + 1].split()
            if len(parts) >= 2:
                try:
                    charge, mult = int(parts[0]), int(parts[1])
                except ValueError:
                    pass
            i += 2
            continue

        m = _RE_ENERGY.search(line)
        if m:
            energies.append(float(m.group(1)))
            i += 1
            continue

        if line.lstrip().startswith("Mode:"):
            i = _parse_freq_block(lines, i, vibrations)
            continue

        m = _RE_EXC.search(line)
        if m:
            ev = float(m.group(2))
            excited.append(ExcitedState(
                index=int(m.group(1)), label="", energy_ev=ev,
                wavelength_nm=1239.841984 / ev if ev > 0 else 0.0,
                osc_strength=0.0))
            i += 1
            continue
        if excited:
            m = _RE_MULT.search(line)
            if m and not excited[-1].label:
                excited[-1].label = m.group(1)
            m = _RE_STRENGTH.search(line)
            if m:
                excited[-1].osc_strength = float(m.group(1))

        if "Thank you very much for using Q-Chem" in line:
            normal_term = True

        i += 1

    if not frames_raw:
        raise FileFormatError(f"{path}: no geometry found in Q-Chem output")

    frames: list[Molecule] = []
    for syms, xyz in frames_raw:
        if frames and len(frames[-1].symbols) == len(syms) and np.allclose(frames[-1].coords, xyz):
            continue
        frames.append(Molecule(syms, xyz, charge=charge or 0, multiplicity=mult or 1))

    return ParseResult(
        frames=frames, program="Q-Chem", source=str(path),
        scf_energies=energies, vibrations=vibrations,
        excited_states=excited, normal_termination=normal_term)


def _parse_freq_block(lines: list[str], i: int, vibrations: list[Vibration]) -> int:
    """Parse one 'Mode: a b c' column group; return next line index."""
    n = len(lines)
    nmodes = len(lines[i].split()) - 1
    freqs = irs = ramans = None
    disp_rows: list[list[float]] = []
    j = i + 1
    while j < n:
        s = lines[j].strip()
        if s.startswith("Frequency:"):
            freqs = [float(t) for t in s.split()[1:]]
        elif s.startswith("IR Intens:"):
            irs = [float(t) for t in s.split()[2:]]
        elif s.startswith("Raman Intens:"):
            ramans = [float(t) for t in s.split()[2:]]
        elif s.startswith("X") and "Y" in s.split() and "Z" in s.split():
            j += 1
            while j < n:
                parts = lines[j].split()
                # displacement rows start with an element symbol ("TransDip" etc. do not)
                if (len(parts) == 1 + 3 * nmodes and parts[0].isalpha()
                        and len(parts[0]) <= 2):
                    try:
                        disp_rows.append([float(x) for x in parts[1:]])
                    except ValueError:
                        break
                    j += 1
                else:
                    break
            break
        elif s.startswith("Mode:") or s.startswith("STANDARD THERMODYNAMIC"):
            break
        j += 1
    if freqs is None:
        return j
    disp = np.array(disp_rows) if disp_rows else None
    for k in range(len(freqs)):
        d = disp[:, 3 * k:3 * k + 3] if disp is not None else None
        vibrations.append(Vibration(
            frequency=freqs[k],
            ir_intensity=irs[k] if irs and k < len(irs) else None,
            raman_activity=ramans[k] if ramans and k < len(ramans) else None,
            displacements=d))
    return j
