"""Gaussian output (.log/.out) and input (.gjf/.com) reader/writer."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..core import elements
from ..core.molecule import Molecule
from ..core.results import ExcitedState, NMRShielding, ParseResult, Vibration
from .errors import FileFormatError

_RE_SCF = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+(-?[\d.]+(?:[DEde][+-]?\d+)?)")
_RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
_RE_EXC = re.compile(
    r"Excited State\s+(\d+):\s+(\S+)\s+(-?[\d.]+)\s*eV\s+(-?[\d.]+)\s*nm\s+f=(-?[\d.]+)"
)
_RE_NMR = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+Isotropic\s+=\s+(-?[\d.]+)")


def looks_like_output(text: str) -> bool:
    return "Entering Gaussian System" in text or "Gaussian, Inc" in text


def read_log(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    n = len(lines)

    frames_std: list[tuple[list[int], np.ndarray]] = []
    frames_inp: list[tuple[list[int], np.ndarray]] = []
    charge, mult = None, None
    route = None
    scf_energies: list[float] = []
    vibrations: list[Vibration] = []
    excited: dict[int, ExcitedState] = {}
    nmr: list[NMRShielding] = []
    normal_terms = 0

    i = 0
    while i < n:
        line = lines[i]

        # "Z-Matrix orientation" is what older Gaussian writes for a job whose
        # input was a Z-matrix; the block below it is laid out like the others
        if "orientation:" in line and ("Standard" in line or "Input" in line
                                       or "Z-Matrix" in line):
            target = frames_std if "Standard" in line else frames_inp
            i += 5  # skip dashes, two header lines, dashes
            zs, xyz = [], []
            while i < n and not lines[i].lstrip().startswith("---"):
                parts = lines[i].split()
                if len(parts) < 5:
                    break
                z = int(parts[1])
                if z > 0:  # skip dummy/ghost centers
                    zs.append(z)
                    xyz.append([float(v) for v in parts[-3:]])
                i += 1
            if zs:
                target.append((zs, np.array(xyz)))
            continue

        if charge is None:
            m = _RE_CHARGE.search(line)
            if m:
                charge, mult = int(m.group(1)), int(m.group(2))
                i += 1
                continue

        if route is None and line.lstrip().startswith("#"):
            route = line.strip()
            i += 1
            continue

        m = _RE_SCF.search(line)
        if m:
            scf_energies.append(float(m.group(1).replace("D", "E").replace("d", "e")))
            i += 1
            continue

        stripped = line.lstrip()
        if stripped.startswith("Frequencies --") and not stripped.startswith("Frequencies ---"):
            i = _parse_freq_block(lines, i, vibrations)
            continue

        m = _RE_EXC.search(line)
        if m:
            idx = int(m.group(1))
            excited[idx] = ExcitedState(
                index=idx,
                label=m.group(2),
                energy_ev=float(m.group(3)),
                wavelength_nm=float(m.group(4)),
                osc_strength=float(m.group(5)),
            )
            i += 1
            continue

        if "Magnetic shielding tensor (ppm)" in line:
            nmr = []  # keep only the last NMR block in the file
            i += 1
            while i < n:
                m = _RE_NMR.match(lines[i])
                if m:
                    nmr.append(NMRShielding(
                        atom_index=int(m.group(1)) - 1,
                        symbol=elements.normalize_symbol(m.group(2)),
                        isotropic=float(m.group(3)),
                    ))
                elif re.match(r"^\s*([XYZ][XYZ]=|Eigenvalues)", lines[i]):
                    pass  # shielding tensor components between atoms
                else:
                    break
                i += 1
            continue

        if "Normal termination" in line:
            normal_terms += 1

        i += 1

    raw_frames = frames_std or frames_inp
    if not raw_frames:
        raise FileFormatError(f"{path}: no geometry found in Gaussian output")

    frames: list[Molecule] = []
    for zs, xyz in raw_frames:
        if frames and len(frames[-1].symbols) == len(zs) and np.allclose(frames[-1].coords, xyz):
            continue  # freq/summary re-prints of the same geometry
        frames.append(Molecule(
            symbols=[elements.SYMBOLS[z] if z < len(elements.SYMBOLS) else "X" for z in zs],
            coords=xyz,
            charge=charge or 0,
            multiplicity=mult or 1,
            charge_known=charge is not None,
        ))

    result = ParseResult(
        frames=frames,
        program="Gaussian",
        source=str(path),
        scf_energies=scf_energies,
        vibrations=vibrations,
        excited_states=[excited[k] for k in sorted(excited)],
        nmr_shieldings=nmr,
        normal_termination=normal_terms > 0,
    )
    if route:
        result.route = route
    return result


def _parse_freq_block(lines: list[str], i: int, vibrations: list[Vibration]) -> int:
    """Parse one 'Frequencies --' column group; return the next line index."""
    n = len(lines)
    freqs = [float(x) for x in lines[i].split("--")[1].split()]
    nmodes = len(freqs)
    ir = raman = None
    syms = []
    if i >= 1:
        prev = lines[i - 1].split()
        if prev and all(not p.replace(".", "").replace("-", "").isdigit() for p in prev):
            syms = prev if len(prev) == nmodes else []
    disp_rows: list[list[float]] = []
    j = i + 1
    while j < n:
        s = lines[j].lstrip()
        if s.startswith("IR Inten"):
            ir = [float(x) for x in lines[j].split("--")[1].split()]
        elif s.startswith("Raman Activ"):
            raman = [float(x) for x in lines[j].split("--")[1].split()]
        elif s.startswith("Atom") and "AN" in s:
            j += 1
            while j < n:
                parts = lines[j].split()
                if len(parts) >= 2 + 3 * nmodes and parts[0].isdigit():
                    disp_rows.append([float(x) for x in parts[2:2 + 3 * nmodes]])
                    j += 1
                else:
                    break
            break
        elif s.startswith("Frequencies"):
            break
        j += 1
    disp = np.array(disp_rows) if disp_rows else None
    for k in range(nmodes):
        d = disp[:, 3 * k:3 * k + 3] if disp is not None else None
        vibrations.append(Vibration(
            frequency=freqs[k],
            ir_intensity=ir[k] if ir and k < len(ir) else None,
            raman_activity=raman[k] if raman and k < len(raman) else None,
            displacements=d,
            symmetry=syms[k] if k < len(syms) else "",
        ))
    return j


def read_gjf(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    n = len(lines)
    i = 0
    while i < n and lines[i].strip().startswith("%"):
        i += 1
    route_lines = []
    while i < n and lines[i].strip():
        route_lines.append(lines[i].strip())
        i += 1
    i += 1
    route_text = " ".join(route_lines).lower()
    if "allcheck" in route_text or "geom=check" in route_text:
        raise FileFormatError(
            f"{path}: geometry is read from a checkpoint file "
            "(geom=check/allcheck); this input contains no coordinates"
        )
    title_lines = []
    while i < n and lines[i].strip():
        title_lines.append(lines[i].strip())
        i += 1
    i += 1
    if i >= n:
        raise FileFormatError(f"{path}: truncated Gaussian input")
    parts = lines[i].split()
    try:
        charge, mult = int(parts[0]), int(parts[1])
    except (ValueError, IndexError) as exc:
        raise FileFormatError(f"{path}: bad charge/multiplicity line {lines[i]!r}") from exc
    i += 1
    symbols, coords = [], []
    while i < n and lines[i].strip():
        parts = lines[i].split()
        if len(parts) < 4:
            break
        m = re.match(r"[A-Za-z]{1,2}|\d{1,3}", parts[0])
        if not m:
            raise FileFormatError(f"{path}: cannot read atom line {lines[i]!r}")
        sym = m.group(0)
        vals = parts[1:]
        # optional freeze-code integer column after the element
        if len(vals) >= 4 and "." not in vals[0] and vals[0].lstrip("-").isdigit():
            vals = vals[1:]
        try:
            xyz = [float(v) for v in vals[:3]]
        except ValueError as exc:
            raise FileFormatError(
                f"{path}: only Cartesian input is supported (Z-matrix line {lines[i]!r})"
            ) from exc
        symbols.append(sym)
        coords.append(xyz)
        i += 1
    if not symbols:
        raise FileFormatError(f"{path}: no atoms found")
    mol = Molecule(symbols, np.array(coords), charge=charge, multiplicity=mult,
                   title=" ".join(title_lines))
    result = ParseResult(frames=[mol], program="Gaussian input", source=str(path))
    if route_lines:
        result.route = " ".join(route_lines)
    return result


def write_gjf(path, molecule: Molecule, route: str = "#p B3LYP/6-31G(d) Opt Freq",
              link0: list[str] | None = None) -> None:
    stem = Path(path).stem
    with open(path, "w") as fh:
        for line in link0 or [f"%chk={stem}.chk", "%mem=4GB", "%nprocshared=8"]:
            fh.write(line + "\n")
        fh.write(route + "\n\n")
        fh.write((molecule.title or stem) + "\n\n")
        fh.write(f"{molecule.charge} {molecule.multiplicity}\n")
        for sym, (x, y, z) in zip(molecule.symbols, molecule.coords):
            fh.write(f" {sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        fh.write("\n")
