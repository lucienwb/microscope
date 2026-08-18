"""Molden format reader ([Atoms], [GEOMETRIES], [FREQ]/[FR-NORM-COORD])."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..core import elements
from ..core.molecule import Molecule
from ..core.results import ParseResult, Vibration
from .errors import FileFormatError

BOHR_TO_ANGSTROM = 0.52917721092

_SECTION_RE = re.compile(r"^\s*\[([A-Za-z0-9_ -]+)\]\s*(.*)$")


def looks_like(text: str) -> bool:
    return "[molden format]" in text.lower()


def _split_sections(lines: list[str]) -> dict[str, tuple[str, list[str]]]:
    sections: dict[str, tuple[str, list[str]]] = {}
    name, arg, body = None, "", []
    for line in lines:
        m = _SECTION_RE.match(line)
        if m:
            if name is not None:
                sections[name] = (arg, body)
            name = m.group(1).strip().lower()
            arg = m.group(2).strip()
            body = []
        elif name is not None:
            body.append(line)
    if name is not None:
        sections[name] = (arg, body)
    return sections


def read(path) -> ParseResult:
    lines = Path(path).read_text(errors="replace").splitlines()
    sections = _split_sections(lines)
    if not sections:
        raise FileFormatError(f"{path}: not a Molden file")

    frames: list[Molecule] = []

    # primary geometry from [Atoms]
    natoms = 0
    if "atoms" in sections:
        arg, body = sections["atoms"]
        to_ang = BOHR_TO_ANGSTROM if arg.lower().startswith("au") else 1.0
        syms, xyz = [], []
        for row in body:
            parts = row.split()
            if len(parts) < 6:
                continue
            try:
                z = int(parts[2])
                coords3 = [float(v) * to_ang for v in parts[3:6]]
            except ValueError:
                continue
            syms.append(elements.SYMBOLS[z] if 0 < z < len(elements.SYMBOLS) else parts[0])
            xyz.append(coords3)
        if syms:
            frames.append(Molecule(syms, np.array(xyz)))
            natoms = len(syms)

    # optimization history from [GEOMETRIES] XYZ
    if "geometries" in sections:
        arg, body = sections["geometries"]
        if arg.upper().startswith("XYZ") or not arg:
            geo_frames = _parse_xyz_frames(body)
            if geo_frames:
                frames = geo_frames  # full trajectory supersedes the single frame
                natoms = frames[-1].natoms

    # frequency job geometry from [FR-COORD] (bohr) if nothing else
    if not frames and "fr-coord" in sections:
        _, body = sections["fr-coord"]
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
            frames.append(Molecule(syms, np.array(xyz)))
            natoms = len(syms)

    if not frames:
        raise FileFormatError(f"{path}: no geometry found in Molden file")

    vibrations: list[Vibration] = []
    if "freq" in sections:
        freqs = [float(t) for row in sections["freq"][1] for t in row.split()]
        intensities = []
        if "int" in sections:
            intensities = [float(t) for row in sections["int"][1] for t in row.split()]
        disps: list[np.ndarray] = []
        if "fr-norm-coord" in sections and natoms:
            current: list[list[float]] = []
            for row in sections["fr-norm-coord"][1]:
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

    return ParseResult(frames=frames, program="Molden", source=str(path),
                       vibrations=vibrations)


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
            frames.append(Molecule(syms, np.array(xyz)))
        i += 2 + count
    return frames
