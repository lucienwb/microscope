"""Gaussian cube file reader (.cube/.cub) — volumetric MO/density grids.

Format: two comment lines; a header with the atom count, grid origin, and
three axis lines (point count + step vector, positive counts = Bohr); the
atom list (Z, charge, x, y, z); and the scalar field with z varying fastest.
A negative atom count marks an orbital cube whose extra header line lists the
MO numbers stored per grid point (cubegen's multi-MO cubes).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core.elements import normalize_symbol
from ..core.molecule import Molecule
from ..core.results import ParseResult
from ..core.volume import VolumeData
from .errors import FileFormatError
from .fchk import BOHR_TO_ANGSTROM


def read(path) -> ParseResult:
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 7:
        raise FileFormatError(f"{path.name}: too short to be a cube file")
    title = lines[0].strip()
    subtitle = lines[1].strip()
    try:
        parts = lines[2].split()
        natoms_signed = int(parts[0])
        origin = np.array([float(x) for x in parts[1:4]])
        nval = int(parts[4]) if len(parts) > 4 else 1
        npts, axes = [], []
        for k in range(3):
            p = lines[3 + k].split()
            npts.append(int(p[0]))
            axes.append([float(x) for x in p[1:4]])
    except (ValueError, IndexError) as exc:
        raise FileFormatError(f"{path.name}: malformed cube header") from exc

    natoms = abs(natoms_signed)
    shape = tuple(abs(n) for n in npts)
    if natoms == 0 or min(shape) < 2:
        raise FileFormatError(f"{path.name}: empty cube grid or atom list")
    # positive point counts mean atomic units (the overwhelmingly common case)
    scale = BOHR_TO_ANGSTROM if npts[0] > 0 else 1.0
    origin *= scale
    axes = np.array(axes) * scale

    try:
        symbols, coords = [], []
        for ln in lines[6:6 + natoms]:
            p = ln.split()
            symbols.append(normalize_symbol(int(float(p[0]))))
            coords.append([float(x) for x in p[2:5]])
        coords = np.array(coords) * scale

        idx = 6 + natoms
        mo_ids: list[int] = []
        if natoms_signed < 0:            # orbital cube: DSET id line(s) follow
            toks = lines[idx].split()
            idx += 1
            count = int(toks[0])
            mo_ids = [int(t) for t in toks[1:]]
            while len(mo_ids) < count:
                mo_ids.extend(int(t) for t in lines[idx].split())
                idx += 1
            nval = count
        data = " ".join(lines[idx:]).split()
        expected = shape[0] * shape[1] * shape[2] * nval
        if len(data) < expected:
            raise FileFormatError(
                f"{path.name}: truncated cube data "
                f"({len(data)} of {expected} values)")
        # single precision: a cube file carries five significant digits,
        # and a 200^3 grid is 64 MB in float64 against 32 in float32
        grid = np.array(data[:expected], dtype=np.float32).reshape(*shape, nval)
    except (ValueError, IndexError) as exc:
        raise FileFormatError(f"{path.name}: malformed cube data") from exc

    molecule = Molecule(symbols, coords, title=title or path.stem,
                        charge_known=False)
    base_label = subtitle or title or path.stem
    volumes = []
    for k in range(nval):
        if mo_ids:
            label = f"MO {mo_ids[k]}"
        elif nval > 1:
            label = f"{base_label} [{k + 1}]"
        else:
            label = base_label
        volumes.append(VolumeData(origin=origin.copy(), axes=axes.copy(),
                                  values=np.ascontiguousarray(grid[..., k]),
                                  label=label))
    return ParseResult(frames=[molecule], program="cube", source=str(path),
                       volumes=volumes)


def write(path, volume: VolumeData, molecule: Molecule, comment: str = "") -> None:
    """Write a grid as a Gaussian cube file, as cubegen lays one out: bohr,
    z running fastest, six values to a line and a new line for each z-row.
    VMD, Avogadro, GaussView and the rest read it; so does :func:`read`."""
    to_bohr = 1.0 / BOHR_TO_ANGSTROM
    values = np.asarray(volume.values, dtype=np.float32)
    nx, ny, nz = values.shape
    out = [_ascii(molecule.title or "microscope"),
           _ascii(comment or volume.label or "grid"),
           _cube_row(molecule.natoms, volume.origin * to_bohr)]
    for n, step in zip(values.shape, volume.axes * to_bohr):
        out.append(_cube_row(n, step))
    for z, xyz in zip(molecule.atomic_numbers, molecule.coords * to_bohr):
        out.append(f"{int(z):5d}{float(z):12.6f}{xyz[0]:12.6f}{xyz[1]:12.6f}{xyz[2]:12.6f}")
    row = ("%13.5E" * 6 + "\n") * (nz // 6) + ("%13.5E" * (nz % 6) + "\n" if nz % 6 else "")
    body = "".join(row % tuple(r) for r in values.reshape(nx * ny, nz).tolist())
    Path(path).write_text("\n".join(out) + "\n" + body, encoding="ascii")


def _ascii(text: str) -> str:
    """Header lines stay ASCII: older readers are Fortran, and a Windows
    locale cannot even encode an alpha."""
    text = " ".join(text.split())                  # one line, whatever the title held
    text = text.replace(" · ", ", ").replace("α", "alpha").replace("β", "beta")
    return text.encode("ascii", "replace").decode().strip()[:80]


def _cube_row(n: int, vector) -> str:
    return f"{n:5d}{vector[0]:12.6f}{vector[1]:12.6f}{vector[2]:12.6f}"
