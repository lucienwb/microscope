"""File loading/saving with format dispatch and optional cclib fallback."""

from __future__ import annotations

from pathlib import Path

from ..core.molecule import Molecule
from ..core.results import ParseResult
from . import cclib_bridge, cdxml, cube, fchk, gaussian, molden, molfile, orca, pdbfile, qchem, xyz
from .errors import FileFormatError, UnsupportedFormatError

OPEN_EXTENSIONS = (".xyz", ".log", ".out", ".fchk", ".fck", ".fch",
                   ".gjf", ".com", ".gau", ".pdb", ".molden", ".cube", ".cub")


def _sniff_output_program(path: Path) -> str:
    head = path.read_text(errors="replace")[:200_000]
    if gaussian.looks_like_output(head):
        return "gaussian"
    if orca.looks_like_output(head):
        return "orca"
    if qchem.looks_like_output(head):
        return "qchem"
    if molden.looks_like(head):
        return "molden"
    return "unknown"


def load(path) -> ParseResult:
    """Parse *path* with the native parsers; fall back to cclib if installed."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    ext = path.suffix.lower()

    try:
        if ext == ".xyz":
            return xyz.read(path)
        if ext in (".fchk", ".fck", ".fch"):
            return fchk.read(path)
        if ext in (".gjf", ".com", ".gau"):
            return gaussian.read_gjf(path)
        if ext == ".pdb":
            return pdbfile.read(path)
        if ext == ".molden":
            return molden.read(path)
        if ext in (".cube", ".cub"):
            return cube.read(path)
        if ext in (".log", ".out"):
            program = _sniff_output_program(path)
            if program == "gaussian":
                return gaussian.read_log(path)
            if program == "orca":
                return orca.read_log(path)
            if program == "qchem":
                return qchem.read_log(path)
            if program == "molden":
                return molden.read(path)
            result = _try_cclib(path)
            if result is not None:
                return result
            raise UnsupportedFormatError(
                f"{path.name}: could not recognize this output file; "
                "install cclib for extended format support."
            )
        raise UnsupportedFormatError(f"{path.name}: unknown file extension {ext!r}")
    except FileFormatError as native_error:
        result = _try_cclib(path)
        if result is not None:
            return result
        raise native_error


def _try_cclib(path: Path) -> ParseResult | None:
    if not cclib_bridge.available():
        return None
    try:
        return cclib_bridge.read(path)
    except Exception:
        return None


def save_molecule(path, molecule: Molecule) -> None:
    """Write *molecule* to *path*, format chosen by extension."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".xyz":
        xyz.write(path, molecule)
    elif ext in (".gjf", ".com", ".gau"):
        gaussian.write_gjf(path, molecule)
    elif ext == ".inp":
        orca.write_inp(path, molecule)
    elif ext in (".in", ".qcin"):
        qchem.write_in(path, molecule)
    elif ext == ".pdb":
        pdbfile.write(path, molecule)
    elif ext in (".mol", ".sdf"):
        molfile.write(path, molecule)
    else:
        raise UnsupportedFormatError(
            f"cannot write {ext!r} files (supported: .xyz, .gjf/.com, "
            ".inp [ORCA], .in/.qcin [Q-Chem], .pdb, .mol/.sdf)"
        )


LEWIS_EXTENSIONS = (".cdxml", ".mol", ".sdf")


def save_lewis(path, structure, coords2d, include=None) -> None:
    """Write a perceived Lewis structure as a file a drawing program can open.

    The coordinates are the ones on the page, so ChemDraw opens the structure
    at the angle it was composed at here. *include* lists the atoms that are
    actually drawn; hydrogens left out become implicit hydrogen counts.
    """
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".cdxml":
        cdxml.write(path, structure, coords2d, include)
    elif ext in (".mol", ".sdf"):
        molfile.write_2d(path, structure, coords2d, include)
    else:
        raise UnsupportedFormatError(
            f"cannot write {ext!r} structures "
            "(supported: .cdxml [ChemDraw], .mol/.sdf [MDL molfile])"
        )
