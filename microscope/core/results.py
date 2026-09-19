"""Containers for parsed calculation results."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .molecule import Molecule
from .orbitals import Orbitals


@dataclass
class Vibration:
    frequency: float                       # cm^-1 (negative = imaginary)
    ir_intensity: float | None = None      # km/mol
    raman_activity: float | None = None
    displacements: np.ndarray | None = None  # (natoms, 3)
    symmetry: str = ""


@dataclass
class ExcitedState:
    index: int
    label: str
    energy_ev: float
    wavelength_nm: float
    osc_strength: float


@dataclass
class NMRShielding:
    atom_index: int                        # 0-based
    symbol: str
    isotropic: float                       # ppm


@dataclass
class ParseResult:
    """Everything read from one file. A field the format does not carry is
    left empty rather than guessed.

    Attributes:
        frames: The geometries, in file order; ``molecule`` is the last.
        program: Which reader produced this ("Gaussian", "ORCA", "Molden", ...).
        source: The path it was read from.
        scf_energies: SCF energy per step, in hartree.
        vibrations: Normal modes with frequencies and IR intensities.
        excited_states: TD-DFT / CIS excitations.
        nmr_shieldings: Isotropic shieldings per atom.
        volumes: Grids from a cube file, ready for an isosurface.
        orbitals: The basis set and molecular orbitals of an fchk or Molden
            file; ``orbitals.volume("homo")`` puts one on a grid.
        normal_termination: Whether the job said it finished, where it says.
        route: The job's route or keyword line, where it has one.
        parser: "native", or "cclib" when the optional fallback read the file.
        warnings: What could not be read although the rest was, e.g. a
            damaged wavefunction in a file whose geometry is fine.
    """

    frames: list[Molecule] = field(default_factory=list)
    program: str = ""
    source: str = ""
    scf_energies: list[float] = field(default_factory=list)  # Hartree
    vibrations: list[Vibration] = field(default_factory=list)
    excited_states: list[ExcitedState] = field(default_factory=list)
    nmr_shieldings: list[NMRShielding] = field(default_factory=list)
    volumes: list = field(default_factory=list)              # VolumeData grids
    orbitals: Orbitals | None = None     # basis set + MOs, from fchk / Molden files
    normal_termination: bool | None = None
    route: str = ""              # the job's route/keyword line, where it has one
    parser: str = "native"       # "native", or "cclib" when the fallback read it
    warnings: list[str] = field(default_factory=list)   # read, but not all of it

    @property
    def molecule(self) -> Molecule | None:
        return self.frames[-1] if self.frames else None

    @property
    def nframes(self) -> int:
        return len(self.frames)

    @property
    def read_by_cclib(self) -> bool:
        """Did the optional fallback read this, rather than a parser here?"""
        return self.parser == "cclib"
