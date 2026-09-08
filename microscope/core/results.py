"""Containers for parsed calculation results."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .molecule import Molecule


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
    frames: list[Molecule] = field(default_factory=list)
    program: str = ""
    source: str = ""
    scf_energies: list[float] = field(default_factory=list)  # Hartree
    vibrations: list[Vibration] = field(default_factory=list)
    excited_states: list[ExcitedState] = field(default_factory=list)
    nmr_shieldings: list[NMRShielding] = field(default_factory=list)
    volumes: list = field(default_factory=list)              # VolumeData grids
    normal_termination: bool | None = None
    route: str = ""              # the job's route/keyword line, where it has one
    parser: str = "native"       # "native", or "cclib" when the fallback read it

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
