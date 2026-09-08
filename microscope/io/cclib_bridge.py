"""Optional fallback parser backend using cclib, if the user has it installed."""

from __future__ import annotations

import numpy as np

from ..core import elements
from ..core.molecule import Molecule
from ..core.results import ExcitedState, NMRShielding, ParseResult, Vibration

WAVENUMBER_PER_EV = 8065.544


def available() -> bool:
    try:
        import cclib  # noqa: F401
        return True
    except ImportError:
        return False


def read(path) -> ParseResult | None:
    """Parse *path* with cclib. Returns None if cclib cannot handle it."""
    import cclib

    data = cclib.io.ccread(str(path))
    if data is None or not hasattr(data, "atomnos"):
        return None

    symbols = [elements.SYMBOLS[z] if 0 < z < len(elements.SYMBOLS) else "X"
               for z in data.atomnos]
    charge = getattr(data, "charge", 0)
    mult = getattr(data, "mult", 1)

    all_coords = getattr(data, "atomcoords", None)
    if all_coords is None or len(all_coords) == 0:
        return None
    frames = [Molecule(symbols, np.array(c), charge=charge, multiplicity=mult)
              for c in all_coords]

    vibrations = []
    freqs = getattr(data, "vibfreqs", None)
    if freqs is not None:
        irs = getattr(data, "vibirs", None)
        disps = getattr(data, "vibdisps", None)
        syms = getattr(data, "vibsyms", None)
        for k, f in enumerate(freqs):
            vibrations.append(Vibration(
                frequency=float(f),
                ir_intensity=float(irs[k]) if irs is not None else None,
                displacements=np.array(disps[k]) if disps is not None else None,
                symmetry=str(syms[k]) if syms is not None else "",
            ))

    excited = []
    etenergies = getattr(data, "etenergies", None)
    if etenergies is not None:
        etoscs = getattr(data, "etoscs", None)
        etsyms = getattr(data, "etsyms", None)
        for k, e_cm in enumerate(etenergies):
            e_ev = float(e_cm) / WAVENUMBER_PER_EV
            excited.append(ExcitedState(
                index=k + 1,
                label=str(etsyms[k]) if etsyms is not None else "",
                energy_ev=e_ev,
                wavelength_nm=1239.841984 / e_ev if e_ev > 0 else 0.0,
                osc_strength=float(etoscs[k]) if etoscs is not None else 0.0,
            ))

    nmr = []
    tensors = getattr(data, "nmrtensors", None)
    if tensors:
        for atom_idx, tdict in sorted(tensors.items()):
            total = tdict.get("total") if isinstance(tdict, dict) else None
            if total is not None:
                nmr.append(NMRShielding(
                    atom_index=int(atom_idx),
                    symbol=symbols[int(atom_idx)],
                    isotropic=float(np.trace(np.array(total)) / 3.0),
                ))

    scf = [float(e) for e in getattr(data, "scfenergies", [])]  # eV in cclib
    metadata = getattr(data, "metadata", {}) or {}

    return ParseResult(
        frames=frames,
        program=str(metadata.get("package", "unknown")) + " (via cclib)",
        source=str(path),
        scf_energies=scf,
        vibrations=vibrations,
        excited_states=excited,
        nmr_shieldings=nmr,
        normal_termination=bool(metadata.get("success")) if "success" in metadata else None,
        parser="cclib",
    )
