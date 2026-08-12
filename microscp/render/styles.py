"""Rendering style presets."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core import elements

# CYLview-like palette overrides on top of CPK colors.
CYLVIEW_COLORS = {
    1: (1.00, 1.00, 1.00),    # H  white
    6: (0.33, 0.33, 0.33),    # C  dark gray
    7: (0.20, 0.33, 0.95),    # N  blue
    8: (0.88, 0.12, 0.10),    # O  red
    9: (0.55, 0.85, 0.35),    # F  light green
    15: (1.00, 0.60, 0.10),   # P  orange
    16: (0.95, 0.85, 0.15),   # S  yellow
    17: (0.25, 0.80, 0.30),   # Cl green
    35: (0.65, 0.28, 0.17),   # Br dark red-brown
    53: (0.58, 0.15, 0.60),   # I  purple
}


@dataclass
class Style:
    name: str = "cylview"
    bond_radius: float = 0.165           # Angstrom
    atom_scale: float = 0.36             # fraction of covalent radius
    min_atom_radius: float = 0.21        # keeps H spheres visible
    palette: dict = field(default_factory=lambda: dict(CYLVIEW_COLORS))
    background: tuple = (1.0, 1.0, 1.0)
    show_hbonds: bool = True
    hbond_radius: float = 0.055
    hbond_color: tuple = (0.45, 0.45, 0.45)
    hbond_dash: float = 0.20             # dash length, Angstrom
    hbond_gap: float = 0.14              # gap between dashes

    def atom_radius(self, z: int) -> float:
        return max(elements.covalent_radius(z) * self.atom_scale, self.min_atom_radius)

    def atom_color(self, z: int) -> tuple[float, float, float]:
        return self.palette.get(int(z), elements.cpk_color(z))
