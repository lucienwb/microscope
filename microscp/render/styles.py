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
    surface_positive: tuple = (0.29, 0.44, 0.86)   # isosurface + lobe (blue)
    surface_negative: tuple = (0.88, 0.38, 0.22)   # isosurface - lobe (orange-red)
    surface_opacity: float = 0.62
    bond_color: tuple | None = None      # uniform bonds (None = split by atom)
    quadrant_color: tuple | None = None  # seam lines on heavy atoms (Houkmol)
    quadrant_width: float = 0.055        # seam line width (fraction of radius)

    def atom_radius(self, z: int) -> float:
        return max(elements.covalent_radius(z) * self.atom_scale, self.min_atom_radius)

    def atom_color(self, z: int) -> tuple[float, float, float]:
        return self.palette.get(int(z), elements.cpk_color(z))


def houk_style() -> Style:
    """The Houk-group look ("Houkmol" in CYLview): glossy ball-and-stick with
    black bonds, near-white carbons and two black great-circle seam lines
    ("quadrants") on every heavy atom, rotating with the molecule."""
    return Style(
        name="houk",
        bond_radius=0.12,
        atom_scale=0.47,
        min_atom_radius=0.24,
        palette={**CYLVIEW_COLORS,
                 6: (0.90, 0.90, 0.90),       # C near-white
                 7: (0.38, 0.40, 0.90)},      # N soft blue-violet
        bond_color=(0.07, 0.07, 0.07),
        quadrant_color=(0.05, 0.05, 0.05),
    )


STYLE_PRESETS = {"cylview": Style, "houk": houk_style}


def make_style(name: str) -> Style:
    return STYLE_PRESETS[name]()
