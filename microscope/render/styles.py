"""Rendering styles: the two presets, and reading or writing your own.

A style is plain data, so a group can keep its own as a JSON file that both
the viewer and `scope -s --style ours.json` use, and every figure in a paper
comes out matching.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, fields
from pathlib import Path

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
    """A preset by name, or a style file by path."""
    if name in STYLE_PRESETS:
        return STYLE_PRESETS[name]()
    if not Path(name).is_file():
        raise StyleError(f"no style preset or file called {name!r} "
                         f"(presets: {', '.join(sorted(STYLE_PRESETS))})")
    return load_style(name)

# ---------------------------------------------------------------- style files

COLOR_FIELDS = ("background", "hbond_color", "surface_positive",
                "surface_negative", "bond_color", "quadrant_color")


class StyleError(ValueError):
    """A style file that cannot be read as one."""


def _to_hex(rgb) -> str:
    return "#" + "".join(f"{round(max(0.0, min(1.0, c)) * 255):02x}" for c in rgb)


def _from_hex(value, where: str):
    """A colour written as #rrggbb, or as three numbers if you prefer."""
    if isinstance(value, (list, tuple)):
        if len(value) != 3:
            raise StyleError(f"{where}: a colour needs three numbers")
        return tuple(float(c) for c in value)
    text = str(value).strip().lstrip("#")
    if len(text) != 6:
        raise StyleError(f"{where}: {value!r} is not #rrggbb")
    try:
        return tuple(int(text[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    except ValueError:
        raise StyleError(f"{where}: {value!r} is not #rrggbb") from None


def style_to_dict(style: Style) -> dict:
    """A style as plain data, with colours as #rrggbb so it can be hand-edited."""
    out = {}
    for f in fields(style):
        value = getattr(style, f.name)
        if f.name == "palette":
            out["palette"] = {elements.SYMBOLS[z]: _to_hex(c)
                              for z, c in sorted(value.items())}
        elif f.name in COLOR_FIELDS:
            out[f.name] = None if value is None else _to_hex(value)
        else:
            out[f.name] = value
    return out


def style_from_dict(data: dict) -> Style:
    """A style from plain data, saying which key is wrong rather than guessing."""
    known = {f.name for f in fields(Style)}
    unknown = set(data) - known
    if unknown:
        raise StyleError("not part of a style: " + ", ".join(sorted(unknown)))

    kwargs = {}
    for key, value in data.items():
        if key == "palette":
            palette = {}
            for symbol, colour in (value or {}).items():
                z = elements.symbol_to_z(symbol)
                if not z:
                    raise StyleError(f"palette: {symbol!r} is not an element")
                palette[z] = _from_hex(colour, f"palette[{symbol}]")
            kwargs["palette"] = palette
        elif key in COLOR_FIELDS:
            kwargs[key] = None if value is None else _from_hex(value, key)
        else:
            kwargs[key] = value
    return Style(**kwargs)


def save_style(path, style: Style) -> None:
    """Write a style as JSON, ready to share or edit by hand."""
    Path(path).write_text(json.dumps(style_to_dict(style), indent=2) + "\n")


def load_style(path) -> Style:
    """Read a style file written by :func:`save_style`, or edited by hand."""
    try:
        data = json.loads(Path(path).read_text())
    except OSError as exc:
        raise StyleError(f"{path}: {exc.strerror}") from None
    except json.JSONDecodeError as exc:
        raise StyleError(f"{path}: not valid JSON ({exc.msg}, line {exc.lineno})") from None
    if not isinstance(data, dict):
        raise StyleError(f"{path}: a style file holds one object")
    try:
        return style_from_dict(data)
    except StyleError as exc:
        raise StyleError(f"{path}: {exc}") from None
