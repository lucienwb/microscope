"""Offscreen render smoke test.

Usage: python scripts/preview.py [input_file] [output.png] [--style cylview|houk]
Renders the file (or a built-in benzene) to a PNG without opening a window.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtGui import QGuiApplication  # noqa: E402

import microscope.io as mio  # noqa: E402
from microscope.core.molecule import Molecule  # noqa: E402
from microscope.render.camera import OrthoCamera  # noqa: E402
from microscope.render.offscreen import render_molecule_image  # noqa: E402
from microscope.render.styles import make_style  # noqa: E402


def benzene() -> Molecule:
    rc, rh = 1.397, 2.481
    symbols, coords = [], []
    for k in range(6):
        a = math.radians(60 * k)
        symbols.append("C")
        coords.append([rc * math.cos(a), rc * math.sin(a), 0.0])
        symbols.append("H")
        coords.append([rh * math.cos(a), rh * math.sin(a), 0.0])
    return Molecule(symbols, np.array(coords), title="benzene")


def main() -> None:
    QGuiApplication(sys.argv)
    args = [a for a in sys.argv[1:]]
    style_name = "cylview"
    if "--style" in args:
        k = args.index("--style")
        style_name = args[k + 1]
        del args[k:k + 2]
    out = "preview.png"
    mol = None
    volume = None
    if args and not args[0].lower().endswith(".png"):
        result = mio.load(args[0])
        mol = result.molecule
        volume = result.volumes[0] if result.volumes else None
        args = args[1:]
    if args:
        out = args[0]
    if mol is None:
        mol = benzene()
    mol.perceive_bonds()

    camera = OrthoCamera()
    center, radius = mol.bounding_sphere()
    camera.fit(center, radius)
    camera.rotate_drag(55, 35)   # tilt so depth is visible

    image = render_molecule_image(mol, make_style(style_name), camera, 900, 700,
                                  supersample=3, volume=volume)
    image.save(out)
    extra = f", isosurface {volume.label!r}" if volume is not None else ""
    print(f"wrote {out} ({mol.formula()}, {mol.natoms} atoms{extra})")


if __name__ == "__main__":
    main()
