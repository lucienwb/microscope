"""Offscreen render smoke test.

Usage: python scripts/preview.py [input_file] [output.png]
Renders the file (or a built-in benzene) to a PNG without opening a window.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtGui import QGuiApplication  # noqa: E402

import microscp.io as mio  # noqa: E402
from microscp.core.molecule import Molecule  # noqa: E402
from microscp.render.camera import OrthoCamera  # noqa: E402
from microscp.render.offscreen import render_molecule_image  # noqa: E402
from microscp.render.styles import Style  # noqa: E402


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
    out = "preview.png"
    mol = None
    if args and not args[0].lower().endswith(".png"):
        mol = mio.load(args[0]).molecule
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

    image = render_molecule_image(mol, Style(), camera, 900, 700, supersample=3)
    image.save(out)
    print(f"wrote {out} ({mol.formula()}, {mol.natoms} atoms)")


if __name__ == "__main__":
    main()
