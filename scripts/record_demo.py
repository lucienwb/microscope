"""Record the animated GIFs the README uses.

    python scripts/record_demo.py            # all of them
    python scripts/record_demo.py orbit      # just one

Everything is rendered offscreen through the same code the viewer and
`scope -s` use, so a clip shows what the program actually draws. No window
opens. The 3-D clips need a working OpenGL 3.3 driver; the Lewis ones are
painted with QPainter and record anywhere.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from PySide6.QtCore import QBuffer
from PySide6.QtGui import QGuiApplication, QPainter

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import microscope  # noqa: E402
from microscope.core import editing  # noqa: E402
from microscope.core import lewis as lewis_core  # noqa: E402
from microscope.core import measure as measure_model  # noqa: E402
from microscope.gui import annotations  # noqa: E402
from microscope.render import lewis2d, lewisdraw  # noqa: E402
from microscope.render.camera import OrthoCamera  # noqa: E402
from microscope.render.offscreen import render_molecule_image  # noqa: E402
from microscope.render.scene import REP_BALL, REP_LINE, REP_STICK  # noqa: E402
from microscope.render.styles import make_style  # noqa: E402

OUT = REPO / "docs"
DATA = REPO / "tests/data"
W, H = 560, 400
COLORS = 64          # enough for shaded spheres, small enough for a README


# ------------------------------------------------------------------ machinery

def _to_pil(qimage) -> Image.Image:
    buf = QBuffer()
    buf.open(QBuffer.OpenModeFlag.ReadWrite)
    qimage.save(buf, "PNG")
    return Image.open(io.BytesIO(buf.data().data())).convert("RGB")


def _shot(molecule, camera, *, style="cylview", reps=None, volume=None,
          iso=None, overlay=None) -> Image.Image:
    """One frame: the molecule as rendered, plus whatever the viewer draws on top."""
    image = render_molecule_image(molecule, make_style(style), camera, W, H,
                                  supersample=2, transparent=False,
                                  volume=volume, isovalue=iso, reps=reps)
    if overlay:
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        annotations.draw_overlay(painter, molecule, camera, W, H, **overlay)
        painter.end()
    return _to_pil(image)


def _save(frames, name: str, ms: int = 90) -> None:
    if len(frames[0].getcolors(maxcolors=8) or [0] * 9) < 3:
        raise SystemExit(f"{name}: the frames came out blank — this needs a "
                         "working OpenGL 3.3 driver.")
    palette = [f.convert("P", palette=Image.ADAPTIVE, colors=COLORS) for f in frames]
    path = OUT / name
    palette[0].save(path, save_all=True, append_images=palette[1:],
                    duration=ms, loop=0, optimize=True)
    print(f"  {path.relative_to(REPO)}  {len(frames)} frames  "
          f"{path.stat().st_size / 1024:.0f} KB")


def _load(name: str):
    result = microscope.load(DATA / name)
    molecule = result.molecule
    if molecule.bonds is None:
        molecule.perceive_bonds()
    return result, molecule


def _camera(molecule, zoom: float = 1.0) -> OrthoCamera:
    camera = OrthoCamera()
    camera.fit(*molecule.bounding_sphere())
    camera.half_height /= zoom
    return camera


def _hold(frames, times: int = 6):
    """Linger on the last frame so a loop reads as a step, not a flicker."""
    frames.extend([frames[-1]] * times)


# ---------------------------------------------------------------------- clips

def orbit():
    """Left-drag: turn the molecule."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    frames = []
    for _ in range(36):
        cam.rotate_drag(18, 0)
        frames.append(_shot(mol, cam))
    _save(frames, "orbit.gif", ms=70)


def styles():
    """V: the CYLview look and the Houk (Houkmol) ball-and-stick style."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-40, -30)
    frames = []
    for name in ("cylview", "houk"):
        frames.append(_shot(mol, cam, style=name))
        _hold(frames, 10)
    _save(frames, "styles.gif", ms=110)


def regions():
    """1 / 2 / 3: give the selected region its own level of detail."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-40, -30)
    ring = [a for a in range(mol.natoms)
            if mol.symbols[a] in ("N", "O")]          # the polar head group
    frames = []
    for code in (REP_BALL, REP_STICK, REP_LINE):
        reps = np.full(mol.natoms, code, dtype=int)
        reps[ring] = REP_BALL
        frames.append(_shot(mol, cam, reps=reps))
        _hold(frames, 9)
    _save(frames, "regions.gif", ms=110)


def labels():
    """L: cycle element, element+number, number, off."""
    _, mol = _load("water_mp2.log")
    cam = _camera(mol, zoom=0.75)
    frames = []
    for mode in ("element", "element+number", "number", "none"):
        frames.append(_shot(mol, cam, overlay={"label_mode": mode}))
        _hold(frames, 8)
    _save(frames, "labels.gif", ms=110)


def measure():
    """Click atoms: distance, then angle, then dihedral."""
    _, mol = _load("dvb_gopt.out")
    cam = _camera(mol)
    cam.rotate_drag(-20, -25)
    frames = []
    for picked in ([0, 1], [0, 1, 2], [0, 1, 2, 3]):
        overlay = {"selection": picked,
                   "selection_text": measure_model.describe(mol, picked)}
        frames.append(_shot(mol, cam, overlay=overlay))
        _hold(frames, 10)
    _save(frames, "measure.gif", ms=110)


def handles():
    """F then the manipulator: slide a fragment along an axis and turn it."""
    result, mol = _load("water_dimer.xyz")
    del result
    cam = _camera(mol, zoom=1.05)
    cam.rotate_drag(-20, -10)
    fragment = editing.connected_fragment(mol.bonds, [3])
    base = mol.coords.copy()
    frames = []
    for step in list(np.linspace(0, 0.9, 12)) + list(np.linspace(0.9, 0, 8)):
        mol.coords = base.copy()
        mol.coords[fragment] += np.array([step, 0.0, 0.0])
        centre = mol.coords[fragment].mean(axis=0)
        frames.append(_shot(mol, cam, overlay={
            "selection": fragment, "gizmo_center": centre}))
    for angle in np.linspace(0, 50, 10):
        mol.coords = base.copy()
        mol.coords = editing.rotate_atoms(
            mol, fragment, np.array([0.0, 0.0, 1.0]), float(angle),
            pivot=base[fragment].mean(axis=0))
        centre = mol.coords[fragment].mean(axis=0)
        frames.append(_shot(mol, cam, overlay={
            "selection": fragment, "gizmo_center": centre}))
    mol.coords = base
    _save(frames, "handles.gif", ms=80)


def axes():
    """Shift+A: the corner XYZ triad, which follows the rotation."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    frames = []
    for _ in range(24):
        cam.rotate_drag(0, 15)
        frames.append(_shot(mol, cam, overlay={"show_axes": True}))
    _save(frames, "axes.gif", ms=80)


def isosurface():
    """I: sweep the isovalue of an orbital from a cube file."""
    result, mol = _load("water_mo.cube")
    volume = result.volumes[0]
    cam = _camera(mol, zoom=0.7)
    cam.rotate_drag(-25, -15)
    levels = list(np.linspace(0.02, 0.11, 10))
    frames = [_shot(mol, cam, volume=volume, iso=float(v)) for v in levels]
    frames += frames[::-1]
    _save(frames, "isosurface.gif", ms=90)


def vibrate():
    """Click an IR band: the molecule walks through that normal mode."""
    result, mol = _load("dvb_ir.out")
    mode = max(result.vibrations, key=lambda v: v.ir_intensity)
    cam = _camera(mol)
    cam.rotate_drag(-30, -20)
    base = mol.coords.copy()
    disp = np.asarray(mode.displacements, dtype=float).reshape(-1, 3)
    frames = []
    for phase in np.linspace(0, 2 * np.pi, 20, endpoint=False):
        mol.coords = base + 0.6 * np.sin(phase) * disp
        frames.append(_shot(mol, cam))
    mol.coords = base
    _save(frames, "vibration.gif", ms=60)


def trajectory():
    """An optimization played back, geometry by geometry."""
    result, _ = _load("g09_dvb_scan.log")
    frames_in = result.frames
    cam = _camera(frames_in[-1], zoom=0.85)
    cam.rotate_drag(-30, -20)
    step = max(1, len(frames_in) // 24)
    frames = []
    for mol in frames_in[::step]:
        if mol.bonds is None:
            mol.perceive_bonds()
        frames.append(_shot(mol, cam))
    _hold(frames, 8)
    _save(frames, "trajectory.gif", ms=110)


def lewis():
    """Shift+L, then turning the drawing to pick the angle it is drawn from."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-80, -80)
    structure = lewis_core.perceive(mol)
    frames = []
    for _ in range(36):
        cam.spin(10)
        frames.append(_to_pil(lewisdraw.render_lewis_image(
            structure, cam, W, H, lewis2d.LewisOptions(), transparent=False)))
    _save(frames, "lewis_spin.gif", ms=80)


def lewis_options():
    """What the drawing can show: hydrogens, carbons, lone pairs, colour."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-80, -80)
    cam.spin(180)
    structure = lewis_core.perceive(mol)
    frames = []
    for options in (lewis2d.LewisOptions(),
                    lewis2d.LewisOptions(color_atoms=True),
                    lewis2d.LewisOptions(color_atoms=True, lone_pairs=True),
                    lewis2d.LewisOptions(color_atoms=True, lone_pairs=True,
                                         show_hydrogens=True)):
        frames.append(_to_pil(lewisdraw.render_lewis_image(
            structure, cam, W, H, options, transparent=False)))
        _hold(frames, 9)
    _save(frames, "lewis_options.gif", ms=110)


CLIPS = {f.__name__: f for f in (orbit, styles, regions, labels, measure,
                                 handles, axes, isosurface, vibrate,
                                 trajectory, lewis, lewis_options)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clips", nargs="*", choices=[*CLIPS, []],
                        help="which to record (default: all)")
    args = parser.parse_args()

    QGuiApplication(sys.argv[:1])       # needed for the GL context and QImage
    for name in args.clips or CLIPS:
        print(f"{name}: {CLIPS[name].__doc__.splitlines()[0]}")
        CLIPS[name]()
    return 0


if __name__ == "__main__":
    sys.exit(main())
