"""Record the animations and pictures the documentation uses.

    python scripts/record_demo.py            # all of them
    python scripts/record_demo.py orbit      # just one

Everything is rendered offscreen through the same code the viewer and
`scope -s` use, so a clip shows what the program actually draws. No window
opens. The 3-D clips need a working OpenGL 3.3 driver; the Lewis ones are
painted with QPainter and record anywhere.

The 3-D clips are animated WebP on a transparent background, so they sit on
the site's light and dark themes alike, in full colour rather than a GIF's 64.
Each is recorded twice: a first pass only notes what every frame will show,
one framing is fitted to all of it - the molecule fills the clip however it
turns - and the second pass draws. Line art (the Lewis clips) stays on white,
like a page.
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
from microscope.cli.parsing import figure_extent  # noqa: E402
from microscope.core import editing  # noqa: E402
from microscope.core import lewis as lewis_core  # noqa: E402
from microscope.core import measure as measure_model  # noqa: E402
from microscope.gui import annotations  # noqa: E402
from microscope.render import lewis2d, lewisdraw  # noqa: E402
from microscope.render.camera import OrthoCamera  # noqa: E402
from microscope.render.offscreen import render_molecule_image  # noqa: E402
from microscope.render.scene import REP_BALL, REP_LINE, REP_STICK  # noqa: E402
from microscope.render.styles import make_style  # noqa: E402

OUT = REPO / "docs" / "assets"
DATA = REPO / "tests/data"
W, H = 560, 400
MARGIN = 0.05        # of the frame, clear round the molecule
_PLAN: list | None = None                     # first pass: what each frame shows
_FIT: tuple[np.ndarray, float] | None = None  # second pass: the framing all frames use


# ------------------------------------------------------------------ machinery

def _to_pil(qimage, alpha: bool = True) -> Image.Image:
    buf = QBuffer()
    buf.open(QBuffer.OpenModeFlag.ReadWrite)
    qimage.save(buf, "PNG")
    return Image.open(io.BytesIO(buf.data().data())).convert("RGBA" if alpha else "RGB")


def _shot(molecule, camera, *, style="cylview", reps=None, volume=None,
          iso=None, overlay=None, caption=None, size=(W, H)) -> Image.Image:
    """One frame: the molecule as rendered, plus whatever the viewer draws on top."""
    style = make_style(style) if isinstance(style, str) else style
    if _PLAN is not None:            # first pass: note it, draw nothing
        points, radii = figure_extent(molecule, style, reps, volume, iso)
        _PLAN.append((camera.rotation.copy(), camera.center.copy(), points, radii))
        return Image.new("RGBA", size)
    if _FIT is not None:
        camera.center, camera.half_height = _FIT[0].copy(), _FIT[1]
    width, height = size
    image = render_molecule_image(molecule, style, camera, width, height,
                                  supersample=2, transparent=True,
                                  volume=volume, isovalue=iso, reps=reps)
    if overlay or caption:
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if overlay:
            annotations.draw_overlay(painter, molecule, camera, width, height, **overlay)
        if caption:           # what the viewer's orbital list says about it
            font = painter.font()
            font.setPixelSize(15)
            painter.setFont(font)
            annotations.draw_halo_text(painter, 16, height - 16, caption)
        painter.end()
    return _to_pil(image)


def _fitting(plan: list, aspect: float, margin: float) -> tuple[np.ndarray, float]:
    """One framing that holds every frame of a clip. A clip seen from one angle
    is centred on what it shows; a turning one keeps its centre of rotation and
    is sized for the widest it gets."""
    rotation, centre = plan[0][0], plan[0][1]
    lows, highs = [], []
    for turn, _, points, radii in plan:
        seen = (points - centre) @ turn[:2].T
        lows.append((seen - radii[:, None]).min(axis=0))
        highs.append((seen + radii[:, None]).max(axis=0))
    low, high = np.min(lows, axis=0), np.max(highs, axis=0)
    if all(np.allclose(turn, rotation) for turn, *_ in plan):
        middle = (low + high) / 2
        centre = centre + rotation[0] * middle[0] + rotation[1] * middle[1]
        half = (high - low) / 2
    else:
        half = np.maximum(-low, high)
    return centre, max(half[1], half[0] / aspect, 0.5) / (1.0 - 2.0 * margin)


def _record(clip, margin: float = MARGIN, aspect: float = W / H) -> None:
    """Run *clip* twice: once to see what it shows, once to draw it framed."""
    global _PLAN, _FIT
    _PLAN, _FIT = [], None
    try:
        clip()
        plan = _PLAN
    finally:
        _PLAN = None
    _FIT = _fitting(plan, aspect, margin) if plan else None
    try:
        clip()
    finally:
        _FIT = None


def _save(frames, name: str, ms: int = 90) -> None:
    if _PLAN is not None:            # the first pass writes nothing
        return
    probe = np.asarray(frames[0])
    visible = probe[..., 3] > 0 if probe.shape[-1] == 4 else probe.min(axis=-1) < 250
    if not visible.any():
        raise SystemExit(f"{name}: the frames came out blank — this needs a "
                         "working OpenGL 3.3 driver.")
    path = OUT / name
    frames[0].save(path, "WEBP", save_all=True, append_images=frames[1:], duration=ms,
                   loop=0, quality=88, method=6)
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


def _landscape(camera: OrthoCamera, coords: np.ndarray) -> None:
    """Turn the picture in the plane of the screen until the molecule's long
    axis runs across it: the clips are wider than they are tall."""
    seen = (coords - coords.mean(axis=0)) @ camera.rotation[:2].T
    values, vectors = np.linalg.eigh(np.cov(seen.T))
    x, y = vectors[:, np.argmax(values)]
    camera.spin(-np.degrees(np.arctan2(y, x)))


def _face_on(camera: OrthoCamera, coords: np.ndarray) -> None:
    """Look down the molecule's thinnest direction, long axis across the
    screen: the flattest, widest view it has."""
    values, vectors = np.linalg.eigh(np.cov((coords - coords.mean(axis=0)).T))
    across, up = vectors[:, 2], vectors[:, 1]
    camera.rotation = np.array([across, up, np.cross(across, up)])


def _hold(frames, times: int = 6):
    """Linger on the last frame so a loop reads as a step, not a flicker."""
    frames.extend([frames[-1]] * times)


# ---------------------------------------------------------------------- clips

HERO = (1200, 520)


def hero():
    """The picture at the top of the README and the site."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    _face_on(cam, mol.coords)
    cam.rotate_drag(25, 35)
    image = _shot(mol, cam, size=HERO)
    if _PLAN is None:
        path = OUT / "hero.png"
        image.save(path)
        print(f"  {path.relative_to(REPO)}  {path.stat().st_size / 1024:.0f} KB")


def orbit():
    """Left-drag: turn the molecule."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    frames = []
    for _ in range(36):
        cam.rotate_drag(18, 0)
        frames.append(_shot(mol, cam))
    _save(frames, "orbit.webp", ms=70)


def styles():
    """V: the CYLview look and the Houk (Houkmol) ball-and-stick style."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-40, -30)
    _landscape(cam, mol.coords)
    frames = []
    for name in ("cylview", "houk"):
        frames.append(_shot(mol, cam, style=name))
        _hold(frames, 10)
    _save(frames, "styles.webp", ms=110)


def regions():
    """1 / 2 / 3: give the selected region its own level of detail."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-40, -30)
    _landscape(cam, mol.coords)
    ring = [a for a in range(mol.natoms)
            if mol.symbols[a] in ("N", "O")]          # the polar head group
    frames = []
    for code in (REP_BALL, REP_STICK, REP_LINE):
        reps = np.full(mol.natoms, code, dtype=int)
        reps[ring] = REP_BALL
        frames.append(_shot(mol, cam, reps=reps))
        _hold(frames, 9)
    _save(frames, "regions.webp", ms=110)


def labels():
    """L: cycle element, element+number, number, off."""
    _, mol = _load("water_mp2.log")
    cam = _camera(mol)
    _face_on(cam, mol.coords)          # the file has water edge-on, bent angle hidden
    frames = []
    for mode in ("element", "element+number", "number", "none"):
        frames.append(_shot(mol, cam, overlay={"label_mode": mode}))
        _hold(frames, 8)
    _save(frames, "labels.webp", ms=110)


def measure():
    """Click atoms: distance, then angle, then dihedral."""
    _, mol = _load("dvb_gopt.out")
    cam = _camera(mol)
    cam.rotate_drag(-20, -25)
    _landscape(cam, mol.coords)
    frames = []
    for picked in ([0, 1], [0, 1, 2], [0, 1, 2, 3]):
        overlay = {"selection": picked,
                   "selection_text": measure_model.describe(mol, picked)}
        frames.append(_shot(mol, cam, overlay=overlay))
        _hold(frames, 10)
    _save(frames, "measure.webp", ms=110)


def handles():
    """F then the manipulator: slide a fragment along an axis and turn it."""
    result, mol = _load("water_dimer.xyz")
    del result
    cam = _camera(mol)
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
    _save(frames, "handles.webp", ms=80)


def axes():
    """Shift+A: the corner XYZ triad, which follows the rotation."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    frames = []
    for _ in range(24):
        cam.rotate_drag(0, 15)
        frames.append(_shot(mol, cam, overlay={"show_axes": True}))
    _save(frames, "axes.webp", ms=80)


def isosurface():
    """I: sweep the isovalue of an orbital from a cube file."""
    result, mol = _load("water_mo.cube")
    volume = result.volumes[0]
    cam = _camera(mol)
    cam.rotate_drag(-25, -15)
    levels = list(np.linspace(0.02, 0.11, 10))
    frames = [_shot(mol, cam, volume=volume, iso=float(v)) for v in levels]
    frames += frames[::-1]
    _save(frames, "isosurface.webp", ms=90)


def orbitals():
    """I: step through the orbitals of a wavefunction file."""
    result, mol = _load("dvb_ir.fchk")
    cam = _camera(mol)
    cam.rotate_drag(40, -70)
    _landscape(cam, mol.coords)
    frames = []
    for which in ("homo-2", "homo-1", "homo", "lumo", "lumo+1"):
        volume = result.orbitals.volume(which)
        frames.append(_shot(mol, cam, volume=volume, iso=0.03, caption=volume.label))
        _hold(frames, 11)          # the viewer holds still too, as you step the list
    _save(frames, "orbitals.webp", ms=100)


def levels():
    """I: the orbital energy-level diagram beside the orbital list."""
    from microscope.render.levels import level_diagram
    from microscope.render.levelsdraw import write_levels

    result, _ = _load("dvb_ir.fchk")
    path = OUT / "levels.png"
    write_levels(path, level_diagram(result.orbitals), picked=("alpha", 34))
    print(f"  {path.relative_to(REPO)}  {path.stat().st_size / 1024:.0f} KB")


def vibrate():
    """Click an IR band: the molecule walks through that normal mode."""
    result, mol = _load("dvb_ir.out")
    mode = max(result.vibrations, key=lambda v: v.ir_intensity)
    cam = _camera(mol)
    cam.rotate_drag(-30, -20)
    _landscape(cam, mol.coords)
    base = mol.coords.copy()
    disp = np.asarray(mode.displacements, dtype=float).reshape(-1, 3)
    frames = []
    for phase in np.linspace(0, 2 * np.pi, 20, endpoint=False):
        mol.coords = base + 0.6 * np.sin(phase) * disp
        frames.append(_shot(mol, cam))
    mol.coords = base
    _save(frames, "vibration.webp", ms=60)


def trajectory():
    """An optimization played back, geometry by geometry."""
    result, _ = _load("g09_dvb_scan.log")
    frames_in = result.frames
    cam = _camera(frames_in[-1])
    cam.rotate_drag(-30, -20)
    _landscape(cam, frames_in[-1].coords)
    step = max(1, len(frames_in) // 24)
    frames = []
    for mol in frames_in[::step]:
        if mol.bonds is None:
            mol.perceive_bonds()
        frames.append(_shot(mol, cam))
    _hold(frames, 8)
    _save(frames, "trajectory.webp", ms=110)


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
            structure, cam, W, H, lewis2d.LewisOptions(), transparent=False), alpha=False))
    _save(frames, "lewis_spin.webp", ms=80)


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
            structure, cam, W, H, options, transparent=False), alpha=False))
        _hold(frames, 9)
    _save(frames, "lewis_options.webp", ms=110)


def style():
    """Ctrl+T: atom size, bond width and element colours, changing live."""
    _, mol = _load("trp.log")
    cam = _camera(mol)
    cam.rotate_drag(-40, -30)
    _landscape(cam, mol.coords)
    start, end = make_style("cylview"), make_style("cylview")
    end.atom_scale, end.bond_radius = 0.52, 0.10
    end.bond_color = (0.10, 0.10, 0.10)
    end.palette = {**end.palette, 6: (0.16, 0.42, 0.55), 7: (0.90, 0.62, 0.15),
                   8: (0.80, 0.25, 0.45)}

    def mix(t):
        step = make_style("cylview")
        step.atom_scale = start.atom_scale + t * (end.atom_scale - start.atom_scale)
        step.bond_radius = start.bond_radius + t * (end.bond_radius - start.bond_radius)
        step.bond_color = None if t < 0.5 else end.bond_color
        step.palette = {z: tuple(a + t * (b - a) for a, b in
                                 zip(start.palette.get(z, c), end.palette.get(z, c)))
                        for z, c in end.palette.items()}
        return step

    frames = []
    for t in list(np.linspace(0, 1, 14)) + [1.0] * 6 + list(np.linspace(1, 0, 8)):
        frames.append(_shot(mol, cam, style=mix(float(t))))
    _save(frames, "style.webp", ms=90)


CLIPS = {f.__name__: f for f in (hero, orbit, styles, style, regions, labels, measure,
                                 handles, axes, orbitals, levels, isosurface, vibrate,
                                 trajectory, lewis, lewis_options)}
# the 3-D clips and how much room each leaves round the molecule: the ones with
# the viewer's labels, measurements or handles on top need more of it
FRAMED = {"hero": MARGIN, "orbit": MARGIN, "styles": MARGIN, "style": MARGIN,
          "regions": MARGIN, "labels": 0.09, "measure": 0.10, "handles": 0.16,
          "axes": 0.10, "orbitals": 0.07, "isosurface": MARGIN, "vibrate": MARGIN,
          "trajectory": MARGIN}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clips", nargs="*", choices=[*CLIPS, []],
                        help="which to record (default: all)")
    args = parser.parse_args()

    QGuiApplication(sys.argv[:1])       # needed for the GL context and QImage
    for name in args.clips or CLIPS:
        print(f"{name}: {CLIPS[name].__doc__.splitlines()[0]}")
        if name in FRAMED:
            aspect = HERO[0] / HERO[1] if name == "hero" else W / H
            _record(CLIPS[name], FRAMED[name], aspect)
        else:
            CLIPS[name]()
    return 0


if __name__ == "__main__":
    sys.exit(main())
