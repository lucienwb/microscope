"""Turning command-line text into values, and setting up the scene.

Kept apart from anything that draws so the tests can exercise every flag
without a display, and so `scope` can reject a bad argument before it has
paid for graphics.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from .. import io as mio
from ..core.geometry import rotation_matrix
from ..core.molecule import Molecule
from ..render.camera import OrthoCamera, orientation_along, orientation_from_plane
from ..render.formats import IMAGE_SUFFIXES, PLOT_SUFFIXES, VECTOR_SUFFIXES
from ..render.scene import REP_BALL, REP_LINE, REP_STICK
from ..render.styles import StyleError, make_style

REP_CODES = {"ball": REP_BALL, "stick": REP_STICK, "line": REP_LINE}
LABEL_MODES = ("none", "element", "element+number", "number")

# a Lewis drawing is line art, and can also leave as a structure file
LEWIS_SUFFIXES = IMAGE_SUFFIXES + VECTOR_SUFFIXES + mio.LEWIS_EXTENSIONS
DEFAULT_FWHM = {"ir": 8.0, "uv": 0.30, "nmr": 0.5}   # cm-1, eV, ppm

# Named orientations as world->view rotations: rows are the world directions
# of the view x (right), y (up) and z (towards the viewer) axes.
NAMED_VIEWS = {
    "front":  np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
    "back":   np.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]]),
    "top":    np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]]),
    "bottom": np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]]),
    "right":  np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]]),
    "left":   np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]]),
}


class CliError(Exception):
    """A bad command line — reported without a traceback."""


def parse_size(text: str) -> tuple[int, int]:
    """'1200x900' -> (1200, 900)."""
    parts = str(text).lower().replace("×", "x").split("x")
    if len(parts) != 2:
        raise CliError(f"--size wants WIDTHxHEIGHT, got {text!r}")
    try:
        w, h = (int(p) for p in parts)
    except ValueError:
        raise CliError(f"--size wants whole numbers, got {text!r}") from None
    if w < 16 or h < 16:
        raise CliError("--size must be at least 16x16")
    return w, h


def parse_indices(text: str, natoms: int | None = None) -> list[int]:
    """'1-6,9' -> [0, 1, 2, 3, 4, 5, 8] (1-based in, 0-based out)."""
    out: list[int] = []
    for chunk in str(text).split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        try:
            if "-" in chunk:
                lo, hi = (int(v) for v in chunk.split("-", 1))
                if hi < lo:
                    raise CliError(f"{chunk!r} counts backwards")
                out.extend(range(lo, hi + 1))
            else:
                out.append(int(chunk))
        except ValueError:
            raise CliError(f"cannot read atom numbers from {text!r}") from None
    if not out:
        raise CliError(f"no atom numbers in {text!r}")
    if any(i < 1 for i in out):
        raise CliError(f"atom numbers are 1-based, got {text!r}")
    if natoms is not None and any(i > natoms for i in out):
        raise CliError(f"atom number out of range in {text!r} "
                       f"(the molecule has {natoms} atoms)")
    return [i - 1 for i in out]


def parse_color(text: str) -> tuple[float, float, float]:
    """Any Qt colour name or #rrggbb -> an (r, g, b) triple of 0-1 floats."""
    from PySide6.QtGui import QColor
    color = QColor(str(text))
    if not color.isValid():
        raise CliError(f"{text!r} is not a colour name or #rrggbb value")
    return (color.redF(), color.greenF(), color.blueF())


def parse_reps(spec: str, natoms: int) -> np.ndarray:
    """'ball' or '1-6,9:ball;10-20:line' -> one REP_* code per atom."""
    reps = np.zeros(natoms, dtype=int)
    for group in str(spec).split(";"):
        group = group.strip()
        if not group:
            continue
        if ":" in group:
            ranges, name = group.rsplit(":", 1)
            atoms = parse_indices(ranges, natoms)
        else:
            name, atoms = group, list(range(natoms))
        name = name.strip().lower()
        if name not in REP_CODES:
            raise CliError(f"unknown representation {name!r} "
                           f"(pick from {', '.join(REP_CODES)})")
        reps[atoms] = REP_CODES[name]
    return reps


def parse_rotation(text: str) -> tuple[float, float]:
    """'30,-15' -> (30.0, -15.0) degrees."""
    parts = str(text).split(",")
    if len(parts) != 2:
        raise CliError(f"expected two angles like 30,-15 — got {text!r}")
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        raise CliError(f"expected two angles like 30,-15 — got {text!r}") from None


def default_output(input_path: str | None, kind: str | None = None) -> str:
    stem = Path(input_path).stem if input_path else "molecule"
    return f"{stem}_{kind}.png" if kind else f"{stem}.png"


def output_suffixes(kind: str | None) -> tuple[str, ...]:
    if kind == "lewis":
        return LEWIS_SUFFIXES
    return PLOT_SUFFIXES if kind else IMAGE_SUFFIXES


def check_output(path: str, kind: str | None = None) -> str:
    suffix = Path(path).suffix.lower()
    allowed = output_suffixes(kind)
    if suffix not in allowed:
        what = suffix or "a file with no extension"
        raise CliError(
            f"cannot write {what}; use " + (
                ", ".join(allowed) if kind else
                ".png or .tif (only formats that keep transparency), "
                "or .cube for the grid of an orbital"))
    return path


def parse_range(text: str) -> tuple[float, float]:
    """'400:1800' -> (400.0, 1800.0)."""
    parts = str(text).replace("..", ":").split(":")
    if len(parts) != 2:
        raise CliError(f"--xrange wants MIN:MAX, got {text!r}")
    try:
        lo, hi = (float(p) for p in parts)
    except ValueError:
        raise CliError(f"--xrange wants numbers, got {text!r}") from None
    if lo == hi:
        raise CliError("--xrange is empty")
    return lo, hi


def resolve_size(args, kind: str | None = None) -> tuple[int, int]:
    """Pixel size: spectra want a wider frame than a molecule does."""
    wide = kind is not None and kind != "lewis"
    return parse_size(args.size or ("1600x900" if wide else "1200x900"))


def spectrum_kind(args) -> str | None:
    """Which spectrum was asked for, if any."""
    for kind in ("ir", "uv", "nmr"):
        if getattr(args, kind):
            return kind
    return None


# ------------------------------------------------------------------ scene setup

def pick_frame(frames: list[Molecule], number: int | None) -> Molecule:
    """1-based frame selection; negatives count back, None means the last."""
    if not frames:
        raise CliError("no structure found in the file")
    if number is None:
        return frames[-1]
    if number == 0 or abs(number) > len(frames):
        raise CliError(f"--frame {number} is out of range "
                       f"(the file has {len(frames)} frames)")
    return frames[number - 1 if number > 0 else number]


def build_camera(molecule: Molecule, view: str | None, align: str | None,
                 rotate: str | None, zoom: float,
                 spin: float = 0.0, fit=None) -> OrthoCamera:
    """The camera the flags describe. *fit*, given the turned camera, frames
    it before --zoom is applied, so a zoom is relative to that framing."""
    camera = OrthoCamera()
    camera.fit(*molecule.bounding_sphere())
    if align:
        idxs = parse_indices(align, molecule.natoms)
        coords = molecule.coords
        if len(idxs) == 2:
            rot = orientation_along(coords[idxs[1]] - coords[idxs[0]])
        elif len(idxs) == 3:
            rot = orientation_from_plane(*(coords[i] for i in idxs))
        else:
            raise CliError("--align wants 2 atoms (look down the bond) "
                           "or 3 atoms (put their plane in the screen)")
        if rot is None:
            raise CliError("--align atoms are collinear or coincident")
        camera.rotation = rot
    elif view:
        if view.lower() in NAMED_VIEWS:
            camera.rotation = NAMED_VIEWS[view.lower()].copy()
        else:
            rx, ry = parse_rotation(view)
            _spin(camera, rx, ry)
    if rotate:
        _spin(camera, *parse_rotation(rotate))
    if spin:
        camera.spin(spin)
    if fit is not None:
        fit(camera)
    if zoom and zoom > 0:
        camera.half_height /= zoom
    elif zoom:
        raise CliError("--zoom must be positive")
    return camera


FIGURE_MARGIN = 0.05          # of the frame, clear on each side of the molecule


def figure_extent(molecule: Molecule, style, reps=None, volume=None,
                  isovalue: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """What a figure shows, as spheres (points, radii): every atom at the size
    it is drawn, and the lobes of an isosurface, which reach well past the
    atoms."""
    from ..render.scene import LINE_RADIUS, REP_LINE, REP_STICK

    zs = molecule.atomic_numbers
    elements, which = np.unique(zs, return_inverse=True)
    radii = np.array([style.atom_radius(z) for z in elements])[which]
    if reps is not None and len(reps) == len(radii):
        reps = np.asarray(reps)
        radii[reps == REP_STICK] = style.bond_radius
        radii[reps == REP_LINE] = LINE_RADIUS
    points, sizes = [molecule.coords], [radii]
    if volume is not None:
        level = abs(isovalue) if isovalue else volume.suggest_isovalue()
        inside = np.argwhere(np.abs(volume.values) >= level)
        if len(inside):
            inside = inside[::max(1, len(inside) // 200_000)]    # the extent, not every point
            points.append(volume.grid_to_world(inside))
            sizes.append(np.full(len(inside), float(np.abs(volume.axes).max())))
    return np.vstack(points), np.concatenate(sizes)


def frame_figure(camera: OrthoCamera, molecule: Molecule, style, aspect: float,
                 reps=None, volume=None, isovalue: float | None = None) -> None:
    """Frame a figure on what it shows from this angle, so it needs no cropping."""
    points, radii = figure_extent(molecule, style, reps, volume, isovalue)
    camera.frame(points, radii, aspect, margin=FIGURE_MARGIN)


def _spin(camera: OrthoCamera, rx: float, ry: float) -> None:
    """Turn the camera by *rx* about the screen's vertical axis and *ry* about
    its horizontal one — the same motion as dragging in the viewer."""
    camera.rotation = rotation_matrix(np.array([0.0, 1.0, 0.0]),
                                      np.radians(rx)) @ camera.rotation
    camera.rotation = rotation_matrix(np.array([1.0, 0.0, 0.0]),
                                      np.radians(ry)) @ camera.rotation


def build_style(args):
    try:
        style = make_style(args.style)
    except StyleError as exc:
        raise CliError(str(exc)) from None
    except OSError:
        raise CliError(f"--style: no preset or file called {args.style!r}") from None
    if args.background:
        style.background = parse_color(args.background)
    if args.bond_color:
        style.bond_color = parse_color(args.bond_color)
    if args.no_hbonds:
        style.show_hbonds = False
    if args.iso_colors:
        parts = args.iso_colors.split(",")
        if len(parts) != 2:
            raise CliError("--iso-colors wants two colours, e.g. blue,red")
        style.surface_positive = parse_color(parts[0])
        style.surface_negative = parse_color(parts[1])
    if args.iso_opacity is not None:
        style.surface_opacity = float(np.clip(args.iso_opacity, 0.0, 1.0))
    return style


def load_volume(args, result):
    """The grid to draw an isosurface from: --cube wins, then an orbital from
    the file's wavefunction, then a grid the file holds itself."""
    if args.cube:
        volumes = list(mio.load(args.cube).volumes)
        if not volumes:
            raise CliError(f"no volumetric data in {args.cube}")
        return _pick_grid(args.mo, volumes)
    if args.mo is not None and result.orbitals is not None:
        if result.orbitals.convention == "unrecognized":
            print("scope: warning: this file's basis-set conventions were not "
                  "recognized, so the orbital may be drawn wrong", file=sys.stderr)
        try:
            return result.orbitals.volume(args.mo)
        except ValueError as exc:
            raise CliError(f"--mo: {exc}") from None
    if not result.volumes:
        if args.mo is not None and result.warnings:
            raise CliError(f"--mo: {Path(args.file).name} has a wavefunction, but "
                           + "; ".join(result.warnings))
        if args.mo is not None:
            raise CliError(f"--mo needs a wavefunction (.fchk or .molden) or a cube "
                           f"file, and {Path(args.file).name} has neither")
        return None
    return _pick_grid(args.mo, result.volumes)


def _pick_grid(spec, volumes):
    """--mo on a cube file: the grid's position, since a cube has no occupations."""
    if spec is None:
        return volumes[0]
    try:
        number = int(spec)
    except ValueError:
        held = "one grid" if len(volumes) == 1 else f"grids 1-{len(volumes)}"
        raise CliError(f"--mo {spec}: a cube file says nothing about occupations, so "
                       f"give the grid's number (it holds {held})") from None
    if not 1 <= number <= len(volumes):
        raise CliError(f"--mo {number} is out of range "
                       f"(the cube holds {len(volumes)} grids)")
    return volumes[number - 1]
