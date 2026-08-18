"""The `scope` command line.

    scope mycalc.log                 open the viewer window
    scope -s mycalc.log              render mycalc.png and exit, no window
    scope -s mycalc.log --ir         plot the IR spectrum instead

Silent mode is the gnuplot-style half of the program: one command line turns a
quantum-chemistry output into a finished figure, so figures can be regenerated
from a shell script or a Makefile after every re-optimization. Everything the
GUI can put in an exported image — style, view, labels, representations,
measurements, cube isosurfaces, transparency — has a flag here, and
--ir/--uv/--nmr switch from drawing the molecule to plotting its spectrum.

Atom numbers on the command line are 1-based, matching the labels shown in the
viewer (C1, C2, …).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

from . import __version__
from . import io as mio
from .core.molecule import Molecule
from .render.camera import (OrthoCamera, orientation_along,
                            orientation_from_plane)
from .core.geometry import rotation_matrix
from .io.errors import FileFormatError, UnsupportedFormatError
from .render.imagefile import IMAGE_SUFFIXES, write_image
from .render.scene import REP_BALL, REP_LINE, REP_STICK
from .render.styles import make_style

REP_CODES = {"ball": REP_BALL, "stick": REP_STICK, "line": REP_LINE}
LABEL_MODES = ("none", "element", "element+number", "number")

# matplotlib writes these; vector formats are what journals usually want
PLOT_SUFFIXES = (".png", ".pdf", ".svg", ".eps", ".tif", ".tiff")
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


# ------------------------------------------------------------------ parsing

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


def check_output(path: str, kind: str | None = None) -> str:
    suffix = Path(path).suffix.lower()
    allowed = PLOT_SUFFIXES if kind else IMAGE_SUFFIXES
    if suffix not in allowed:
        what = suffix or "a file with no extension"
        raise CliError(
            f"cannot write {what}; use " + (
                ", ".join(allowed) if kind else
                ".png or .tif (only formats that keep transparency)"))
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
    return parse_size(args.size or ("1600x900" if kind else "1200x900"))


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
                 rotate: str | None, zoom: float) -> OrthoCamera:
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
    if zoom and zoom > 0:
        camera.half_height /= zoom
    elif zoom:
        raise CliError("--zoom must be positive")
    return camera


def _spin(camera: OrthoCamera, rx: float, ry: float) -> None:
    """Turn the camera by *rx* about the screen's vertical axis and *ry* about
    its horizontal one — the same motion as dragging in the viewer."""
    camera.rotation = rotation_matrix(np.array([0.0, 1.0, 0.0]),
                                      np.radians(rx)) @ camera.rotation
    camera.rotation = rotation_matrix(np.array([1.0, 0.0, 0.0]),
                                      np.radians(ry)) @ camera.rotation


def build_style(args):
    style = make_style(args.style)
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
    """The grid to draw an isosurface from: --cube wins, else one in the file."""
    volumes = list(result.volumes)
    if args.cube:
        volumes = list(mio.load(args.cube).volumes)
        if not volumes:
            raise CliError(f"no volumetric data in {args.cube}")
    if not volumes:
        return None
    if args.mo is None:
        return volumes[0]
    if not 1 <= args.mo <= len(volumes):
        raise CliError(f"--mo {args.mo} is out of range "
                       f"(the cube holds {len(volumes)} grids)")
    return volumes[args.mo - 1]


# ------------------------------------------------------------------ spectra

def pick_nucleus(shieldings, wanted: str | None) -> str:
    """The element to plot: what was asked for, else H, else C, else any."""
    available = sorted({s.symbol for s in shieldings},
                       key=lambda e: {"H": 0, "C": 1}.get(e, 2))
    if not available:
        raise CliError("no NMR shieldings in this file")
    if wanted is None:
        return available[0]
    if wanted not in available:
        raise CliError(f"no {wanted} shieldings in this file "
                       f"(it has {', '.join(available)})")
    return wanted


def build_spectrum(kind: str, result, args):
    """Turn the parsed result into a Spectrum, with the flags applied."""
    from .spectra import ir_spectrum, nmr_spectrum, uvvis_spectrum

    fwhm = args.fwhm if args.fwhm is not None else DEFAULT_FWHM[kind]
    if fwhm <= 0:
        raise CliError("--fwhm must be positive")
    try:
        if kind == "ir":
            if not result.vibrations:
                raise CliError("no vibrational frequencies in this file "
                               "(it needs a freq job)")
            return ir_spectrum(result.vibrations, scale=args.freq_scale,
                               fwhm=fwhm)
        if kind == "uv":
            if not result.excited_states:
                raise CliError("no excited states in this file "
                               "(it needs a TD-DFT/CIS job)")
            return uvvis_spectrum(result.excited_states, fwhm_ev=fwhm,
                                  unit=args.unit)
        nucleus = pick_nucleus(result.nmr_shieldings, args.nucleus)
        return nmr_spectrum(result.nmr_shieldings, element=nucleus,
                            reference=args.reference, fwhm=fwhm)
    except ValueError as exc:            # the builders' own complaints
        raise CliError(str(exc)) from None


def plot_spectrum(args) -> str:
    """Silent mode with --ir/--uv/--nmr: write the spectrum figure."""
    kind = spectrum_kind(args)
    if not args.file:
        raise CliError(f"give a file to plot, e.g. scope -s --{kind} mycalc.log")
    if not Path(args.file).is_file():
        raise FileNotFoundError(args.file)
    out = check_output(args.output or default_output(args.file, kind), kind)

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    from .spectra.plot import draw_spectrum, set_xrange

    result = mio.load(args.file)
    spec = build_spectrum(kind, result, args)

    width, height = resolve_size(args, kind)
    dpi = max(int(args.dpi), 1)
    figure = Figure(figsize=(width / dpi, height / dpi), dpi=dpi,
                    tight_layout=True)
    FigureCanvasAgg(figure)              # savefig needs a canvas to start from
    ax = figure.add_subplot(111)
    # UV-Vis oscillator strengths get their own axis, as in the viewer
    stick_axis = ax.twinx() if kind == "uv" else None
    draw_spectrum(ax, spec, sticks=not args.no_sticks, stick_axis=stick_axis,
                  stick_label="oscillator strength f")
    if args.xrange:
        set_xrange(ax, spec, *parse_range(args.xrange))
    if args.title:
        ax.set_title(args.title)

    figure.savefig(out, transparent=not args.opaque)
    if args.csv:
        import pandas as pd
        pd.DataFrame({spec.xlabel: spec.x, spec.ylabel: spec.y}).to_csv(
            args.csv, index=False)
    return out



def _ensure_qt_app():
    """A QGuiApplication is needed even with no window (GL context, QImage)."""
    global _qt_app
    # a headless Linux box has no display; macOS needs its real platform plugin
    if sys.platform.startswith("linux") and not os.environ.get("DISPLAY") \
            and not os.environ.get("WAYLAND_DISPLAY"):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtGui import QGuiApplication
    _qt_app = QGuiApplication.instance() or QGuiApplication([sys.argv[0]])
    return _qt_app


def render(args) -> str:
    """Silent mode: render *args.file* to an image and return the path."""
    if not args.file:
        raise CliError("give a file to render, e.g. scope -s mycalc.log")
    if not Path(args.file).is_file():          # before Qt starts up
        raise FileNotFoundError(args.file)
    out = check_output(args.output or default_output(args.file))

    _ensure_qt_app()
    from PySide6.QtGui import QPainter
    from .gui.annotations import draw_annotations
    from .render.offscreen import render_molecule_image

    result = mio.load(args.file)
    molecule = pick_frame(result.frames, args.frame)
    if molecule.bonds is None:
        molecule.perceive_bonds()

    style = build_style(args)
    camera = build_camera(molecule, args.view, args.align, args.rotate, args.zoom)
    volume = load_volume(args, result)
    reps = parse_reps(args.rep, molecule.natoms) if args.rep else None
    measured = [parse_indices(spec, molecule.natoms) for spec in args.measure or []]
    for idxs in measured:
        if not 2 <= len(idxs) <= 4:
            raise CliError("--measure wants 2 atoms (distance), 3 (angle) "
                           "or 4 (dihedral)")
    width, height = resolve_size(args)

    try:
        image = render_molecule_image(
            molecule, style, camera, width, height,
            supersample=args.supersample, transparent=not args.opaque,
            volume=volume, isovalue=args.iso, reps=reps)
    except RuntimeError as exc:
        raise CliError(
            f"{exc}\nOffscreen rendering needs a working OpenGL 3.3 driver; "
            "on a headless Linux machine try QT_QPA_PLATFORM=offscreen or run "
            "under xvfb-run.") from None

    if measured or args.labels != "none" or args.axes:
        painter = QPainter(image)
        # the export path scales annotation sizes with the render resolution
        draw_annotations(painter, molecule, camera, image.width(), image.height(),
                         measured, args.labels,
                         scale=max(1.0, height / 800.0), axes=args.axes)
        painter.end()

    write_image(image, out)
    return out


def launch_gui(args) -> int:
    from .gui.app import main as gui_main
    return gui_main(args.file, style=args.style,
                    label_mode=args.labels, axes=args.axes)


# ------------------------------------------------------------------ entry point

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scope",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Molecular structure and spectroscopy viewer.",
        epilog="""examples:
  scope mycalc.log                       open the viewer
  scope -s mycalc.log                    render mycalc.png and exit
  scope -s opt.log -o fig.tif --style houk --size 2000x1500
  scope -s mol.xyz --view top --labels number --measure 3,4
  scope -s mol.log --rep '1-12:ball;13-40:line' --view 30,-15
  scope -s homo.cube --iso 0.03 --iso-colors purple,gold --axes
  scope -s freq.log --ir -o ir.pdf --fwhm 12 --freq-scale 0.965
  scope -s td.log --uv --unit eV --xrange 2:7
  scope -s nmr.log --nmr --nucleus C --reference 186.4 --csv shifts.csv

Atom numbers are 1-based, as shown in the viewer.""")
    parser.add_argument("file", nargs="?", help="structure/output file to open")
    parser.add_argument("-s", "--silent", action="store_true",
                        help="do not open a window: render to an image and exit")
    parser.add_argument("-o", "--output", metavar="PATH",
                        help="image to write (.png or .tif; default: <file>.png)")
    parser.add_argument("--version", action="version",
                        version=f"microscope {__version__}")

    look = parser.add_argument_group("appearance")
    look.add_argument("--style", default="cylview", choices=("cylview", "houk"),
                      help="rendering style (default: cylview)")
    look.add_argument("--rep", metavar="SPEC",
                      help="representation per atom range, e.g. "
                           "'stick' or '1-12:ball;13-40:line'")
    look.add_argument("--labels", default="none", choices=LABEL_MODES,
                      help="atom labels (default: none)")
    look.add_argument("--background", metavar="COLOR",
                      help="background colour (name or #rrggbb)")
    look.add_argument("--bond-color", metavar="COLOR",
                      help="draw every bond in one colour instead of split")
    look.add_argument("--no-hbonds", action="store_true",
                      help="hide the dashed hydrogen bonds")
    look.add_argument("--axes", action="store_true",
                      help="draw the XYZ orientation triad")

    view = parser.add_argument_group("view")
    view.add_argument("--view", metavar="SPEC",
                      help="front|back|top|bottom|left|right, or two angles "
                           "like 30,-15 (degrees about the screen axes)")
    view.add_argument("--rotate", metavar="RX,RY",
                      help="turn the view by these angles after --view/--align")
    view.add_argument("--align", metavar="ATOMS",
                      help="2 atoms: look down that bond; 3 atoms: put their "
                           "plane in the screen")
    view.add_argument("--zoom", type=float, default=1.0, metavar="FACTOR",
                      help="magnify by this factor (default: 1.0)")

    image = parser.add_argument_group("image (silent mode)")
    image.add_argument("--size", metavar="WxH",
                       help="pixel size (default: 1200x900 for a structure, "
                            "1600x900 for a spectrum)")
    image.add_argument("--supersample", type=int, default=3, metavar="N",
                       help="anti-aliasing oversampling (default: 3)")
    image.add_argument("--opaque", action="store_true",
                       help="fill the background instead of leaving it transparent")
    image.add_argument("--measure", action="append", metavar="ATOMS",
                       help="annotate a distance/angle/dihedral, e.g. 3,4 "
                            "(repeatable)")

    spectra = parser.add_argument_group(
        "spectra (silent mode; picks a spectrum instead of a structure)")
    which = spectra.add_mutually_exclusive_group()
    which.add_argument("--ir", action="store_true",
                       help="plot the IR spectrum of a frequency job")
    which.add_argument("--uv", "--uvvis", action="store_true", dest="uv",
                       help="plot the UV-Vis spectrum of a TD-DFT/CIS job")
    which.add_argument("--nmr", action="store_true",
                       help="plot the NMR spectrum of a shielding job")
    spectra.add_argument("--fwhm", type=float, metavar="VALUE",
                         help="line width: cm⁻¹ (IR, default 8), eV (UV, "
                              "default 0.3) or ppm (NMR, default 0.5)")
    spectra.add_argument("--freq-scale", type=float, default=1.0,
                         metavar="FACTOR",
                         help="scale IR frequencies (default: 1.0)")
    spectra.add_argument("--unit", default="nm", choices=("nm", "eV"),
                         help="UV-Vis x axis (default: nm)")
    spectra.add_argument("--nucleus", metavar="ELEMENT",
                         help="NMR nucleus to plot (default: H, else C)")
    spectra.add_argument("--reference", type=float, default=0.0, metavar="PPM",
                         help="NMR reference shielding σ₀; 0 plots raw σ "
                              "(default: 0)")
    spectra.add_argument("--xrange", metavar="MIN:MAX",
                         help="limit the x axis, e.g. 400:1800")
    spectra.add_argument("--no-sticks", action="store_true",
                         help="draw only the broadened curve")
    spectra.add_argument("--dpi", type=int, default=300, metavar="N",
                         help="figure resolution (default: 300)")
    spectra.add_argument("--title", metavar="TEXT", help="title above the plot")
    spectra.add_argument("--csv", metavar="PATH",
                         help="also write the plotted curve as CSV")

    data = parser.add_argument_group("data")
    data.add_argument("--frame", type=int, metavar="N",
                      help="trajectory frame, 1-based (default: the last one)")
    data.add_argument("--cube", metavar="FILE",
                      help="cube file to draw an isosurface from")
    data.add_argument("--mo", type=int, metavar="N",
                      help="which grid of a multi-MO cube to use (1-based)")
    data.add_argument("--iso", type=float, metavar="VALUE",
                      help="isosurface level (default: chosen from the grid)")
    data.add_argument("--iso-colors", metavar="POS,NEG",
                      help="isosurface lobe colours, e.g. blue,red")
    data.add_argument("--iso-opacity", type=float, metavar="VALUE",
                      help="isosurface opacity, 0-1")
    return parser


# flags that only mean something in silent mode, grouped by what they act on
_STRUCTURE_ONLY = ("supersample", "measure", "rotate", "view", "align", "frame",
                   "cube", "mo", "iso", "iso_colors", "iso_opacity", "rep",
                   "background", "bond_color", "no_hbonds", "zoom", "style",
                   "labels", "axes")
_SPECTRUM_ONLY = ("ir", "uv", "nmr", "fwhm", "freq_scale", "unit", "nucleus",
                  "reference", "xrange", "no_sticks", "dpi", "title", "csv")
# these three configure the viewer too, so they are allowed without -s
_GUI_ALLOWED = ("style", "labels", "axes")
_DEFAULTS = {"zoom": 1.0, "supersample": 3, "style": "cylview",
             "labels": "none", "freq_scale": 1.0, "unit": "nm", "reference": 0.0,
             "dpi": 300}


def _used(args, names) -> list[str]:
    return [n for n in names
            if getattr(args, n, None) not in (None, False, _DEFAULTS.get(n))]


def _flags(names) -> str:
    return ", ".join("--" + n.replace("_", "-") for n in names)


def check_gui_flags(args, parser) -> None:
    """Refuse silently-ignored flags rather than pretending they worked."""
    silent_only = ("output", "size", "opaque")     # shared by both silent modes
    used = _used(args, [n for n in _STRUCTURE_ONLY if n not in _GUI_ALLOWED]
                 + list(_SPECTRUM_ONLY) + list(silent_only))
    if used:
        verb = "only applies" if len(used) == 1 else "only apply"
        parser.error(f"{_flags(used)} {verb} with -s/--silent "
                     "(without it the viewer window opens)")


def check_mode_flags(args, parser) -> None:
    """In silent mode, keep structure flags and spectrum flags apart."""
    if spectrum_kind(args):
        stray = _used(args, _STRUCTURE_ONLY)
        if stray:
            verb = "does not apply" if len(stray) == 1 else "do not apply"
            parser.error(f"{_flags(stray)} {verb} to a spectrum plot")
    else:
        stray = _used(args, _SPECTRUM_ONLY)
        if stray:
            verb = "needs" if len(stray) == 1 else "need"
            parser.error(f"{_flags(stray)} {verb} one of --ir, --uv or --nmr")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.silent:
            check_mode_flags(args, parser)
            print(plot_spectrum(args) if spectrum_kind(args) else render(args))
            return 0
        check_gui_flags(args, parser)
        return launch_gui(args)
    except CliError as exc:
        print(f"scope: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"scope: no such file: {exc.filename or exc}", file=sys.stderr)
        return 2
    except (FileFormatError, UnsupportedFormatError, OSError) as exc:
        print(f"scope: could not read the file: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
