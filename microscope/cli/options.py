"""The argument parser, and which flags may be used together.

Rejecting a render-only flag in GUI mode, or a spectrum flag while drawing a
structure, beats silently ignoring it.
"""

from __future__ import annotations

import argparse

from .. import __version__
from .lewis import LEWIS_OPTION_FLAGS
from .parsing import LABEL_MODES, spectrum_kind


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
  scope -s mol.log --lewis -o scheme.svg --align 2,3,4
  scope -s mol.log --lewis -o mol.cdxml     (opens in ChemDraw at that angle)
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
    view.add_argument("--spin", type=float, default=0.0, metavar="DEG",
                      help="turn the picture in the plane of the page")
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

    flat = parser.add_argument_group(
        "Lewis structure (the flat, ChemDraw-style drawing)")
    flat.add_argument("--lewis", action="store_true",
                      help="draw the Lewis structure instead of the 3D model; "
                           "-o may be an image (.png/.tif), vector art "
                           "(.svg/.pdf) or a structure file (.cdxml/.mol)")
    flat.add_argument("--show-hydrogens", action="store_true",
                      help="draw every hydrogen instead of folding it into a label")
    flat.add_argument("--carbon-labels", action="store_true",
                      help="write C on carbons instead of leaving a bare vertex")
    flat.add_argument("--lone-pairs", action="store_true",
                      help="mark non-bonding electron pairs with dots")
    flat.add_argument("--color-atoms", action="store_true",
                      help="tint labels by element instead of drawing in black")

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
_LEWIS_ONLY = LEWIS_OPTION_FLAGS
_STRUCTURE_ONLY = ("supersample", "measure", "rotate", "view", "align", "frame",
                   "cube", "mo", "iso", "iso_colors", "iso_opacity", "rep",
                   "background", "bond_color", "no_hbonds", "zoom", "spin",
                   "style", "labels", "axes", "lewis") + _LEWIS_ONLY
_SPECTRUM_ONLY = ("ir", "uv", "nmr", "fwhm", "freq_scale", "unit", "nucleus",
                  "reference", "xrange", "no_sticks", "dpi", "title", "csv")
# a Lewis drawing is line art: none of the 3D renderer's knobs reach it
_LEWIS_INCOMPATIBLE = ("style", "rep", "bond_color", "no_hbonds", "labels",
                       "axes", "measure", "supersample", "cube", "mo", "iso",
                       "iso_colors", "iso_opacity")
# these configure the viewer too, so they are allowed without -s
_GUI_ALLOWED = ("style", "labels", "axes", "lewis") + _LEWIS_ONLY
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
        return
    stray = _used(args, _SPECTRUM_ONLY)
    if stray:
        verb = "needs" if len(stray) == 1 else "need"
        parser.error(f"{_flags(stray)} {verb} one of --ir, --uv or --nmr")
    check_lewis_flags(args, parser)


def check_lewis_flags(args, parser) -> None:
    """The flat drawing has its own switches, and ignores the renderer's."""
    if args.lewis:
        stray = _used(args, _LEWIS_INCOMPATIBLE)
        if stray:
            verb = "does not apply" if len(stray) == 1 else "do not apply"
            parser.error(f"{_flags(stray)} {verb} to a Lewis drawing — it is "
                         "line art, not a render")
        return
    stray = _used(args, _LEWIS_ONLY)
    if stray:
        verb = "needs" if len(stray) == 1 else "need"
        parser.error(f"{_flags(stray)} {verb} --lewis")
