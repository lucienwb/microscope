"""`scope -s file.log --lewis`: the flat ChemDraw-style drawing.

Writing a .cdxml or .mol needs no graphics at all, and this deliberately
returns before Qt is ever started for those.
"""

from __future__ import annotations

from pathlib import Path

from .. import io as mio
from ..render.formats import VECTOR_SUFFIXES
from .parsing import (
    CliError,
    build_camera,
    check_output,
    default_output,
    parse_color,
    pick_frame,
    resolve_size,
)

# the argparse destinations of the four Lewis toggles, named once so the
# option checker and the drawing cannot drift apart
LEWIS_OPTION_FLAGS = ("show_hydrogens", "carbon_labels", "lone_pairs",
                      "color_atoms")

def build_lewis_options(args):
    """The four switches that shape the drawing."""
    from ..render.lewis2d import LewisOptions
    return LewisOptions(show_hydrogens=args.show_hydrogens,
                        carbon_labels=args.carbon_labels,
                        lone_pairs=args.lone_pairs,
                        color_atoms=args.color_atoms)


def lewis_option_dict(args) -> dict:
    return {name: bool(getattr(args, name)) for name in LEWIS_OPTION_FLAGS}


def render_lewis(args) -> str:
    """Silent mode with --lewis: write the flat drawing, or a structure file.

    The drawing is the molecule seen through the camera the view flags set up,
    so --view/--align/--rotate/--spin choose the angle it is drawn from, the
    same choice the interactive Lewis view is there to make.
    """
    if not args.file:
        raise CliError("give a file to draw, e.g. scope -s --lewis mycalc.log")
    if not Path(args.file).is_file():          # before Qt starts up
        raise FileNotFoundError(args.file)
    out = check_output(args.output or default_output(args.file, "lewis"), "lewis")

    from ..core.lewis import perceive
    from ..render.lewis2d import hidden_hydrogens

    result = mio.load(args.file)
    molecule = pick_frame(result.frames, args.frame)
    if molecule.bonds is None:
        molecule.perceive_bonds()
    structure = perceive(molecule)
    camera = build_camera(molecule, args.view, args.align, args.rotate,
                          args.zoom, args.spin)
    options = build_lewis_options(args)
    width, height = resolve_size(args, "lewis")
    points = camera.project(molecule.coords, width, height)

    suffix = Path(out).suffix.lower()
    if suffix in mio.LEWIS_EXTENSIONS:         # a structure file needs no Qt
        hidden = hidden_hydrogens(structure, options)
        mio.save_lewis(out, structure, points,
                       [a for a in range(molecule.natoms) if not hidden[a]])
        return out

    from .main import _ensure_qt_app
    _ensure_qt_app()
    from ..render.imagefile import write_image
    from ..render.lewisdraw import render_lewis_image, write_lewis_vector

    background = parse_color(args.background) if args.background else (1.0, 1.0, 1.0)
    if suffix in VECTOR_SUFFIXES:
        write_lewis_vector(out, structure, camera, width, height, options,
                           None if not args.opaque else background)
    else:
        write_image(render_lewis_image(structure, camera, width, height, options,
                                       transparent=not args.opaque,
                                       background=background), out)
    return out
