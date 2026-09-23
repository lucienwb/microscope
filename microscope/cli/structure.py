"""`scope -s file.log`: render the molecule itself, offscreen through OpenGL."""

from __future__ import annotations

from pathlib import Path

from .. import io as mio
from .parsing import (
    CliError,
    build_camera,
    build_style,
    check_output,
    default_output,
    frame_figure,
    load_volume,
    parse_indices,
    parse_reps,
    pick_frame,
    resolve_size,
)


def render(args) -> str:
    """Silent mode: render *args.file* to an image and return the path."""
    if not args.file:
        raise CliError("give a file to render, e.g. scope -s mycalc.log")
    if not Path(args.file).is_file():          # before Qt starts up
        raise FileNotFoundError(args.file)
    out = args.output or default_output(args.file)
    if Path(out).suffix.lower() in mio.CUBE_EXTENSIONS:
        return write_cube(args, out)
    out = check_output(out)

    from .main import _ensure_qt_app
    _ensure_qt_app()
    from PySide6.QtGui import QPainter

    from ..gui.annotations import draw_annotations
    from ..render.imagefile import write_image
    from ..render.offscreen import render_molecule_image

    result = mio.load(args.file)
    molecule = pick_frame(result.frames, args.frame)
    if molecule.bonds is None:
        molecule.perceive_bonds()

    style = build_style(args)
    volume = load_volume(args, result)
    reps = parse_reps(args.rep, molecule.natoms) if args.rep else None
    measured = [parse_indices(spec, molecule.natoms) for spec in args.measure or []]
    for idxs in measured:
        if not 2 <= len(idxs) <= 4:
            raise CliError("--measure wants 2 atoms (distance), 3 (angle) "
                           "or 4 (dihedral)")
    width, height = resolve_size(args)
    # framed on what the figure shows from the chosen angle; --zoom then scales that
    camera = build_camera(
        molecule, args.view, args.align, args.rotate, args.zoom,
        fit=lambda cam: frame_figure(cam, molecule, style, width / height, reps,
                                     volume, args.iso))

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


def write_cube(args, out: str) -> str:
    """`-o homo.cube`: write the grid itself for another program, drawing
    nothing - so no Qt and no OpenGL."""
    result = mio.load(args.file)
    volume = load_volume(args, result)
    if volume is None:
        raise CliError("a .cube output is the grid of an orbital: choose one "
                       "with --mo, e.g. --mo homo")
    mio.save_cube(out, volume, result.molecule)
    return out
