"""`scope -s file.log`: render the molecule itself, offscreen through OpenGL."""

from __future__ import annotations

from pathlib import Path

from .. import io as mio
from ..render.imagefile import write_image
from .parsing import (
    CliError,
    build_camera,
    build_style,
    check_output,
    default_output,
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
    out = check_output(args.output or default_output(args.file))

    from .main import _ensure_qt_app
    _ensure_qt_app()
    from PySide6.QtGui import QPainter

    from ..gui.annotations import draw_annotations
    from ..render.offscreen import render_molecule_image

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
