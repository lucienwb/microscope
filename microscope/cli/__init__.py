"""The `scope` command line.

    scope mycalc.log                 open the viewer window
    scope -s mycalc.log              render mycalc.png and exit, no window
    scope -s mycalc.log --ir         plot the IR spectrum instead
    scope -s mycalc.log --lewis      draw the flat Lewis structure

Silent mode is the gnuplot-style half of the program: one command line turns a
quantum-chemistry output into a finished figure, so figures can be regenerated
from a shell script or a Makefile after every re-optimization. Everything the
GUI can put in an exported image — style, view, labels, representations,
measurements, cube isosurfaces, transparency — has a flag here, and
--ir/--uv/--nmr switch from drawing the molecule to plotting its spectrum,
and --lewis switches to the flat ChemDraw-style drawing (which can be written
as a picture or as a .cdxml/.mol structure file to carry on editing).

Atom numbers on the command line are 1-based, matching the labels shown in the
viewer (C1, C2, …).

The work is split by what it needs: :mod:`parsing` turns text into values and
needs nothing, :mod:`options` declares the flags, and :mod:`structure`,
:mod:`spectrum` and :mod:`lewis` are the three things `-s` can produce.
"""

from __future__ import annotations

from .main import main
from .options import build_parser
from .parsing import CliError

__all__ = ["main", "build_parser", "CliError"]

# Everything else is reached through the submodule that owns it. Resolving it
# on demand (PEP 562) keeps `from microscope.cli import ...` convenient while
# leaving the modules that need Qt or matplotlib unimported until something
# actually asks to draw: `scope -s --lewis -o x.cdxml` is pure text and must
# not pay for graphics.
_WHERE = {
    "lewis": ("build_lewis_options", "lewis_option_dict", "render_lewis",
              "LEWIS_OPTION_FLAGS"),
    "main": ("launch_gui",),
    "options": ("check_gui_flags", "check_lewis_flags", "check_mode_flags"),
    "parsing": ("DEFAULT_FWHM", "LABEL_MODES", "LEWIS_SUFFIXES", "NAMED_VIEWS",
                "PLOT_SUFFIXES", "REP_CODES", "build_camera", "build_style",
                "check_output", "default_output", "load_volume",
                "output_suffixes", "parse_color", "parse_indices",
                "parse_range", "parse_reps", "parse_rotation", "parse_size",
                "pick_frame", "resolve_size", "spectrum_kind"),
    "spectrum": ("build_spectrum", "pick_nucleus", "plot_spectrum"),
    "structure": ("render",),
}
_HOME = {name: module for module, names in _WHERE.items() for name in names}


def __getattr__(name):
    module = _HOME.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module
    value = getattr(import_module(f".{module}", __name__), name)
    globals()[name] = value                     # only look it up once
    return value


def __dir__():
    return sorted(set(globals()) | set(_HOME))
