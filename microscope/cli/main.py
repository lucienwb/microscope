"""The `scope` entry point: work out what was asked for, then do it.

Each of the three things it can do lives in its own module; this one only
decides which, and starts Qt at the last possible moment.
"""

from __future__ import annotations

import os
import sys

from ..io.errors import FileFormatError, UnsupportedFormatError
from .options import build_parser, check_gui_flags, check_lewis_flags, check_mode_flags
from .parsing import CliError, spectrum_kind

# Kept alive for the length of the process: a QGuiApplication that is garbage
# collected takes the GL context with it.
_qt_app = None


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


def launch_gui(args) -> int:
    from ..gui.app import main as gui_main
    from .lewis import lewis_option_dict
    return gui_main(args.file, style=args.style, label_mode=args.labels,
                    axes=args.axes, lewis=args.lewis,
                    lewis_options=lewis_option_dict(args))


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.silent:
            check_mode_flags(args, parser)
            # imported here, not at the top: plotting a spectrum should not
            # pay for OpenGL, and writing a .cdxml should not pay for Qt
            if spectrum_kind(args):
                from .spectrum import plot_spectrum
                print(plot_spectrum(args))
            elif args.lewis:
                from .lewis import render_lewis
                print(render_lewis(args))
            else:
                from .structure import render
                print(render(args))
            return 0
        check_gui_flags(args, parser)
        check_lewis_flags(args, parser)
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
