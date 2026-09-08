"""Starting the viewer.

The window itself lives in :mod:`microscope.gui.mainwindow`; this is only the
entry point the CLI calls, kept apart so opening a window is one readable
function.
"""

from __future__ import annotations

import sys

from PySide6.QtGui import QSurfaceFormat
from PySide6.QtWidgets import QApplication

from .mainwindow import MainWindow

__all__ = ["MainWindow", "main"]


def main(path: str | None = None, style: str = "cylview",
         label_mode: str = "none", axes: bool = False, lewis: bool = False,
         lewis_options: dict | None = None) -> int:
    """Open the viewer window. Returns the Qt exit code.

    The optional arguments are the handful of `scope` flags that make sense
    for an interactive session; the rest are silent-mode only.
    """
    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    fmt.setSamples(8)
    QSurfaceFormat.setDefaultFormat(fmt)

    app = QApplication(sys.argv[:1])
    app.setApplicationName("microscope")
    window = MainWindow()
    window.resize(1000, 720)
    window.show()
    if style != "cylview":
        window._set_representation(style)
    if label_mode != "none":
        window.viewport.set_label_mode(label_mode)
        window._label_actions[label_mode].setChecked(True)
    if axes:
        window._axes_action.setChecked(True)      # toggled -> viewport
    if path:
        window.open_file(path)
    for name, value in (lewis_options or {}).items():
        window._lewis_options[name].setChecked(bool(value))
    if lewis:
        window._lewis_action.setChecked(True)     # needs the file open first
    return app.exec()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
