"""Opening a file off the GUI thread.

Reading a protein and finding its bonds takes a couple of seconds, and done on
the GUI thread the window cannot so much as repaint meanwhile. The part that
is plain Python and numpy - reading, bond perception, and loading matplotlib
when the file has a spectrum to show - happens on a worker thread; the window
gets the result back through a queued signal and builds the scene, which needs
the GL context, where it always did.
"""

from __future__ import annotations

from PySide6.QtCore import QThread, Signal

from .. import io as mio
from ..core.results import ParseResult


def prepare(path) -> ParseResult:
    """Everything about opening *path* that does not touch Qt."""
    result = mio.load(path)
    molecule = result.molecule
    if molecule is not None and molecule.bonds is None:
        molecule.perceive_bonds()
    if result.vibrations or result.excited_states or result.nmr_shieldings:
        # the spectra dock is about to need it; importing is not drawing
        import matplotlib.backends.backend_qtagg  # noqa: F401
        import matplotlib.figure  # noqa: F401
    return result


class Loader(QThread):
    """Runs :func:`prepare` for one path; reports back on the GUI thread."""

    loaded = Signal(object)          # the ParseResult
    failed = Signal(str)             # what went wrong, for the user

    def __init__(self, path: str, then=None, parent=None):
        super().__init__(parent)
        self.path = path
        self.then = then             # run once the file is on screen

    def run(self) -> None:
        try:
            result = prepare(self.path)
        except Exception as exc:     # reported, not raised: this is not the GUI thread
            self.failed.emit(str(exc) or type(exc).__name__)
        else:
            self.loaded.emit(result)
