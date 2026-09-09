"""microscope — molecular structure and spectroscopy viewer.

The names below are the supported surface: read a calculation, look at what
came out of it, write it somewhere else. Everything under ``microscope.core``
is plain numpy and safe to use headless; ``microscope.render`` and
``microscope.gui`` need graphics and are not imported here, so a script that
only wants the data never pays for Qt.

    >>> import microscope
    >>> result = microscope.load("opt.log")
    >>> result.molecule.formula()
    'C10H10'

Anything not listed in ``__all__`` is internal and may move between releases.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .core.lewis import LewisStructure, perceive
from .core.molecule import Molecule
from .core.results import ExcitedState, NMRShielding, ParseResult, Vibration
from .core.volume import VolumeData
from .io import load, save_lewis, save_molecule
from .io.errors import FileFormatError, UnsupportedFormatError
from .render.styles import Style, StyleError, load_style, save_style

__all__ = [
    "__version__",
    # reading and writing
    "load", "save_molecule", "save_lewis",
    # what you get back
    "Molecule", "ParseResult", "Vibration", "ExcitedState", "NMRShielding",
    "VolumeData",
    # Lewis-structure perception
    "perceive", "LewisStructure",
    # rendering styles, savable as files a group can share
    "Style", "load_style", "save_style",
    # errors worth catching
    "FileFormatError", "UnsupportedFormatError", "StyleError",
]
