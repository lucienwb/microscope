"""What the file dialogs offer, and how a chosen filter names a file.

One place to add a format: the filter string the dialog shows and the
extension to append when the user types a bare name.
"""

from __future__ import annotations

from pathlib import Path

OPEN_FILTER = (
    "Molecular files (*.xyz *.log *.out *.fchk *.fck *.fch *.gjf *.com *.gau "
    "*.pdb *.molden *.cube *.cub);;All files (*)"
)
SAVE_FILTER = ("XYZ (*.xyz);;Gaussian input (*.gjf *.com);;"
               "ORCA input (*.inp);;Q-Chem input (*.in *.qcin);;PDB (*.pdb);;"
               "MDL molfile (*.mol)")
_FILTER_DEFAULT_EXT = {"XYZ": ".xyz", "Gaussian": ".gjf", "ORCA": ".inp",
                       "Q-Chem": ".in", "PDB": ".pdb", "MDL": ".mol"}

EXPORT_FILTER = "PNG image (*.png);;TIFF image, uncompressed (*.tif *.tiff)"
# the flat drawing is line art, so it can also be kept as vector art
LEWIS_EXPORT_FILTER = ("PNG image (*.png);;SVG drawing (*.svg);;PDF (*.pdf);;"
                       "TIFF image, uncompressed (*.tif *.tiff)")
_EXPORT_DEFAULT_EXT = {"PNG": ".png", "TIFF": ".tif", "SVG": ".svg", "PDF": ".pdf"}

CHEMDRAW_FILTER = "ChemDraw (*.cdxml);;MDL molfile (*.mol)"
_CHEMDRAW_DEFAULT_EXT = {"ChemDraw": ".cdxml", "MDL": ".mol"}


def with_extension(path: str, chosen: str, defaults: dict) -> str:
    """Append the extension of the selected filter when none was typed."""
    if Path(path).suffix:
        return path
    for name, ext in defaults.items():
        if chosen.startswith(name):
            return path + ext
    return path

STYLE_FILTER = "microscope style (*.json);;All files (*)"
