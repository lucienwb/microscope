"""Which extensions each kind of output is written to.

Knowing the names is not the same as being able to paint them: the command
line has to check an output path before it decides whether it needs graphics
at all, so this module stays free of Qt and matplotlib and the writers import
their own suffix list from here.
"""

from __future__ import annotations

IMAGE_SUFFIXES = (".png", ".tif", ".tiff")      # raster, written through Qt
VECTOR_SUFFIXES = (".svg", ".pdf")              # line art, also through Qt
PLOT_SUFFIXES = (".png", ".pdf", ".svg", ".eps", ".tif", ".tiff")   # matplotlib
