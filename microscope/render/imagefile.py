"""Writing rendered images to disk.

Only formats that can keep a transparent background are supported, which is
what publication figures need: PNG, or uncompressed TIFF as a lossless master.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QMarginsF, QRect, QSize, QSizeF
from PySide6.QtGui import QImage, QImageWriter, QPageSize, QPdfWriter


def write_image(image: QImage, path) -> None:
    """Save as PNG or TIFF — both keep the transparent background intact."""
    writer = QImageWriter(str(path))
    if Path(path).suffix.lower() in (".tif", ".tiff"):
        writer.setCompression(0)                # lossless master copy
    if not writer.write(image):
        raise RuntimeError(writer.errorString() or f"could not write {path}")


# Line art is worth keeping as line art: SVG for an illustration program, PDF
# for a manuscript. One unit is one point in both.

def vector_device(path, width: int, height: int, title: str = "figure"):
    """A QPainter device writing *path* as SVG or PDF, *width* x *height* points."""
    if Path(path).suffix.lower() == ".svg":
        from PySide6.QtSvg import QSvgGenerator  # QtSvg is a separate module

        generator = QSvgGenerator()
        generator.setFileName(str(path))
        generator.setSize(QSize(int(width), int(height)))
        generator.setViewBox(QRect(0, 0, int(width), int(height)))
        generator.setTitle(title)
        generator.setDescription("Drawn by microscope")
        return generator
    writer = QPdfWriter(str(path))
    writer.setResolution(72)
    writer.setPageSize(QPageSize(QSizeF(width, height), QPageSize.Unit.Point,
                                 "figure", QPageSize.SizeMatchPolicy.ExactMatch))
    writer.setPageMargins(QMarginsF(0, 0, 0, 0))
    return writer
