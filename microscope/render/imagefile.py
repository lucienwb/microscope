"""Writing rendered images to disk.

Only formats that can keep a transparent background are supported, which is
what publication figures need: PNG, or uncompressed TIFF as a lossless master.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QImage, QImageWriter

IMAGE_SUFFIXES = (".png", ".tif", ".tiff")


def write_image(image: QImage, path) -> None:
    """Save as PNG or TIFF — both keep the transparent background intact."""
    writer = QImageWriter(str(path))
    if Path(path).suffix.lower() in (".tif", ".tiff"):
        writer.setCompression(0)                # lossless master copy
    if not writer.write(image):
        raise RuntimeError(writer.errorString() or f"could not write {path}")
