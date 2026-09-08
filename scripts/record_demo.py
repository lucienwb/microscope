"""Record the animated GIFs the README uses.

    python scripts/record_demo.py lewis      # no OpenGL needed
    python scripts/record_demo.py orbit vibrate

The Lewis clips are painted with QPainter and record anywhere, including a
headless box. The 3-D clips need a real OpenGL 3.3 context, so they have to
be recorded on a machine with a working driver.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from PIL import Image  # noqa: E402
from PySide6.QtCore import QBuffer  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from microscope.gui.app import MainWindow  # noqa: E402

OUT = REPO / "docs"
COLORS = 32          # line art needs very few; keeps the files small


def _frame(widget) -> Image.Image:
    buf = QBuffer()
    buf.open(QBuffer.OpenModeFlag.ReadWrite)
    widget.grab().save(buf, "PNG")
    image = Image.open(io.BytesIO(buf.data().data())).convert("RGB")
    return image.convert("P", palette=Image.ADAPTIVE, colors=COLORS)


def _save(frames, name: str, ms: int = 80) -> None:
    if len(frames[0].convert("RGB").getcolors(maxcolors=8) or [0] * 9) < 3:
        raise SystemExit(
            f"{name}: every frame came out blank. The 3-D clips need a working "
            "OpenGL 3.3 driver — record them on a machine that has one.")
    path = OUT / name
    frames[0].save(path, save_all=True, append_images=frames[1:],
                   duration=ms, loop=0, optimize=True)
    print(f"{path.relative_to(REPO)}  {len(frames)} frames  "
          f"{path.stat().st_size / 1024:.0f} KB")


def lewis(app, window):
    """Turn a Lewis structure in the plane of the page."""
    window.open_file(REPO / "tests/data/trp.log")
    window._lewis_action.setChecked(True)
    app.processEvents()
    window._reset_view()
    window.lewis.camera.rotate_drag(-80, -80)
    frames = []
    for _ in range(36):
        window.lewis.camera.spin(10)
        app.processEvents()
        frames.append(_frame(window.lewis))
    _save(frames, "lewis_spin.gif")


def orbit(app, window):
    """Turn the 3-D model. Needs a working OpenGL driver."""
    window.open_file(REPO / "tests/data/trp.log")
    app.processEvents()
    frames = []
    for _ in range(48):
        window.viewport.camera.rotate_drag(14, 0)
        window.viewport.update()
        app.processEvents()
        frames.append(_frame(window.viewport))
    _save(frames, "orbit.gif", ms=60)


def vibrate(app, window):
    """Play a normal mode, the way clicking an IR band does. Needs OpenGL."""
    window.open_file(REPO / "tests/data/dvb_ir.out")
    app.processEvents()
    result = window.result
    mode = max(result.vibrations, key=lambda v: v.ir_intensity)
    window.viewport.animate_mode(mode.displacements)
    frames = []
    for _ in range(40):
        window.viewport._anim_tick()
        app.processEvents()
        frames.append(_frame(window.viewport))
    window.viewport.stop_animation()
    _save(frames, "vibration.gif", ms=50)


CLIPS = {"lewis": lewis, "orbit": orbit, "vibrate": vibrate}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clips", nargs="*", default=["lewis"],
                        choices=[*CLIPS, []], help="which clips to record")
    args = parser.parse_args()

    app = QApplication(sys.argv[:1])
    window = MainWindow()
    window.resize(660, 460)
    window.show()
    for name in args.clips or ["lewis"]:
        CLIPS[name](app, window)
    return 0


if __name__ == "__main__":
    sys.exit(main())
