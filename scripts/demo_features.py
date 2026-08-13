"""Guided demo of the viewer features.

Opens divinylbenzene with:
  - the benzene ring aligned flat to the screen,
  - atom labels on (element+number),
  - a pinned C–C distance and a pinned ring angle,
then leaves the app running for you to explore.

Usage:
    python scripts/demo_features.py                   # interactive
    python scripts/demo_features.py --screenshot x.png  # render setup and exit
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from PySide6.QtCore import QTimer  # noqa: E402
from PySide6.QtGui import QSurfaceFormat  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from microscp.gui.app import MainWindow  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "tests" / "data" / "dvb_ir.out"


def _neighbors(mol) -> dict[int, list[int]]:
    neighbors: dict[int, list[int]] = {}
    for i, j in mol.bonds:
        neighbors.setdefault(int(i), []).append(int(j))
        neighbors.setdefault(int(j), []).append(int(i))
    return neighbors


def ring_triple(mol):
    """Find a bonded C–C–C triple (vertex with two carbon neighbors)."""
    zs = mol.atomic_numbers
    for k, nbrs in _neighbors(mol).items():
        cn = [n for n in nbrs if zs[n] == 6]
        if zs[k] == 6 and len(cn) >= 2:
            return [cn[0], k, cn[1]]
    return [0, 1, 2]


def main():
    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    fmt.setSamples(8)
    QSurfaceFormat.setDefaultFormat(fmt)

    app = QApplication(sys.argv)
    win = MainWindow()
    win.resize(1100, 800)
    win.show()
    win.open_file(str(DATA))

    def setup():
        vp = win.viewport
        mol = vp.molecule
        if mol is None:
            return
        # ring plane -> screen
        triple = ring_triple(mol)
        vp.selection = list(triple)
        vp.align_to_selection()
        # pin the ring angle, then a C-C bond distance far from it (no overlap)
        vp.pin_selection()
        zs = mol.atomic_numbers
        heavy = {k: sum(1 for n in nb if zs[n] > 1)
                 for k, nb in _neighbors(mol).items()}
        vertex = mol.coords[triple[1]]
        i, j = max(((int(i), int(j)) for i, j in mol.bonds
                    if zs[i] == 6 and zs[j] == 6
                    and heavy[int(i)] >= 2 and heavy[int(j)] >= 2),
                   key=lambda b: np.linalg.norm(
                       (mol.coords[b[0]] + mol.coords[b[1]]) / 2.0 - vertex))
        vp.selection = [i, j]
        vp.pin_selection()
        # labels on
        vp.set_label_mode("element+number")
        win._label_actions["element+number"].setChecked(True)
        vp.update()
        win.statusBar().showMessage(
            "Demo ready — try: L (labels) · M (pin) · A (align) · C (center) · Esc",
            15000)

        if "--screenshot" in sys.argv:
            out = sys.argv[sys.argv.index("--screenshot") + 1]

            def grab():
                win.viewport.grabFramebuffer().save(out)
                print(f"wrote {out}")
                win.close()
                app.quit()

            QTimer.singleShot(400, grab)

    QTimer.singleShot(800, setup)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
