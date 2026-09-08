"""The menu bar, as one readable table.

Kept out of the window because it is a different kind of code: what the
program can do and which key does it, rather than how any of it works. Adding
a command means one entry here and one method on the window.
"""

from __future__ import annotations

from PySide6.QtGui import QAction, QActionGroup, QKeySequence

from ..render.scene import REP_BALL, REP_LINE, REP_STICK


def build_menus(window):
    file_menu = window.menuBar().addMenu("&File")
    add_action(window, file_menu, "&Open…", QKeySequence.StandardKey.Open, window.open_dialog)
    add_action(window, file_menu, "&Save As…", QKeySequence.StandardKey.Save, window.save_dialog)
    add_action(window, file_menu, "&Export Image…", "Ctrl+E", window.export_image)
    add_action(window, file_menu, "Save for &ChemDraw…", "Ctrl+Shift+S",
                     window.save_chemdraw)
    file_menu.addSeparator()
    add_action(window, file_menu, "&Quit", QKeySequence.StandardKey.Quit, window.close)

    edit_menu = window.menuBar().addMenu("&Edit")
    # everything in here is greyed out while the Lewis structure is showing
    window._edit_actions = [
        add_action(window, edit_menu, "&Undo", QKeySequence.StandardKey.Undo,
                         window._undo),
        add_action(window, edit_menu, "&Redo", QKeySequence.StandardKey.Redo,
                         window._redo),
    ]
    edit_menu.addSeparator()
    window._edit_actions += [
        add_action(window, edit_menu, "&Adjust Selection…", "E",
                         window._adjust_selection),
        add_action(window, edit_menu, "Select Connected &Fragment", "F",
                         window._select_fragment),
    ]
    window._gizmo_action = QAction("Show &Move/Rotate Handles", window, checkable=True)
    window._gizmo_action.setChecked(window.viewport.show_gizmo)
    window._gizmo_action.setShortcut(QKeySequence("G"))
    window._gizmo_action.setToolTip(
        "Drag the coloured arrows to slide the selection, the rings to turn "
        "it, the centre dot to move it in the screen plane (Shift snaps)")
    window._gizmo_action.toggled.connect(window._toggle_gizmo)
    edit_menu.addAction(window._gizmo_action)
    window._edit_actions.append(window._gizmo_action)
    edit_menu.addSeparator()
    window._edit_actions += [
        add_action(window, edit_menu, "Delete Selected Atoms", "X",
                         window._delete_atoms),
        add_action(window, edit_menu, "Recompute &Bonds", "B",
                         window._recompute_bonds),
    ]

    view_menu = window.menuBar().addMenu("&View")
    add_action(window, view_menu, "&Reset View", "Ctrl+R", window._reset_view)
    view_menu.addSeparator()

    repr_menu = view_menu.addMenu("Re&presentation")
    repr_group = QActionGroup(window)
    window._repr_actions = {}
    for name, text in (("cylview", "&CYLview"), ("houk", "&Houk (Houkmol)")):
        action = QAction(text, window, checkable=True)
        action.triggered.connect(lambda _=False, n=name: window._set_representation(n))
        repr_group.addAction(action)
        repr_menu.addAction(action)
        window._repr_actions[name] = action
    window._repr_actions["cylview"].setChecked(True)
    add_action(window, repr_menu, "&Toggle Representation", "V",
                     window._toggle_representation)
    repr_menu.addSeparator()
    for rep, text, key in ((REP_BALL, "Selection → &Ball && Stick", "1"),
                           (REP_STICK, "Selection → &Stick", "2"),
                           (REP_LINE, "Selection → &Line", "3")):
        # swallow QAction.triggered's checked argument (it would land in r)
        add_action(window, repr_menu, text, key,
                         lambda _=False, r=rep: window._apply_atom_rep(r))

    labels_menu = view_menu.addMenu("Atom &Labels")
    group = QActionGroup(window)
    window._label_actions = {}
    for mode, text in (("none", "&None"), ("element", "&Element"),
                       ("element+number", "Element + Num&ber"),
                       ("number", "N&umber")):
        action = QAction(text, window, checkable=True)
        action.triggered.connect(lambda _=False, m=mode: window.viewport.set_label_mode(m))
        group.addAction(action)
        labels_menu.addAction(action)
        window._label_actions[mode] = action
    window._label_actions["none"].setChecked(True)
    add_action(window, labels_menu, "&Cycle Labels", "L", window._cycle_labels)

    window._hbond_action = QAction("Show &Hydrogen Bonds", window, checkable=True)
    window._hbond_action.setChecked(True)
    window._hbond_action.setShortcut(QKeySequence("H"))
    window._hbond_action.toggled.connect(window._toggle_hbonds)
    view_menu.addAction(window._hbond_action)

    window._axes_action = QAction("Show &XYZ Axes", window, checkable=True)
    window._axes_action.setChecked(window.viewport.show_axes)
    window._axes_action.setShortcut(QKeySequence("Shift+A"))
    window._axes_action.setToolTip(
        "Corner triad showing how the world x/y/z axes point (also exported)")
    window._axes_action.toggled.connect(window._toggle_axes)
    view_menu.addAction(window._axes_action)

    add_action(window, view_menu, "&Isosurface…", "I", window._show_surface_dialog)

    lewis_menu = view_menu.addMenu("&Lewis Structure")
    window._lewis_action = QAction("Show &Lewis Structure (2D)", window,
                                 checkable=True)
    window._lewis_action.setShortcut(QKeySequence("Shift+L"))
    window._lewis_action.setToolTip(
        "Draw the molecule flat, ChemDraw style. Turning it picks the angle "
        "the drawing is made from; the structure cannot be edited here")
    window._lewis_action.toggled.connect(window._toggle_lewis)
    lewis_menu.addAction(window._lewis_action)
    lewis_menu.addSeparator()
    window._lewis_options = {}
    for option, text in (("show_hydrogens", "Show All &Hydrogens"),
                         ("carbon_labels", "Label &Carbons"),
                         ("lone_pairs", "Show Lone &Pairs"),
                         ("color_atoms", "Colour &Atoms")):
        action = QAction(text, window, checkable=True)
        action.setChecked(getattr(window.lewis.options, option))
        action.toggled.connect(
            lambda on, name=option: window.lewis.set_option(name, on))
        lewis_menu.addAction(action)
        window._lewis_options[option] = action

    view_menu.addSeparator()
    add_action(window, view_menu, "&Align View to Selection", "A", window._align_view)
    add_action(window, view_menu, "&Center on Selected Atom", "C", window._center_atom)
    add_action(window, view_menu, "Center on &Molecule", "Home", window._center_molecule)

    measure_menu = window.menuBar().addMenu("&Measure")
    add_action(window, measure_menu, "&Pin Measurement", "M", window._pin_measurement)
    add_action(window, measure_menu, "&Clear Pinned Measurements", "Shift+M",
                     window._clear_pinned)

    spectra_menu = window.menuBar().addMenu("&Spectra")
    add_action(window, spectra_menu, "&Show/Hide Spectra", "S", window._toggle_spectra)
    add_action(window, spectra_menu, "Stop &Animation", "Space",
                     window.viewport.stop_animation)

    help_menu = window.menuBar().addMenu("&Help")
    add_action(window, help_menu, "&About", None, window.about)

def add_action(window, menu, text, shortcut, slot):
    action = QAction(text, window)
    if shortcut is not None:
        action.setShortcut(QKeySequence(shortcut))
    action.triggered.connect(slot)
    menu.addAction(action)
    return action
