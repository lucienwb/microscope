# CLAUDE.md — project guide

**microscope** — cross-platform molecular structure & spectroscopy
viewer with CYLview-style rendering. Package/import name `microscope`, CLI
command `scope`. Repo: https://github.com/lucienwb/microscope
(renamed from `microscp` on 2026-08-18)

## Environment & commands

```bash
conda activate microscope          # python 3.12 env; pip install -e ".[dev]" done
python -m pytest tests/ -q       # 273 tests, headless-safe (no GL/display needed)
ruff check microscope tests scripts    # lint; mypy   # types, gui excluded
scope <file>                     # launch the GUI
scope -s <file> -o fig.png ...   # silent mode: render a figure, no window
python scripts/preview.py <file> out.png    # offscreen render (visual checks)
python scripts/demo_features.py  # GUI feature showcase
```

Absolute env python: `/opt/homebrew/Caskroom/miniforge/base/envs/microscope/bin/python`

## Hard constraints (agreed with maintainer)

- Dependencies limited to: numpy, scipy, matplotlib, pandas, PySide6, PyOpenGL.
  **cclib / Open Babel must stay optional** (fallback only, never imported at top level).
- All parsers are in-house; study cclib source for format knowledge if needed.
- **Pure visualization tool**: it draws what QC programs produced and never
  computes chemistry itself. Evaluating MO/density grids from basis sets was
  explicitly rejected by the maintainer — cube files stay read-only input.
- Stays a pip-installable CLI/GUI app — no bundled installers.
- The maintainer runs `git commit`/`push` themselves — stage changes and
  provide the command lines instead of committing.

## Architecture

- Result arrays are stored as narrowly as the data allows: bonds int32,
  atomic numbers int16, bond orders / formal charges / lone pairs / radicals /
  hydrogen counts int8, cube grids and isosurface meshes float32 (coordinates
  stay float64 — geometry is where precision matters). LewisStructure.neighbors
  is an Adjacency: CSR-style offsets+flat int32 arrays, i.e. the sparse form of
  the symmetric bond matrix, indexed like the list of lists it replaced. On
  photosystem II that is 13 MB of result arrays down to 2.8.
- Perception is array work, not per-atom Python: bond lengths in one call,
  bonds grouped by element pair so the reference lengths are looked up once
  per pair rather than once per bond, valence/charge/lone-pair counts through
  arrays indexed by atomic number, and bincount wherever something is summed
  onto atoms. That is 4x on a protein. Two things were tried and rejected
  with numbers: scipy's cKDTree finds the same bonds 1.4x faster but costs
  more to import (0.25 s) than it saves, and visiting only the forward half
  of each cell's neighbourhood made bond perception 3.5x *slower*, because
  one large numpy block per cell beats fourteen small ones.
- `core/` — Molecule (bond perception through a grid of cells one bond wide,
  not an N^2 pairwise pass — a protein is tens of thousands of atoms and the
  difference array alone would be hundreds of gigabytes; elements.py memoizes
  normalize_symbol because Molecule.atomic_numbers rebuilds from the symbols
  on every read), lewis.py (bond orders/formal
  charges/lone pairs/radicals perceived from geometry: per-element-pair
  reference bond lengths, then valence saturation taking the most contracted
  bonds first — plain greedy strands two atoms of a six-ring and draws
  benzene with two double bonds, so ties go to the most constrained bond,
  and whatever greed still strands (fused rings from tetracene on) is
  rescued by flipping an alternating chain of candidate bonds (_augment).
  A bond that plainly wants a triple claims its valence before the doubles
  are handed out, but only where both atoms can afford every triple they
  want - that is what keeps CO2 two doubles while a nitrile next to a ring
  keeps its triple instead of leaving a nitrogen anion;
  LewisStructure.matches_file compares the perceived charge with the file's
  and the UI warns when they differ, unless Molecule.charge_known is False
  because the format never states one (xyz/pdb/cube/molden); a charge the
  file declares but the drawing has not placed is reconciled by _reconcile —
  onto the one unbonded atom if there is exactly one, else onto the species,
  where LewisStructure.net_charge/net_radicals put it on square brackets the
  way a delocalized radical cation is drawn. Metal complexes are left out of
  that: their mismatch comes from drawing dative bonds as plain lines, and
  the drawing says so instead. Around a metal the charges follow the ionic
  (oxidation-state) convention: a ligand is counted without its bond to the
  metal, so a phosphine keeps its lone pair and stays neutral instead of
  becoming a phosphonium, an atom bonded only to metals is that ligand as its
  own ion (Cl-, H-, oxo 2-), and whatever the ligands do not account for goes
  on the metal, where it reads as the oxidation state - Rh(+1) on a
  hydroformylation catalyst, Mn(+1) on a pincer carbonyl, Mo(+4) on MoOCl4 2-.
  It needs one metal and a stated charge to have somewhere unambiguous to put
  it; with two metals the drawing says it could not. A file whose carbons are
  mostly short of four
  bonds has had its hydrogens left out, as X-ray structures do, and then no
  formal charges are claimed at all — otherwise a protein draws as tens of
  thousands of carbanions), geometry, editing (fragment-aware
  set_distance/angle/dihedral with ring fallback, plus rigid translate_atoms/
  rotate_atoms and connected_fragment used by the manipulator),
  contacts (H-bond detection), history.py (Snapshot + EditHistory: undo/redo
  as plain data, no Qt), measure.py (what a 2-4 atom selection measures and
  how it is worded), vibration.py (ModeAnimation: a normal mode swinging
  about a geometry; the viewport keeps only the timer),
  volume.py (VolumeData grids), isosurface.py (vectorized marching
  tetrahedra + gradient normals), results dataclasses, elements data.
  Pure numpy, no Qt.
- `io/` — one module per format (gaussian, orca, qchem, molden, fchk, xyz,
  pdbfile, cube) + writers (xyz/gjf/pdb/ORCA inp/Q-Chem in/molfile) plus
  cdxml.py + molfile.py::write_2d behind `io.save_lewis()` (ChemDraw CDXML and
  MDL V2000 written at the *screen* positions, so a structure opens in
  ChemDraw at the angle it was composed at; hidden hydrogens become implicit
  counts). `io.load()`
  sniffs .log/.out content to pick the parser; cclib_bridge is the optional
  fallback. Cube reader handles multi-MO cubes (negative natom convention).
- `render/` — custom OpenGL renderer: **ray-cast impostor** spheres/cylinders
  (GLSL in shaders.py, orthographic camera only), plus a translucent triangle
  -mesh path for isosurfaces (depth pre-pass then blend, so only the nearest
  layer shades — no back-face facet noise). scene.py builds instance buffers
  (incl. dashed H-bond segments, per-sphere seam flags, uniform bond color
  when style.bond_color set, per-atom REP_BALL/REP_STICK/REP_LINE codes via
  build_scene(..., reps=) — stick shrinks atoms to bond_radius, line to
  LINE_RADIUS, inter-region bonds take the thinner radius) and
  build_surface_meshes(volume, iso, style); offscreen.py renders publication
  PNGs (accepts volume=, reps=); styles.py holds presets: CYLview default +
  houk_style() ("Houkmol", matched to a real CYLview render: glossy, black
  bonds, C 0.90 gray, black great-circle seam lines on heavy atoms drawn in
  the sphere shader from world-x/y plane normals via
  renderer.set_style_params(quad_color, quad_width)), make_style(name) which
  takes a preset name *or* a path: a Style is plain data, and save_style /
  load_style keep it as JSON with elements by symbol and colours as #rrggbb so
  a group can hand-edit and share one. A style file only needs the keys it
  changes, and a bad one raises StyleError naming the key rather than a
  traceback. Both front ends read the same file (`--style ours.json`).
  lewis2d.py is the flat ChemDraw-style layout — pure numpy, no Qt: projects
  through the same camera (so turning the molecule picks the drawing's angle),
  folds hydrogens into labels, trims bonds to the label box, offsets the
  second line of a double bond into the ring, places charges/lone-pair dots,
  and reports overlapping atoms. Sizes come from the *unprojected* median bond
  length so the lettering does not shrink as a bond swings edge-on.
  lewisdraw.py paints a layout with QPainter (real font metrics feed the
  layout's label boxes) and writes PNG/TIFF, SVG and PDF — no GL involved,
  it is line art.
- `cli/` — the `scope` entry point (`[project.scripts]` → `microscope.cli:main`),
  split by what each part needs: parsing.py turns text into values and needs no
  graphics, options.py declares the flags and their compatibility, and
  structure.py / spectrum.py / lewis.py are the three things `-s` can produce.
  `cli/__init__.py` resolves everything but main/build_parser/CliError through
  a PEP 562 `__getattr__`, and the renderers import their drawing libraries
  inside the functions, so `import microscope.cli` and `scope -s --lewis -o
  x.cdxml` pull in no Qt, no OpenGL and no matplotlib at all (there is a
  subprocess test for exactly that), and `--ir` loads only matplotlib.
  Without `-s/--silent` it calls gui.app.main(path, style=, label_mode=, axes=);
  with it, it renders headless via render.offscreen + annotations and writes
  through render/imagefile.py (write_image, shared with the GUI export). Pure
  parsing helpers (parse_size/indices/color/reps/rotation, build_camera,
  pick_frame, NAMED_VIEWS) are separated from the GL work so tests stay
  headless. --ir/--uv/--nmr switch it to plotting a spectrum (matplotlib Agg,
  PNG/PDF/SVG/EPS/TIFF) via spectra/plot.py; --lewis switches to the flat
  drawing, where -o may be an image, .svg/.pdf, or a .cdxml/.mol structure
  file (that last path deliberately never starts Qt).
  Render-only flags are rejected in
  GUI mode instead of ignored, and structure vs spectrum flags are kept apart
  (check_mode_flags); on Linux without DISPLAY it sets QT_QPA_PLATFORM=offscreen.
- `spectra/` — broadening + ir/uvvis/nmr Spectrum builders, plus plot.py
  (draw_spectrum/style_axes/set_xrange on a caller-supplied matplotlib axes)
  shared by the GUI dock and the CLI. No Qt, no pyplot. Note: y limits are set
  explicitly from the data — pinning the bottom to 0 disables autoscale, so an
  inferred top clips tall peaks; and a twin axis needs tick_right()/
  set_label_position("right") re-applied after clear().
- `gui/` — app.py is only `main()`; the window is mainwindow.py, menus.py is
  the menu bar as one table (a new command is one entry there and one method
  on the window), the modal dialogs are dialogs.py, and filetypes.py holds the
  file-dialog filters with the extension each one implies (one place to add a
  format). viewport (picking via ID-buffer render, unlimited click-selection
  with measurement readout at 2–4 atoms, labels,
  snapshot-based undo/redo, QTimer vibration animation, volume/isosurface
  state with set_volume/set_isosurface incl. lobe colors), annotations.py
  (QPainter measurement/label drawing shared by viewport overlay and image
  export: atom labels centered on atoms, distance text parallel to bond,
  compact angle arcs with the value beside them, dihedral rotation arrows
  spanning between solid arms that overlay the outer bonds in the
  down-the-bond view, via geometry.angle_arc_points/dihedral_arc_points/
  dihedral_arm_points; draw_axis_indicator = toggleable corner XYZ triad,
  also exported), gizmo.py (LightCone-style manipulator on the selection:
  X/Y/Z arrows + rotation rings + centre dot, all screen-space hit-testing
  and QPainter drawing on top of the orthographic projection — arrows win
  over rings, edge-on axes below MIN_AXIS_PX are skipped; viewport owns the
  drag state and re-derives coords from a base snapshot each move),
  lewisview.py (read-only 2D page sharing the viewport's camera object, so
  both views stay on one orientation; drag turns, Shift+drag spins in the
  page, and the Edit menu is disabled while it shows),
  spectra dock
  (matplotlib tabs; IR click-band → animate mode), mainwindow.py (menus/shortcuts,
  SurfaceDialog on `I` with SURFACE_PALETTES presets + QColorDialog,
  Export Image → transparent PNG / uncompressed TIFF via _write_image
  (user wants transparency-capable formats only), exports include
  isosurface + annotations).

## Tooling

`ruff` and `mypy` are configured in pyproject.toml and must both pass.
Ruff's `F` rules are the ones that earn their keep — a name used but never
imported is exactly what slips through when code moves between modules, and
it has already caught a `self` left behind in an extracted method and a
re-exported constant an autofix removed. mypy skips `gui/` (PySide6's stubs
make `self.style` resolve to `QWidget.style()`, so a viewport with its own
style attribute draws dozens of identical false complaints) and checks the
rest clean; it found a module-level global dropped during a refactor that the
tests had not.

## Testing conventions

- tests/test_imports.py imports every module in the package. The suite drives
  the viewer through offscreen widgets but imports little of `gui`, so a bad
  import there used to survive a green run.
- tests/molecules.py holds the idealized geometries the Lewis tests are built
  from; the Lewis tests themselves are split by what they cover
  (perception / layout / files).
- Parser tests run against real program outputs in `tests/data/` (sourced from
  the cclib test suite — keep the attribution in tests/data/README.md).
- GUI/GL work is verified by offscreen rendering to PNG and *looking at the
  image* (scripts in the session scratchpad follow the pattern of
  scripts/preview.py and scripts/demo_features.py with --screenshot; since the
  CLI exists, `scope -s <file> …` is often the quickest way to get one).
- macOS GL is 4.1 core; shaders are `#version 330 core`; picking/sample-shading
  guards exist for GL < 4.0.
- Lewis perception is checked against textbook answers (benzene Kekule, fused
  rings kekulized all the way along, CO2 vs a triple bond, nitrate's charge
  separation, hypervalent sulfate, CO's C-/O+, a doublet becoming a radical
  rather than an anion, a trigonal boronate ester staying neutral while a
  four-coordinate borate does not, a metal carbonyl keeping its triple bond,
  a hydrogen-less file claiming no charges while nitrate keeps its, a bare
  chloride taking the file's -1 while a delocalized radical cation goes in
  brackets, Q-Chem ghost atoms and fragment-only geometry blocks skipped)
  — those are the regressions that matter, not the pixels. The reference
  lengths in lewis.py::_REFERENCE were checked against real catalyst
  geometries (examples/catalysis — the maintainer's own Mn pincer and Rh
  hydroformylation runs, gitignored on purpose, not in the repo): B-O is
  single at every length
  that occurs, and the C-O triple is 1.15 rather than 1.128 so that a metal
  carbonyl stretched by back-donation still reads as C(triple)O.
- CI (.github/workflows/tests.yml): tests import PySide6 since test_cli.py and
  test_gizmo.py exist, and on Linux runners `import PySide6.QtGui` fails until
  the Qt runtime libs are apt-installed (libegl1 libgl1 libx11-6 libxkbcommon0
  libfontconfig1 libfreetype6 libdbus-1-3 — the DT_NEEDED list of libQt6Gui.so.6
  minus base-image libs). QT_QPA_PLATFORM=offscreen is set for the whole
  workflow. No test may construct a QGuiApplication/QApplication.

## Status & next steps

Done: M1 viewer, M2 formats (Gaussian/ORCA/Q-Chem/Molden native), M3 spectra
(IR + mode animation, UV-Vis, NMR), M4 editing (adjust/delete/undo), extras
(input writers, H-bond dashes, cube-file isosurface display with color
presets/custom colors, publication-style measurement annotations,
transparent PNG / uncompressed TIFF export, Houk/Houkmol style toggled with
V, per-region ball&stick/stick/line representations on keys 1/2/3,
flat ChemDraw-style Lewis mode on Shift+L with CDXML/molfile/SVG/PDF output,
style editor on Ctrl+T with styles savable as shareable JSON files,
LightCone-style move/rotate manipulator on the selection with `F` fragment
select and `G` toggle, `Shift+A` XYZ axis triad, gnuplot-style silent-mode
CLI `scope -s`, incl. --ir/--uv/--nmr spectrum plots). See CHANGELOG.md.

Agreed next candidates: PyPI release (both `microscope` and
`scope` are taken on PyPI, so the distribution name needs a decision).
Computing MO/density grids from basis sets was proposed and **rejected** by the
maintainer — see the pure-visualization constraint above.
