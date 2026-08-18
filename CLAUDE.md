# CLAUDE.md — project guide

**microscope** — cross-platform molecular structure & spectroscopy
viewer with CYLview-style rendering. Package/import name `microscope`, CLI
command `scope`. Repo: https://github.com/lucienwb/microscope
(renamed from `microscp` on 2026-08-18)

## Environment & commands

```bash
conda activate microscope          # python 3.12 env; pip install -e ".[dev]" done
python -m pytest tests/ -q       # 112 tests, headless-safe (no GL/display needed)
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

- `core/` — Molecule (+bond perception), geometry, editing (fragment-aware
  set_distance/angle/dihedral with ring fallback, plus rigid translate_atoms/
  rotate_atoms and connected_fragment used by the manipulator),
  contacts (H-bond detection),
  volume.py (VolumeData grids), isosurface.py (vectorized marching
  tetrahedra + gradient normals), results dataclasses, elements data.
  Pure numpy, no Qt.
- `io/` — one module per format (gaussian, orca, qchem, molden, fchk, xyz,
  pdbfile, cube) + writers (xyz/gjf/pdb/ORCA inp/Q-Chem in). `io.load()`
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
  renderer.set_style_params(quad_color, quad_width)), make_style(name).
- `cli.py` — the `scope` entry point (`[project.scripts]` → `microscope.cli:main`).
  Without `-s/--silent` it calls gui.app.main(path, style=, label_mode=, axes=);
  with it, it renders headless via render.offscreen + annotations and writes
  through render/imagefile.py (write_image, shared with the GUI export). Pure
  parsing helpers (parse_size/indices/color/reps/rotation, build_camera,
  pick_frame, NAMED_VIEWS) are separated from the GL work so tests stay
  headless. --ir/--uv/--nmr switch it to plotting a spectrum (matplotlib Agg,
  PNG/PDF/SVG/EPS/TIFF) via spectra/plot.py. Render-only flags are rejected in
  GUI mode instead of ignored, and structure vs spectrum flags are kept apart
  (check_mode_flags); on Linux without DISPLAY it sets QT_QPA_PLATFORM=offscreen.
- `spectra/` — broadening + ir/uvvis/nmr Spectrum builders, plus plot.py
  (draw_spectrum/style_axes/set_xrange on a caller-supplied matplotlib axes)
  shared by the GUI dock and the CLI. No Qt, no pyplot. Note: y limits are set
  explicitly from the data — pinning the bottom to 0 disables autoscale, so an
  inferred top clips tall peaks; and a twin axis needs tick_right()/
  set_label_position("right") re-applied after clear().
- `gui/` — viewport (picking via ID-buffer render, unlimited click-selection
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
  spectra dock
  (matplotlib tabs; IR click-band → animate mode), app.py (menus/shortcuts,
  SurfaceDialog on `I` with SURFACE_PALETTES presets + QColorDialog,
  Export Image → transparent PNG / uncompressed TIFF via _write_image
  (user wants transparency-capable formats only), exports include
  isosurface + annotations).

## Testing conventions

- Parser tests run against real program outputs in `tests/data/` (sourced from
  the cclib test suite — keep the attribution in tests/data/README.md).
- GUI/GL work is verified by offscreen rendering to PNG and *looking at the
  image* (scripts in the session scratchpad follow the pattern of
  scripts/preview.py and scripts/demo_features.py with --screenshot; since the
  CLI exists, `scope -s <file> …` is often the quickest way to get one).
- macOS GL is 4.1 core; shaders are `#version 330 core`; picking/sample-shading
  guards exist for GL < 4.0.

## Status & next steps

Done: M1 viewer, M2 formats (Gaussian/ORCA/Q-Chem/Molden native), M3 spectra
(IR + mode animation, UV-Vis, NMR), M4 editing (adjust/delete/undo), extras
(input writers, H-bond dashes, cube-file isosurface display with color
presets/custom colors, publication-style measurement annotations,
transparent PNG / uncompressed TIFF export, Houk/Houkmol style toggled with
V, per-region ball&stick/stick/line representations on keys 1/2/3,
LightCone-style move/rotate manipulator on the selection with `F` fragment
select and `G` toggle, `Shift+A` XYZ axis triad, gnuplot-style silent-mode
CLI `scope -s`, incl. --ir/--uv/--nmr spectrum plots). See CHANGELOG.md.

Agreed next candidates: style editor, PyPI release (both `microscope` and
`scope` are taken on PyPI, so the distribution name needs a decision).
Computing MO/density grids from basis sets was proposed and **rejected** by the
maintainer — see the pure-visualization constraint above.
