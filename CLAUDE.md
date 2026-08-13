# CLAUDE.md — project guide

**microscp** ("Microscope") — cross-platform molecular structure & spectroscopy
viewer with CYLview-style rendering. Repo: https://github.com/lucienwb/microscp

## Environment & commands

```bash
conda activate microscp          # python 3.12 env; pip install -e ".[dev]" done
python -m pytest tests/ -q       # 55 tests, headless-safe (no GL/display needed)
microscp <file>                  # launch the GUI
python scripts/preview.py <file> out.png    # offscreen render (visual checks)
python scripts/demo_features.py  # GUI feature showcase
```

Absolute env python: `/opt/homebrew/Caskroom/miniforge/base/envs/microscp/bin/python`

## Hard constraints (agreed with maintainer)

- Dependencies limited to: numpy, scipy, matplotlib, pandas, PySide6, PyOpenGL.
  **cclib / Open Babel must stay optional** (fallback only, never imported at top level).
- All parsers are in-house; study cclib source for format knowledge if needed.
- Stays a pip-installable CLI/GUI app — no bundled installers.
- The maintainer runs `git commit`/`push` themselves — stage changes and
  provide the command lines instead of committing.

## Architecture

- `core/` — Molecule (+bond perception), geometry, editing (fragment-aware
  set_distance/angle/dihedral with ring fallback), contacts (H-bond detection),
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
- `spectra/` — broadening + ir/uvvis/nmr Spectrum builders (no Qt).
- `gui/` — viewport (picking via ID-buffer render, unlimited click-selection
  with measurement readout at 2–4 atoms, labels,
  snapshot-based undo/redo, QTimer vibration animation, volume/isosurface
  state with set_volume/set_isosurface incl. lobe colors), annotations.py
  (QPainter measurement/label drawing shared by viewport overlay and image
  export: atom labels centered on atoms, distance text parallel to bond,
  compact angle arcs with the value beside them, dihedral rotation arrows
  spanning between solid arms that overlay the outer bonds in the
  down-the-bond view, via geometry.angle_arc_points/dihedral_arc_points/
  dihedral_arm_points), spectra dock
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
  scripts/preview.py and scripts/demo_features.py with --screenshot).
- macOS GL is 4.1 core; shaders are `#version 330 core`; picking/sample-shading
  guards exist for GL < 4.0.

## Status & next steps

Done: M1 viewer, M2 formats (Gaussian/ORCA/Q-Chem/Molden native), M3 spectra
(IR + mode animation, UV-Vis, NMR), M4 editing (adjust/delete/undo), extras
(input writers, H-bond dashes, cube-file isosurface display with color
presets/custom colors, publication-style measurement annotations,
transparent PNG / uncompressed TIFF export, Houk/Houkmol style toggled with
V, per-region ball&stick/stick/line representations on keys 1/2/3). See
CHANGELOG.md.

Agreed next candidates: computing MO/density grids natively from fchk/Molden
(parse MO coefficients + basis, evaluate on a grid, feed the existing
isosurface pipeline — no cubegen needed), style editor, PyPI release.
