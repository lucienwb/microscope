# Changelog

## Unreleased

- Spectra from the command line: `scope -s file.log --ir/--uv/--nmr` plots the
  spectrum instead of the molecule, to PNG/PDF/SVG/EPS/TIFF at `--dpi`, with
  `--fwhm`, `--freq-scale`, `--unit nm|eV`, `--nucleus`, `--reference`,
  `--xrange`, `--no-sticks`, `--title` and `--csv` (writes the plotted curve
  as data). The drawing code is shared with the spectra dock, so batch figures
  match the screen; the UV-Vis oscillator-strength label now stays on the
  right-hand axis after a redraw
- Silent mode: `scope -s file.log …` renders a figure straight to PNG/TIFF
  without opening a window, so figures can be regenerated from a shell script
  or Makefile the way gnuplot makes a plot. Flags cover style, per-range
  representations, labels, view (named orientations, screen-axis angles,
  `--align` down a bond or plane), zoom, size, supersampling, background and
  bond colours, measurements, the XYZ axes, cube isosurfaces with level,
  opacity and lobe colours, and trajectory frame selection. Atom numbers are
  1-based as in the viewer, and `--style`/`--labels`/`--axes` also work
  without `-s` to open the viewer already configured
- Interactive move/rotate manipulator on the selection (LightCone-style):
  coloured X/Y/Z arrows translate it, three rings rotate it about the
  corresponding axis, and a centre dot moves it in the screen plane. Handles
  light up on hover, `Shift` snaps to 0.1 Å / 15°, every drag is a single
  undo step, and `G` (Edit → Show Move/Rotate Handles) turns them off
- `F` grows the selection to the whole connected fragment, so a substituent
  or ligand can be grabbed and moved in one go
- `Shift+A` (View → Show XYZ Axes) shows a corner triad of the world axes,
  faded for axes pointing away from the viewer; exported images include it
- **Renamed**: the project, distribution and import package are now
  `microscope` (was `microscp`), and the console command is now the shorter
  `scope` (was `microscp`). Update existing installs with a fresh
  `pip install -e .`; `import microscp` becomes `import microscope`
- Houk (Houkmol) style: the classic Houk-group figure look — glossy
  ball-and-stick with black bonds, near-white carbons and two black
  great-circle "quadrant" seam lines on every heavy atom (they rotate with
  the molecule), readable even in black-and-white print. `V` (or View →
  Representation) toggles between the CYLview and Houk looks; exports and
  `scripts/preview.py --style houk` use the same styles
- Mixed representations: select a region and press `1`/`2`/`3` to draw it
  as ball-and-stick / sticks / thin lines (nothing selected = whole
  molecule) — highlight the chemistry in full detail and keep the
  environment lightweight; bonds between regions taper to the thinner style
- Click-selection is no longer capped at 4 atoms: keep clicking to select
  any number (for region operations); the live distance/angle/dihedral
  readout still appears whenever exactly 2–4 atoms are selected

- Publication-style measurement annotations: distances are written parallel
  to the bond and centered on it, angles are marked with a compact arc in the
  plane of the three atoms with the value beside it, dihedrals with a
  rotation arrow around the central bond spanning between the two outer
  bonds (Newman-style) — in the viewport and in exported images
- Atom labels are drawn centered on their atoms
- Isosurface color controls in the `I` dialog: curated preset pairs
  (blue/red, teal/orange, …) plus a full color picker (RGB/HEX) per lobe
- Export Image now writes PNG or uncompressed TIFF (both keep the
  transparent background) and includes the displayed isosurface, atom
  labels and pinned measurements (toggleable)
- Cube file (`.cube`/`.cub`) reader: geometry plus volumetric MO/density
  grids, including cubegen-style multi-MO cubes
- Orbital/density isosurfaces rendered as smooth translucent meshes
  (in-house marching-tetrahedra extraction, gradient normals, depth-resolved
  transparency); surfaces appear automatically when a cube is opened and
  `I` opens live isovalue/opacity/MO controls; also available headless via
  `render_molecule_image(..., volume=...)` and `scripts/preview.py`
- ORCA (`.inp`) and Q-Chem (`.in`/`.qcin`) input writers, available from
  File → Save As and `microscope.io.save_molecule`
- Hydrogen bonds drawn as CYLview-style dashed lines (D–H···A criteria:
  H···A ≤ 2.6 Å, angle ≥ 120°); toggle with `H`
- Save dialog appends the correct extension when none is typed
- Test suite grown with more real program outputs (from the cclib
  collection): Gaussian 09 unrestricted SP, relaxed PES scan and a
  Mo₄OCl₄²⁻ transition-metal single point, ORCA 5.0 optimization,
  Q-Chem 6.0 SMD solvent job, unrestricted `.fchk`
- Fixed a missing Qt import that crashed the Adjust dialog (`E`)

## 0.1.0 — 2026-08-12

Initial release.

- **Viewer**: CYLview-style OpenGL renderer (ray-cast impostor spheres and
  cylinders, split-color bonds, orthographic camera), atom picking,
  distance/angle/dihedral measurements with pinnable in-scene labels,
  zoom-adaptive atom labels, view alignment (along bond / plane to screen),
  rotation centering on an atom, high-resolution transparent PNG export
- **File formats**: native parsers for Gaussian (`.log`/`.out`, `.gjf`,
  `.fchk`), ORCA, Q-Chem, Molden, XYZ (multi-frame) and PDB; content-based
  program detection; XYZ/GJF/PDB writers; optional cclib fallback
- **Spectra**: IR (Lorentzian broadening, frequency scaling, click a band to
  animate the vibrational mode in 3D), UV-Vis (Gaussian broadening, nm/eV),
  NMR (per-nucleus, σ→δ referencing); CSV and PNG/SVG/PDF export
- **Editing**: fragment-aware distance/angle/dihedral adjustment with live
  preview, atom deletion, bond re-perception, 100-level undo/redo
- **Trajectories**: frame slider with SCF energies and playback button
