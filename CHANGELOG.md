# Changelog

## Unreleased

- ruff and mypy are configured and both pass. They immediately found three
  things the tests could not: a `self` left behind in a method extracted to a
  module-level function, a re-exported constant an autofix had removed (which
  broke every GUI import), and a module-level global dropped in a refactor
  that would have failed on the first `scope -s` render
- `ParseResult.extras` is gone; the route line and which parser read the file
  are named fields
- Undo/redo is `core.history` — Snapshot and EditHistory, plain data with no
  Qt, tested on its own. What a selection measures and how it is worded is
  `core.measure` for the same reason
- The menu bar is `gui/menus.py`, one table of what the program can do and
  which key does it. Adding a command is one entry there and one method on the
  window
- Every module is now import-tested, and the Lewis tests are split by what
  they cover (perception, layout, files) with the shared geometries in
  tests/molecules.py

- `microscope` now has a public API: `load`, `save_molecule`, `save_lewis`,
  `Molecule`, `ParseResult`, `perceive` and the errors worth catching, listed
  in `__all__`. Importing it pulls in no Qt, so a script that only wants the
  numbers never pays for graphics
- The two largest modules are split along what they do. `gui/app.py` (922
  lines, a window class with 42 methods plus three dialogs plus the entry
  point) is now app.py + mainwindow.py + dialogs.py + filetypes.py, and
  `cli.py` (741 lines) is a package: parsing, options, and one module per
  thing `-s` can produce
- `scope -s file.log --ir` no longer needs the Qt runtime libraries it never
  used, and writing a `.cdxml` loads no graphics library at all

- Result arrays are stored as narrowly as the data allows. The adjacency was
  a Python list per atom holding boxed integers — 8 MB for photosystem II,
  more than every array in the structure put together — and is now the sparse
  (CSR) form of the bond matrix in two int32 arrays, 0.65 MB, indexed exactly
  as before. Bonds are int32, atomic numbers int16, bond orders and formal
  charges and lone pairs and radicals and hydrogen counts int8; coordinates
  stay float64, since geometry is where precision matters. 13 MB down to 2.8
- Cube grids and isosurface meshes are float32: a cube file carries five
  significant digits and single precision holds seven. With int32 grid indices
  and the gradient taken one axis at a time, extracting an isosurface from a
  128^3 grid peaks at 80 MB instead of 150, with vertices still landing within
  3e-3 A of an analytic sphere
- Line broadening accumulates in blocks instead of building the whole
  (grid x peaks) matrix, which for a few thousand normal modes was hundreds of
  megabytes for a sum that never needed it. Results are bit-identical

- Ions and radicals are drawn instead of merely flagged. The charge a file
  declares was only ever compared against, never used: it now goes on the one
  atom it can belong to when there is exactly one — a bare chloride, an ion in
  a fragment calculation — and otherwise on **square brackets round the whole
  drawing**, with the charge and a dot per unpaired electron outside, which is
  how a delocalized radical cation is drawn anyway. `[divinylbenzene]•+` comes
  out right rather than as a neutral molecule with a warning
- Formats that never state a charge (xyz, pdb, cube, molden) say so
  (`Molecule.charge_known`) instead of defaulting to zero and being contra-
  dicted by every ion in the file
- A drawing containing a metal says plainly that the charges around it are
  bookkeeping — a dative bond drawn as a plain line has to put a charge on the
  donor, and the balancing one belongs on the metal, which is deliberately
  left out of the counting
- Q-Chem ghost atoms (counterpoise centres, written `GH`) are no longer read
  as real atoms, and a fragment job's per-fragment geometry blocks no longer
  win over the whole system — `CH3---Na+` was arriving as one sodium atom

- **Big structures now open.** Bond perception compared every pair of atoms,
  so a protein needed an N^2 difference array — hundreds of gigabytes for a
  crystal structure, and the process was killed. It walks a grid of cells one
  bond wide instead, finding exactly the same bonds: photosystem II (54k
  atoms) goes from out-of-memory to 0.1 s, and the whole draw — parse,
  perceive, lay out — is under 3 s for every entry tried, up to 98k atoms.
  The Lewis overlap count and `elements.normalize_symbol` had the same shape
  of problem and got the same treatment
- A file whose carbons are mostly short of four bonds has had its hydrogens
  left out, the way an X-ray structure does. Formal charges are no longer
  invented for it — a protein was drawing as tens of thousands of carbanions
- Gaussian outputs written from a Z-matrix (older versions label the geometry
  block "Z-Matrix orientation") are read instead of refused

- **Lewis structure mode** (`Shift+L`, View → Lewis Structure): draws the
  molecule flat and skeletal the way ChemDraw does — bare carbon vertices,
  heteroatom labels with their hydrogens folded in (OH, NH₂, and HO– when the
  bond leaves to the right), double and triple bonds as parallel lines with
  the inner line inside the ring, formal charges, and lone-pair dots and
  radical dots on request. The drawing is the molecule seen through the same
  camera as the 3-D view, so **turning it chooses the angle the drawing is
  made from**: line the structure up, switch over, drag to turn (Shift+drag
  spins it in the plane of the page), then save. The mode is read-only — the
  Edit menu is greyed out — because picking the angle is all it is for
- Bond orders, formal charges, lone pairs and radicals are perceived from the
  geometry (reference bond lengths per element pair, then valence saturation
  handing out the most contracted bonds first, which is what makes a benzene
  ring come out as a Kekulé structure; where that still strands an atom, as
  it does on fused rings from tetracene on, an alternating chain of candidate
  bonds is flipped to free it). Perception cannot be right for
  everything — a carbocation and a carbanion have the same connectivity — so
  the perceived total charge is compared against the charge in the file and
  the status bar says so when they disagree, rather than quietly showing a
  wrong drawing. Atoms stacked behind each other in the current view are
  flagged the same way, since the fix is to turn the structure
- Save the drawing for ChemDraw (File → Save for ChemDraw…, `Ctrl+Shift+S`):
  **ChemDraw CDXML** or an **MDL molfile**, both written at the angle on
  screen with the perceived bond orders, charges and implicit hydrogens, so
  the structure arrives ready to rearrange rather than as a bag of atoms
- Export the drawing as line art: PNG and TIFF as before, plus **SVG and PDF**
  — a Lewis structure is vector art, and now stays that way in a manuscript
- `scope -s mol.log --lewis` does all of it headless: `-o` takes an image
  (`.png`/`.tif`), vector art (`.svg`/`.pdf`) or a structure file
  (`.cdxml`/`.mol`), with `--show-hydrogens`, `--carbon-labels`,
  `--lone-pairs` and `--color-atoms`. Writing a `.cdxml`/`.mol` needs no
  graphics at all. `--lewis` without `-s` opens the viewer in that mode
- `--spin DEG` turns the picture in the plane of the page, for the command
  line and (as Shift+drag) in the viewer
- File → Save As… also writes MDL molfiles (`.mol`/`.sdf`)

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
