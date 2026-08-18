# microscope

**microscope** is a cross-platform molecular **structure and
spectroscopy viewer** for quantum chemistry, with CYLview-style
publication-quality rendering. It installs the short command `scope`.

It reads the everyday files of computational chemistry — Gaussian, ORCA and
Q-Chem outputs, `xyz`, `pdb`, `fchk`, Molden — and turns them into clean 3D
structures, interactive IR/UV-Vis/NMR spectra, editable geometries, and
figures ready for a paper.

![tryptophan rendered by microscope](docs/screenshot_trp.png)

## Highlights

- **CYLview-look rendering** — ray-cast impostor spheres/cylinders (never
  faceted), split-color bonds, orthographic camera, high-resolution export
  to transparent PNG or uncompressed TIFF
- **Houk (Houkmol) style** — one key (`V`) switches to the classic Houk-group
  look: glossy ball-and-stick with black bonds, near-white carbons and the
  signature quadrant seam lines on every heavy atom, readable even in
  black-and-white printouts
- **Mixed representations** — select a region and give it its own level of
  detail: ball-and-stick for the important part (`1`), sticks for the
  surroundings (`2`), thin lines for the rest (`3`)
- **One viewer for all programs** — Gaussian / ORCA / Q-Chem outputs are
  recognized by content and parsed natively (no cclib dependency)
- **Interactive spectra** — click an IR band and *watch the molecule vibrate*;
  UV-Vis and NMR with adjustable broadening, scaling and referencing
- **Orbital & density isosurfaces** — open a cube file and the surface appears,
  with live isovalue/opacity/color controls (in-house marching-tetrahedra
  mesher; preset color pairs or any custom RGB/HEX color)
- **Drag to move and rotate** — select an atom or a whole fragment and a
  manipulator appears on it: coloured X/Y/Z arrows slide it, three rings turn
  it, the centre dot moves it in the screen plane (hold Shift to snap to
  0.1 Å / 15°). Everything is undoable
- **Measure & edit** — publication-style measurement annotations (distance
  written along the bond, angles marked with an arc, dihedrals with a
  rotation arrow around the central bond); fragment-aware geometry
  adjustment with live preview and undo/redo
- **Trajectories** — optimization/IRC/scan playback with energies
- **Batch figures from the shell** — `scope -s mycalc.log --style houk
  --view 30,-15 -o fig.png` renders without opening a window, so figures
  regenerate from a script like a gnuplot plot
- **Scriptable** — the same parsers and writers work headless from Python

## Install

Runs anywhere Python 3.10+ runs (Windows, macOS, Linux):

```bash
git clone https://github.com/lucienwb/microscope
cd microscope
conda create -n microscope python=3.12   # or any venv
conda activate microscope
pip install -e .

scope mycalc.log                         # done
```

Dependencies are deliberately minimal and pip-installable everywhere:
`numpy`, `scipy`, `matplotlib`, `pandas`, `PySide6`, `PyOpenGL`. All file
parsers are written in-house; if `cclib` happens to be installed it is used
automatically as a fallback for formats the native parsers do not cover.

## Usage

```bash
scope                 # open the GUI
scope mycalc.log      # open a file directly
scope -s mycalc.log   # no window: render mycalc.png and exit
```

### Making figures from the command line

`-s`/`--silent` turns Microscope into a batch plotter, the way `gnuplot` makes
a plot from a script: no window opens, the figure is written and the command
returns. Put the line in a shell script or a Makefile and every figure
regenerates itself after a re-optimization.

```bash
scope -s opt.log -o fig.tif --style houk --size 2000x1500
scope -s mol.xyz --view top --labels number --measure 3,4
scope -s mol.log --rep '1-12:ball;13-40:line' --view 30,-15 --axes
scope -s homo.cube --iso 0.03 --iso-colors purple,gold --zoom 1.3
scope -s scan.log --frame 12 --align 3,4 -o frame12.png
```

Output is a transparent PNG (or uncompressed TIFF) named after the input
unless `-o` says otherwise; `--opaque` fills the background instead. Atom
numbers are 1-based, exactly as the viewer labels them. `scope --help` lists
every flag — style, representations, labels, view/rotation/alignment, zoom,
size, supersampling, measurements, cube isosurfaces and colours.

`--ir`, `--uv` and `--nmr` plot the spectrum instead of the molecule, with the
same broadening and axes as the spectra dock — and vector output for journals:

```bash
scope -s freq.log --ir -o ir.pdf --fwhm 12 --freq-scale 0.965
scope -s freq.log --ir --xrange 600:1800 --title 'fingerprint region'
scope -s td.log   --uv --unit eV --xrange 2:7 -o uv.svg
scope -s nmr.log  --nmr --nucleus C --reference 186.4 --csv shifts.csv
```

IR and δ axes run high → low as spectra are conventionally drawn, sticks sit
under the broadened curve (`--no-sticks` drops them), UV-Vis oscillator
strengths get their own right-hand axis, and `--csv` writes the plotted curve
as data. Figures go to PNG, PDF, SVG, EPS or TIFF at `--dpi` (300 by default).

`--style`, `--labels` and `--axes` also work without `-s`, to open the viewer
already set up that way; the render-only flags refuse to run silently ignored.

- **Rotate**: left-drag &nbsp;·&nbsp; **Pan**: right-drag &nbsp;·&nbsp; **Zoom**: scroll
- **Representation**: `V` toggles between the CYLview look and the Houk
  (Houkmol) ball-and-stick style (View → Representation).
  Select atoms and press `1` (ball & stick), `2` (stick) or `3` (line) to
  mix levels of detail — e.g. active site in full detail, spectator ligands
  as sticks, everything else as lines; with nothing selected the whole
  molecule switches
- **Select & measure**: click atoms to build a selection of any size (click
  again to deselect one; empty space or Esc clears). While 2, 3 or 4 atoms
  are selected the distance / angle / dihedral is shown live — distances
  written parallel to the bond, angles with an arc at the vertex, dihedrals
  with a rotation arrow around the central bond (Newman-style); larger
  selections just show the atom count, ready for region operations like
  `1`/`2`/`3`
- **Move & rotate by hand**: a manipulator sits on whatever is selected —
  drag the red/green/blue **arrows** to slide it along x/y/z, the matching
  **rings** to turn it about that axis, or the **centre dot** to move it in
  the screen plane; hold `Shift` to snap to 0.1 Å / 15°. `F` first grows the
  selection to the whole connected fragment, so a substituent or ligand is
  picked up in one keystroke. `G` hides the handles when they are in the way

![the move/rotate manipulator and the XYZ axis indicator](docs/screenshot_gizmo.png)

- **XYZ axes**: `Shift+A` shows a corner triad of the world axes (View → Show
  XYZ Axes); it is included in exported images too
- **Pin measurements**: `M` keeps the current measurement displayed in the scene (several at once; `Shift+M` clears)
- **Atom labels**: `L` cycles element / element+number / number (View → Atom Labels)
- **Align view**: `A` with 2 atoms selected looks down the bond; with 3 atoms selected puts their plane in the screen (GaussView-style)
- **Rotation center**: `C` rotates about the selected atom; `Home` re-centers on the molecule
- **Edit geometry**: select 2–4 atoms and press `E` — a dialog adjusts the
  distance/angle/dihedral **with live preview**; the attached fragment moves
  along (GaussView-style; ring bonds move only the end atom)
- **Delete atoms**: select and press `X` &nbsp;·&nbsp; **Undo/Redo**: `Ctrl+Z` / `Ctrl+Shift+Z`
- **Recompute bonds**: `B` (after large geometry changes)
- **Hydrogen bonds**: shown automatically as CYLview-style dashed lines (D–H···A
  with H···A ≤ 2.6 Å and angle ≥ 120°); `H` toggles them
- **Isosurfaces**: open a `.cube` file (from `cubegen`, `orca_plot`, …) and the
  orbital/density surface shows automatically; `I` opens live controls
  (isovalue, opacity, colors — curated presets or any RGB/HEX color per lobe —
  and MO selection for multi-MO cubes)

![orbital isosurface from a cube file](docs/screenshot_isosurface.png)
- **Export**: File → Export Image… for high-resolution figures — transparent
  PNG or uncompressed TIFF — with isosurfaces, atom labels, pinned
  measurements and the XYZ axes included
- **Save**: File → Save As… to write `xyz`, Gaussian input (`.gjf`), or `pdb` — including edited structures

![labels, pinned measurements and plane alignment](docs/screenshot_features.png)

![the same molecule in the Houk (Houkmol) style](docs/screenshot_houk.png)

## Supported files

| Format | Read | Write |
|---|---|---|
| XYZ (single & multi-frame) | ✅ | ✅ |
| Gaussian output (`.log`/`.out`) | ✅ geometries, SCF energies, frequencies + IR + normal modes, TD-DFT, NMR shieldings | — |
| Gaussian input (`.gjf`/`.com`) | ✅ | ✅ |
| Gaussian formatted checkpoint (`.fchk`) | ✅ | — |
| ORCA output (`.out`) | ✅ geometries, energies, frequencies + IR + normal modes, TD-DFT absorption, NMR shieldings | — |
| ORCA input (`.inp`) | — | ✅ |
| Q-Chem output (`.out`) | ✅ geometries, energies, frequencies + IR + normal modes, TD-DFT | — |
| Q-Chem input (`.in`/`.qcin`) | — | ✅ |
| Molden (`.molden`) | ✅ geometries, frequencies + normal modes | — |
| Cube (`.cube`/`.cub`) | ✅ geometry + volumetric MO/density grids (incl. multi-MO) | — |
| PDB | ✅ | ✅ |

Program outputs are recognized by content, not extension — any `.log`/`.out`
file is sniffed and routed to the right parser. If `cclib` is installed it
serves as a fallback for anything the native parsers miss.

## Spectra

Open any frequency / TD-DFT / NMR job and the spectra panel appears
automatically (`S` toggles it):

- **IR** — Lorentzian-broadened spectrum with sticks; adjustable FWHM and
  frequency scaling factor. **Click a band and the 3D structure animates that
  vibrational mode** (click again or press Space to stop).
- **UV-Vis** — Gaussian broadening in the energy domain; nm or eV axis;
  oscillator-strength sticks on a secondary axis.
- **NMR** — per-nucleus spectrum (¹H, ¹³C, …); enter a reference shielding σ₀
  to switch from raw shieldings to the chemical-shift (δ) scale.

Every spectrum exports as CSV (data) or PNG/SVG/PDF (figure) for publications.

![IR spectrum with the animated C–H stretch selected](docs/screenshot_ir_spectrum.png)

## Roadmap

- **M1 — viewer** ✅ CYLview-style rendering, measurements, labels, image export
- **M2 — formats** ✅ native ORCA/Q-Chem/Molden parsers, trajectory playback
- **M3 — spectra** ✅ IR with mode animation, UV-Vis, NMR
- **M4 — editing** ✅ distance/angle/dihedral adjustment with live preview, atom deletion, undo/redo
- **Extras** ✅ ORCA/Q-Chem input writers, dashed hydrogen-bond display,
  cube-file orbital/density isosurfaces
- **Ideas for later** — computing MO/density grids natively from fchk/Molden
  (no cubegen needed), style editor, PyPI release

## Python API

Everything the GUI does with files is scriptable:

```python
import microscope.io as mio

result = mio.load("mycalc.log")            # any supported format
mol = result.molecule                       # final geometry
print(mol.formula(), result.scf_energies[-1])
for v in result.vibrations[:3]:
    print(v.frequency, v.ir_intensity)
mio.save_molecule("next_job.gjf", mol)     # ready to resubmit
```

See `examples/api_demo.py` for a full tour (trajectories, IR, TD-DFT, NMR,
format conversion).

## Development

```bash
pip install -e ".[dev]"
pytest                       # 55 tests against real QM output files
python scripts/preview.py    # offscreen render smoke test -> preview.png
python scripts/demo_features.py   # GUI feature showcase
```

Tests run headless (no display needed) and cover the parsers, geometry
engine, editing operations, and spectra models on all three platforms via
GitHub Actions.

## Acknowledgments

- The rendering style is inspired by **CYLview** (C. Y. Legault, Université
  de Sherbrooke) — this project is an independent implementation.
- Parser test fixtures in `tests/data/` come from the **cclib** project's
  test suite (BSD-3-Clause); see `tests/data/README.md`.

## License

MIT
