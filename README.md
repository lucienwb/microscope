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
- **Lewis structure mode** — `Shift+L` redraws the molecule flat and
  skeletal the way ChemDraw does, with double bonds, charges and (on request)
  lone pairs. Turning it picks the angle the drawing is made from; save the
  picture, or save a **ChemDraw `.cdxml`/`.mol`** to keep editing there
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
- **Whole proteins, not just molecules** — bonds are found through a grid of
  cells rather than by comparing every pair, so a 98,000-atom PDB entry parses
  and draws in a few seconds instead of exhausting memory
- **Scriptable** — `import microscope` gives you the parsers and writers with
  no graphics attached, so a script that only wants numbers never starts Qt

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
scope                 # open the viewer
scope mycalc.log      # open a file in it
scope -s mycalc.log   # no window: render mycalc.png and exit
```

![turning a Lewis structure to pick the angle it is drawn from](docs/lewis_spin.gif)

*Lewis mode: the drawing is the molecule seen through the same camera as the
3-D view, so turning it is how you choose the angle it is drawn from.*

## Keyboard reference

Everything the viewer does, in one place. Every entry is also in the menus.

### Looking

| Key | What it does |
|---|---|
| left-drag / right-drag / scroll | rotate · pan · zoom |
| `V` | CYLview look ⇄ Houk (Houkmol) ball-and-stick |
| `1` `2` `3` | ball & stick · stick · line, for the selection (or all of it) |
| `L` | cycle atom labels: element → element+number → number |
| `Shift+A` | corner XYZ axis triad |
| `A` | align the view: down the bond (2 atoms) or into the plane (3 atoms) |
| `C` · `Home` · `Ctrl+R` | rotate about the selected atom · re-centre on the molecule · reset the view |
| `H` | hydrogen bonds (shown by default as dashed lines) |
| `S` · `Space` | show/hide the spectra panel · stop a vibration animation |

### Selecting, measuring, editing

| Key | What it does |
|---|---|
| click · click again · `Esc` | add an atom to the selection · remove it · clear |
| *(2–4 atoms selected)* | the distance / angle / dihedral is shown live |
| `M` · `Shift+M` | pin the measurement into the scene · clear the pinned ones |
| `F` | grow the selection to the whole connected fragment |
| `G` | show/hide the move-rotate handles |
| `E` | adjust the selected distance/angle/dihedral, with live preview |
| `X` | delete the selected atoms |
| `B` | recompute bonds after a large geometry change |
| `Ctrl+Z` · `Ctrl+Y` | undo · redo (the platform's standard keys) |

### Drawing and saving

| Key | What it does |
|---|---|
| `Shift+L` | Lewis structure mode (read-only; drag turns, `Shift`+drag spins) |
| `I` | isosurface controls: isovalue, opacity, colours, which MO |
| `Ctrl+E` | Export Image — PNG/TIFF, and SVG/PDF in Lewis mode |
| `Ctrl+Shift+S` | Save for ChemDraw — `.cdxml` or `.mol` at the angle on screen |
| `Ctrl+O` · `Ctrl+S` | open · save as (`xyz`, `.gjf`, `.inp`, `.in`, `pdb`, `.mol`) |

## The command line

`-s`/`--silent` turns microscope into a batch plotter, the way `gnuplot` makes
a plot from a script: no window opens, the figure is written, the command
returns. Put the line in a Makefile and every figure regenerates itself after
a re-optimization.

```bash
scope -s opt.log -o fig.tif --style houk --size 2000x1500
scope -s mol.xyz --view top --labels number --measure 3,4
scope -s mol.log --rep '1-12:ball;13-40:line' --view 30,-15 --axes
scope -s homo.cube --iso 0.03 --iso-colors purple,gold --zoom 1.3
scope -s scan.log --frame 12 --align 3,4 -o frame12.png
```

`--ir`, `--uv` and `--nmr` plot the spectrum instead of the molecule, with the
same broadening and axes as the spectra dock, and vector output for journals:

```bash
scope -s freq.log --ir -o ir.pdf --fwhm 12 --freq-scale 0.965
scope -s freq.log --ir --xrange 600:1800 --title 'fingerprint region'
scope -s td.log   --uv --unit eV --xrange 2:7 -o uv.svg
scope -s nmr.log  --nmr --nucleus C --reference 186.4 --csv shifts.csv
```

`--lewis` draws the flat structure, and `-o` decides what comes out — a
picture, vector art, or a structure file:

```bash
scope -s mol.log --lewis --align 2,3,4 -o scheme.svg
scope -s mol.log --lewis -o mol.cdxml        # opens in ChemDraw at that angle
scope -s mol.log --lewis --lone-pairs --color-atoms -o lewis.png
```

Writing a `.cdxml` or `.mol` needs no graphics library at all, and `--ir`
loads only matplotlib — neither starts Qt.

### Flags worth knowing

| Flag | |
|---|---|
| `-o FILE` | output; the extension picks the format. Default: the input's name |
| `--style cylview\|houk` · `--rep SPEC` | look, and per-region detail (`'1-12:ball;13-40:line'`) |
| `--view NAME` · `--rotate X,Y` · `--align A,B[,C]` | named view · turn by degrees · look down a bond or into a plane |
| `--zoom` · `--size WxH` · `--supersample N` | framing and resolution |
| `--labels MODE` · `--measure A,B[,C[,D]]` · `--axes` | annotations drawn into the figure |
| `--cube FILE` · `--iso VALUE` · `--iso-colors A,B` · `--iso-opacity` · `--mo N` | a cube's isosurface: which file, what level, lobe colours, and which orbital |
| `--frame N` | which geometry of a trajectory |
| `--opaque` · `--background COLOR` | fill the background instead of leaving it transparent |
| `--bond-color COLOR` · `--no-hbonds` | one colour for every bond · hide the dashed hydrogen bonds |
| `--fwhm` · `--freq-scale` · `--unit` · `--xrange` · `--nucleus` · `--reference` | spectrum shape and axes |
| `--no-sticks` · `--csv FILE` · `--dpi` · `--title` | spectrum extras |
| `--show-hydrogens` · `--carbon-labels` · `--lone-pairs` · `--color-atoms` | what the Lewis drawing shows |
| `--spin DEG` | turn the picture in the plane of the page |

`scope --help` is the full list. Atom numbers are 1-based, exactly as the
viewer labels them. `--style`, `--labels`, `--axes` and `--lewis` also work
without `-s`, to open the viewer already set up that way; the render-only
flags are refused in GUI mode rather than silently ignored.

## In the viewer

**Mixed representations.** Select a region and give it its own level of
detail — the active site in full ball-and-stick, spectator ligands as sticks,
everything else as thin lines.

**Move and rotate by hand.** A manipulator sits on whatever is selected: drag
the red/green/blue arrows to slide it along x/y/z, the matching rings to turn
it about that axis, the centre dot to move it in the screen plane. Hold
`Shift` to snap to 0.1 Å / 15°. `F` grows the selection to the whole connected
fragment first, so a substituent is picked up in one keystroke.

![the move/rotate manipulator and the XYZ axis indicator](docs/screenshot_gizmo.png)

**Measurements** are drawn the way a paper draws them: the distance written
along the bond, angles marked with an arc at the vertex, dihedrals with a
rotation arrow around the central bond. `M` pins one into the scene, and
pinned measurements are included in exported images.

![labels, pinned measurements and plane alignment](docs/screenshot_features.png)

**Lewis structure mode** (`Shift+L`) redraws the molecule flat and skeletal
the way ChemDraw does — bare carbon vertices, `OH`/`NH₂` labels with the
hydrogens folded in, double bonds with the second line inside the ring, formal
charges, lone-pair dots on request. It shares the camera with the 3-D view, so
turning it chooses the angle the drawing is made from. The mode is
deliberately read-only, because picking the angle is all it is for.

![the Lewis structure mode](docs/screenshot_lewis.png)

Bond orders and charges are worked out from the geometry, and that cannot
always be right — a carbocation and a carbanion have the same connectivity,
and a metal complex has no Lewis structure at all. A charge the file declares
but the drawing cannot put on any atom goes on square brackets round the whole
structure, the way a delocalized radical cation is drawn. Where it still
cannot be sure, the drawing says so: the perceived charge disagreeing with the
file, a metal whose surrounding charges are only bookkeeping, atoms stacked
behind each other in this view (which just means: turn it).

**Isosurfaces.** Open a `.cube` file (from `cubegen`, `orca_plot`, …) and the
orbital or density surface appears; `I` opens live controls.

![orbital isosurface from a cube file](docs/screenshot_isosurface.png)

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
| MDL molfile (`.mol`/`.sdf`) | — | ✅ 3-D, or the flat drawing with bond orders |
| ChemDraw (`.cdxml`) | — | ✅ the flat drawing at the angle on screen |

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

Done: the viewer, native Gaussian/ORCA/Q-Chem/Molden parsers, trajectory
playback, IR with mode animation, UV-Vis and NMR, geometry editing with undo,
input writers, hydrogen bonds, cube isosurfaces, the Houk style, mixed
representations, the manipulator, Lewis structure mode, and the `scope -s`
command line.

Next: a style editor, and a PyPI release once the distribution name is
settled (`microscope` and `scope` are both taken).

Computing MO/density grids from basis sets was considered and **rejected** —
this draws what quantum-chemistry programs produced, and never computes
chemistry itself.

## Python API

Reading files and everything derived from them is scriptable, and importing
`microscope` starts no graphics — a script that only wants numbers never pays
for Qt:

```python
import microscope

result = microscope.load("mycalc.log")      # any supported format
mol = result.molecule                        # final geometry
print(mol.formula(), result.scf_energies[-1])

for v in result.vibrations[:3]:
    print(v.frequency, v.ir_intensity)

lewis = microscope.perceive(mol)             # bond orders, charges, lone pairs
print(lewis.total_charge, lewis.matches_file)

microscope.save_molecule("next_job.gjf", mol)
```

`microscope.__all__` is the supported surface; anything else is internal.
`examples/api_demo.py` is a fuller tour (trajectories, IR, TD-DFT, NMR,
format conversion).

## Development

```bash
pip install -e ".[dev]"
pytest                            # 247 tests, headless, no display needed
ruff check microscope tests scripts
mypy
python scripts/preview.py         # offscreen render smoke test
python scripts/record_demo.py lewis   # regenerate the README animation
```

Tests run headless on all three platforms via GitHub Actions and cover the
parsers, the geometry engine, editing, Lewis perception and the spectra
models. `scripts/record_demo.py orbit vibrate` records the 3-D clips too, but
needs a working OpenGL 3.3 driver.

## Acknowledgments

- The rendering style is inspired by **CYLview** (C. Y. Legault, Université
  de Sherbrooke) — this project is an independent implementation.
- Parser test fixtures in `tests/data/` come from the **cclib** project's
  test suite (BSD-3-Clause); see `tests/data/README.md`.

## License

MIT
