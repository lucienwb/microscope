---
hide:
  - navigation
---

# microscope

**A molecular structure and spectroscopy viewer for quantum chemistry**, with
CYLview-style rendering. It reads the files your calculations write — Gaussian,
ORCA and Q-Chem outputs, `xyz`, `pdb`, `fchk`, Molden, cube — and turns them
into clean 3-D structures, interactive spectra, editable geometries, and
figures ready for a paper. It installs the short command `scope`.

![tryptophan rendered by microscope, with its intramolecular hydrogen bond](assets/hero.png)

## What it does

<div class="grid cards" markdown>

-   :material-rotate-3d-variant: **Look at a structure**

    ---

    Turn it, frame it, label it — in the CYLview render or the Houk style.

    [:octicons-arrow-right-24: Looking at a structure](guide/viewing.md)

-   :material-palette: **Your own style**

    ---

    Edit colours and sizes live, save them, and share one look across a group.

    [:octicons-arrow-right-24: Your own style](guide/style.md)

-   :material-ruler: **Measure it**

    ---

    Distances, angles and dihedrals, drawn the way a paper draws them.

    [:octicons-arrow-right-24: Selecting and measuring](guide/measuring.md)

-   :material-cursor-move: **Move it**

    ---

    Drag a fragment along an axis or turn it about one, with undo.

    [:octicons-arrow-right-24: Moving and editing](guide/editing.md)

-   :material-molecule: **Draw it flat**

    ---

    ChemDraw-style Lewis structures, saved as `.cdxml`, `.mol`, SVG or PDF.

    [:octicons-arrow-right-24: Lewis structures](guide/lewis.md)

-   :material-blur: **Orbitals**

    ---

    Cube files as isosurfaces, with live isovalue and colours.

    [:octicons-arrow-right-24: Orbitals and densities](guide/orbitals.md)

-   :material-play-circle: **Watch it move**

    ---

    Optimizations, IRCs and scans played back; normal modes animated.

    [:octicons-arrow-right-24: Trajectories and vibrations](guide/trajectories.md)

-   :material-chart-bell-curve: **Spectra**

    ---

    IR, UV-Vis and NMR — click a band and watch that mode.

    [:octicons-arrow-right-24: Spectra](guide/spectra.md)

</div>

## Three ways to use it

=== "The viewer"

    ```bash
    scope mycalc.log
    ```

    Opens the window. Everything is on the keyboard — see the
    [keyboard reference](reference/keyboard.md).

=== "From a script"

    ```bash
    scope -s mycalc.log --style houk --view 30,-15 -o fig.png
    scope -s freq.log --ir -o ir.pdf
    scope -s mol.log --lewis -o scheme.cdxml
    ```

    No window opens: the figure is written and the command returns, the way
    `gnuplot` makes a plot. Put it in a Makefile and every figure regenerates
    itself after a re-optimization. See the [command line](reference/command-line.md).

=== "From Python"

    ```python
    import microscope

    result = microscope.load("mycalc.log")
    print(result.molecule.formula(), result.scf_energies[-1])
    ```

    The parsers and writers with no graphics attached, so a script that only
    wants the numbers never starts Qt. See the [Python API](reference/python-api.md).

## Install

```bash
git clone https://github.com/lucienwb/microscope
cd microscope
pip install -e .
scope mycalc.log
```

It runs anywhere Python 3.10 or later does, and needs only `numpy`, `scipy`,
`matplotlib`, `pandas`, `PySide6` and `PyOpenGL`. More in
[Getting started](getting-started.md).
