# Getting started

## Install

microscope runs anywhere Python 3.10 or later runs — Windows, macOS, Linux —
and needs only six packages, all of which install with `pip`.

```bash
git clone https://github.com/lucienwb/microscope
cd microscope
conda create -n microscope python=3.12      # or any virtual environment
conda activate microscope
pip install -e .
```

That installs the `scope` command. The dependencies are deliberately few:
`numpy`, `scipy`, `matplotlib`, `pandas`, `PySide6` and `PyOpenGL`. Every file
parser is written in-house. If [cclib](https://cclib.github.io) happens to be
installed it is used as a fallback for formats the native parsers do not cover,
but nothing requires it.

!!! note "Not on PyPI yet"
    Both `microscope` and `scope` are taken on PyPI, so for now it installs from
    the repository. `pip install -e .` keeps it editable, so a `git pull` is all
    an update takes.

## Try it now

No calculation of your own to hand? The repository comes with real program
outputs in `tests/data/`, and every one of these opens from the `microscope`
folder you just cloned:

| Command | What to try |
|---|---|
| `scope tests/data/trp.log` | tryptophan: drag to turn it, ++v++ for the other style, ++shift+l++ for the Lewis structure |
| `scope tests/data/dvb_ir.out` | a frequency job: click a band in the IR spectrum to watch that vibration |
| `scope tests/data/dvb_ir.fchk` | ++i++ for the orbitals; the arrow keys walk up and down the list |
| `scope tests/data/g09_dvb_scan.log` | a coordinate scan: the frame bar steps through its geometries |
| `scope tests/data/water_mo.cube` | an orbital from a cube file; ++i++ to change the isovalue and colours |
| `scope tests/data/1crn.pdb` | a small protein: ++2++ for sticks, then ++d++ for depth cueing |

The same files make figures without opening a window:

```bash
scope -s tests/data/trp.log -o trp.png                     # the molecule
scope -s tests/data/dvb_ir.fchk --mo homo -o homo.png      # its HOMO
scope -s tests/data/dvb_ir.out --ir -o ir.png              # its IR spectrum
scope -s tests/data/trp.log --lewis -o trp.svg             # a Lewis structure
```

## Open something

```bash
scope mycalc.log
```

Any Gaussian, ORCA or Q-Chem output works — the program is recognized from the
file's contents, not its extension — as do `xyz`, `pdb`, `fchk`, Gaussian
inputs, Molden and cube files. The [file formats](reference/formats.md) page
lists what is read from each.

Drag to turn it, scroll to zoom, and press ++v++ to see the other style. The
[keyboard reference](reference/keyboard.md) has everything else on one page.

## Make a figure without opening anything

```bash
scope -s mycalc.log -o figure.png
```

`-s` renders and exits. Every choice the viewer offers has a flag, so a figure
can be regenerated from a shell script whenever the calculation changes:

```bash
scope -s opt.log -o fig.tif --style houk --size 2000x1500 --view 30,-15
```

The output is a transparent PNG, or an uncompressed TIFF as a lossless master.
The [command line](reference/command-line.md) page lists every flag.

## Where to go next

- The **[Guide](guide/viewing.md)** walks through each thing the viewer does,
  with an animation of it doing it.
- The **[Reference](reference/keyboard.md)** is for looking things up: keys,
  flags, style-file keys, formats, the Python API.
- **[How it works](about/how-it-works.md)** explains the decisions behind the
  rendering and the Lewis perception, for the curious.
