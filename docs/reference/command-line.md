# Command line

```bash
scope                  # open the viewer
scope mycalc.log       # open a file in it
scope -s mycalc.log    # no window: render mycalc.png and exit
```

`-s` / `--silent` turns microscope into a batch plotter. No window opens; the
figure is written and the command returns, so a figure can live in a shell
script or a Makefile and regenerate itself after every re-optimization.

Without `-o`, the output is named after the input. Atom numbers are 1-based,
exactly as the viewer labels them.

!!! tip "A value that starts with a minus sign"
    Write `--rotate=-40,-70`, with an `=`. Without it, `-40,-70` looks like a
    flag of its own and the command stops with *expected one argument*.

## What `-s` can make

=== "A figure of the molecule"

    ```bash
    scope -s opt.log -o fig.tif --style houk --size 2000x1500
    scope -s mol.xyz --view top --labels number --measure 3,4
    scope -s mol.log --rep '1-12:ball;13-40:line' --view 30,-15 --axes
    scope -s homo.cube --iso 0.03 --iso-colors purple,gold --zoom 1.3
    scope -s scan.log --frame 12 --align 3,4 -o frame12.png
    ```

    A transparent PNG, or an uncompressed TIFF. `--opaque` fills the background.

=== "A spectrum"

    ```bash
    scope -s freq.log --ir -o ir.pdf --fwhm 12 --freq-scale 0.965
    scope -s td.log   --uv --unit eV --xrange 2:7 -o uv.svg
    scope -s nmr.log  --nmr --nucleus C --reference 186.4 --csv shifts.csv
    ```

    PNG, PDF, SVG, EPS or TIFF. Only matplotlib is loaded.

=== "A Lewis structure"

    ```bash
    scope -s mol.log --lewis -o scheme.svg
    scope -s mol.log --lewis -o mol.cdxml
    scope -s mol.log --lewis --lone-pairs --color-atoms -o lewis.png
    ```

    A picture, vector art, or a `.cdxml` / `.mol` for ChemDraw. The structure
    files need no graphics library at all.

`--style`, `--labels`, `--axes` and `--lewis` also work without `-s`, to open the
viewer already set up that way. Flags that only make sense when rendering are
refused in the viewer rather than silently ignored, and flags for a structure
and for a spectrum cannot be mixed.

## Every flag

This table is generated from the program's own argument parser, so it cannot
fall out of date. `scope --help` prints the same.

<!-- command-line-reference -->
