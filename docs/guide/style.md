# Your own style

++ctrl+t++ (View → Style…) opens the style editor. Every change applies to the
view as you make it.

![editing atom size, bond width and element colours](../assets/style.gif)

What you can change:

- **Start from** — the CYLview or Houk preset, or a style you saved earlier
- **Atom size** and **bond width**
- **Background**
- **Bonds** — one colour for every bond, or each bond split by its two atoms
- **Hydrogen bonds** — on or off
- **Elements** — a colour button for each element that is actually in the open
  structure, so tryptophan shows H, C, N and O rather than a periodic table to
  hunt through

## A group standard

The reason to save a style is so that every figure matches. Save it once and
the same file drives both the viewer and the command line:

```bash
scope -s reactant.log --style ourgroup.json -o fig1.png
scope -s ts.log       --style ourgroup.json -o fig2.png
scope -s product.log  --style ourgroup.json -o fig3.png
```

Keep the file with the manuscript and the look travels with it — a co-author,
or you in six months, gets the same figures from the same command.

## The file

A style file is JSON, written to be read and edited by hand: elements are named
by symbol and colours are `#rrggbb`. It only needs the keys you want to change;
everything else falls back to the CYLview preset.

```json
{
  "name": "our group",
  "atom_scale": 0.42,
  "bond_color": "#1a1a1a",
  "palette": {"C": "#4d4d4d", "N": "#2f5fd0", "O": "#d13b2e"}
}
```

A colour may also be written as three numbers from 0 to 1, as in
`[0.3, 0.3, 0.3]`, if that is easier to produce from a script.

If a file is wrong, it says what is wrong rather than raising a traceback:

```text
scope: ours.json: not part of a style: colour
scope: ours.json: palette: 'Xx' is not an element
scope: ours.json: background: 'red' is not #rrggbb
```

Every key, with its default, is listed under
[style files](../reference/style-files.md).

## From Python

```python
import microscope

style = microscope.load_style("ourgroup.json")
style.atom_scale = 0.5
microscope.save_style("bigger.json", style)
```
