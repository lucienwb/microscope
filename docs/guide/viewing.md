# Looking at a structure

## Turning and framing

Left-drag turns the molecule, right-drag pans, and the wheel zooms.

![turning a molecule](../assets/orbit.webp){ .render }

To get a particular view rather than a dragged one, select atoms first:

| Keys | Selection | What happens |
|---|---|---|
| ++a++ | two atoms | look straight down the bond between them |
| ++a++ | three atoms | put their plane in the screen |
| ++c++ | one atom | rotate about that atom |
| ++home++ | — | re-centre on the whole molecule |
| ++ctrl+r++ | — | reset the view entirely |

The camera is orthographic, so there is no perspective to distort lengths:
a bond the same length on screen is the same length in the molecule.

## Two looks

++v++ switches between the two built-in styles.

![switching between the CYLview look and the Houk style](../assets/styles.webp){ .render }

**CYLview** is the default: ray-cast spheres and cylinders, so nothing is ever
faceted however far you zoom, with bonds split half-and-half in the colours of
the atoms they join.

**Houk** ("Houkmol") is the Houk group's ball-and-stick: glossy, with black
bonds, near-white carbons, and two black seam lines on every heavy atom that
turn with the molecule. It stays readable in a black-and-white printout.

Neither has to be the last word — the [style editor](style.md) starts from
either and saves your own.

## Labels and axes

++l++ cycles the atom labels: element, element and number, number, and off.
Atom numbers are 1-based, matching the numbering the command line uses.

![cycling the atom labels](../assets/labels.webp){ .render }

++shift+a++ puts a triad of the world axes in the corner. It turns with the
molecule, and it is included in exported images, so a figure can say which
way is which.

![the XYZ axis triad following the rotation](../assets/axes.webp){ .render }

## Depth cueing

In a big structure, the atoms at the back get in the way of the ones at the
front. ++d++ fades the far side of the molecule toward the background, so
what is near you stands out.

![a small protein in sticks, with depth cueing off and then on](../assets/depth_cue.webp)

It is off until you ask for it, because plenty of people prefer every atom
at full strength, and a figure for print is often clearer that way. How
strongly it fades is the **Depth cueing** slider in the
[style editor](style.md): the front of the molecule always keeps its own
colour, and the slider sets how far toward the background the back goes.
The fade is toward the style's background colour even in a transparent
export, so for a figure going on a dark slide, set the background to match
the slide. The same effect is `--depth-cue 0.5` on the
[command line](../reference/command-line.md), and `"depth_cue"` in a
[style file](../reference/style-files.md), so a group style can switch it on
for everyone.

## Exporting

++ctrl+e++ exports the view as a high-resolution image, with the isosurface,
labels, pinned measurements and axes all included. Two formats, both chosen
because they keep a transparent background:

- **PNG** for everything, transparent by default
- **uncompressed TIFF** as a lossless master for a publisher

For the same thing without the window, see
[the command line](../reference/command-line.md).
