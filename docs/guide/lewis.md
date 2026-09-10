# Lewis structures

++shift+l++ (View → Lewis Structure) swaps the 3-D model for a flat, skeletal
drawing, the way ChemDraw draws it: bare carbon vertices, `OH` and `NH₂` labels
with their hydrogens folded in, the second line of a ring double bond drawn
inside the ring, and formal charges.

## Picking the angle

The drawing shares its camera with the 3-D view, so **turning it is how you
choose the angle it is drawn from**. Drag to turn, ++shift++-drag to spin it in
the plane of the page, scroll to zoom.

![turning a Lewis structure to pick the angle it is drawn from](../assets/lewis_spin.gif)

The mode is deliberately read-only — the Edit menu greys out — because picking
the angle is all it is for. When two atoms land on top of each other, a note in
the corner says so, which only ever means: turn it.

![the Lewis structure mode beside the 3-D view](../assets/screenshot_lewis.png)

## What it shows

Four switches in View → Lewis Structure decide how much is drawn:

![the Lewis drawing options](../assets/lewis_options.gif)

- **Colour atoms** — heteroatoms in their element colour
- **Show lone pairs** — as dots around each atom
- **Show all hydrogens** — rather than folding them into the labels
- **Label carbons** — rather than leaving them as bare vertices

## Saving it

- ++ctrl+e++ saves a **picture**: PNG or TIFF, or — because this is line art —
  **SVG or PDF**, with the lettering kept as real text
- ++ctrl+shift+s++ saves for **ChemDraw**: a `.cdxml`, or an MDL molfile, at the
  angle on screen, carrying the bond orders, charges and implicit hydrogens, so
  it opens ready to rearrange rather than as a pile of atoms

From the command line, `-o` decides which:

```bash
scope -s mol.log --lewis -o scheme.svg        # vector art
scope -s mol.log --lewis -o mol.cdxml         # opens in ChemDraw at this angle
scope -s mol.log --lewis --lone-pairs --color-atoms -o lewis.png
```

Writing a `.cdxml` or `.mol` needs no graphics at all, so it works on a cluster
node with no display.

## How the charges are worked out

A Lewis structure is perceived from the geometry — nothing is recomputed, the
file is read exactly as the program wrote it. That cannot always be right, so
the drawing is careful to say what it is claiming and where it gave up.

**Bond orders** come from bond lengths. Every pair of elements has known single,
double and triple lengths, so the nearest reference is what a bond wants to be;
those wishes are then reconciled with the valence each atom has. Aromatic rings
come out as Kekulé structures, including fused ones — pentacene, perylene,
azulene, corannulene.

**An ion** has the same connectivity as the neutral molecule, so its charge
cannot be seen in the geometry. But the file says what the charge is, so it is
placed: on the atom if there is only one it can belong to — a bare chloride —
and otherwise on **square brackets round the whole structure**, the way a
delocalized radical cation is drawn.

**A metal complex** is drawn with ionic counting. A dative bond drawn as a plain
line would make a phosphine into a phosphonium, which is not how anyone draws a
complex, so each ligand is counted without its bond to the metal: the phosphine
keeps its lone pair and stays neutral, a chloride bonded only to the metal is a
chloride ion, and whatever is left over goes on the metal. Counted that way, the
number on the metal **is its oxidation state** — Rh(+1) on a hydroformylation
catalyst, Mn(+1) on a pincer carbonyl, Mo(+4) in MoOCl₄²⁻. A note in the corner
says which convention the drawing is using.

**A crystal structure** usually has no hydrogens. Counting the bonds they would
have made as formal charges buries a protein in carbanions, so when most
carbons are short of four bonds the drawing claims no charges at all, and says
why.

If, after all that, the drawing's total charge still disagrees with the file,
a note says so rather than quietly showing something wrong.
