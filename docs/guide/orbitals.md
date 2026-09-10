# Orbitals and densities

Open a `.cube` file — from Gaussian's `cubegen`, ORCA's `orca_plot`, or anything
else that writes the format — and the orbital or density appears as an
isosurface.

![sweeping the isovalue of an orbital](../assets/isosurface.gif)

++i++ opens live controls:

- **Isovalue** — where the surface is drawn. Orbitals start at the conventional
  |ψ| = 0.02 au
- **Opacity**
- **Colours** — curated pairs for the + and − lobes, or any RGB or hex colour
  for each
- **Which orbital**, for a cube holding several

Surfaces are translucent with a depth pre-pass, so only the nearest layer is
shaded and back faces do not show through as noise. They are included when the
view is exported.

microscope draws what a cube file contains; it does not compute orbitals or
densities from a basis set. Generate the cube with the program that ran the
calculation.

## From the command line

```bash
scope -s homo.cube --iso 0.03 --iso-colors purple,gold --zoom 1.3 -o homo.png
scope -s opt.log --cube orbitals.cube --mo 3 -o mo3.png
```
