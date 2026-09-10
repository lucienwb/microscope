# Selecting and measuring

## Selecting

Click an atom to add it to the selection, click it again to drop it, and click
empty space or press ++esc++ to clear. There is no limit on how many are
selected — larger selections are what the [representation](representations.md)
and [editing](editing.md) commands act on.

## Measuring

With two, three or four atoms selected, the distance, angle or dihedral is shown
live, in the status bar and in the scene. It is drawn the way a paper draws it:

- a **distance** written along the bond it measures
- an **angle** marked with an arc at the vertex, the value beside it
- a **dihedral** with a rotation arrow around the central bond — seen straight
  down that bond, it reads as a Newman projection

![measuring a distance, an angle and a dihedral](../assets/measure.gif)

The status bar names the atoms the way the viewer labels them — `d(C1–C2) =
1.421 Å`, `∠(C1–C2–C3) = 120.9°` — so a value can be copied straight into a
table.

## Pinning

++m++ pins the current measurement into the scene, so it stays while you select
the next one. Pin as many as you like; ++shift+m++ clears them. Pinned
measurements are included when the view is exported.

![labels, pinned measurements and plane alignment](../assets/screenshot_features.png)

## From the command line

```bash
scope -s mol.xyz --measure 3,4 --measure 3,4,5 -o fig.png
```

Each `--measure` takes two, three or four 1-based atom numbers.
