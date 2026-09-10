# Mixed representations

Select a region and give it its own level of detail. The part that matters can
stay in full ball-and-stick while its surroundings drop to sticks, and
everything else to thin lines.

![giving a region its own level of detail](../assets/regions.gif)

| Keys | Representation |
|---|---|
| ++1++ | ball and stick |
| ++2++ | stick |
| ++3++ | line |

With nothing selected, the key applies to the whole molecule. A bond between
two regions takes the thinner of the two, so the boundary does not show a step.

This is the usual way to show an active site: the catalyst and substrate in
ball-and-stick, the ligand framework as sticks, the solvent shell as lines.

## From the command line

```bash
scope -s complex.log --rep '1-12:ball;13-40:stick;41-200:line' -o fig.png
```

`--rep` takes semicolon-separated `atoms:representation` pairs, with 1-based
atom ranges, or a single word for the whole molecule.
