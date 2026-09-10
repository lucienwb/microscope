# Python API

Everything the viewer does with files works from Python, and importing
`microscope` starts no graphics — a script that only wants the numbers never
pays for Qt.

```python
import microscope

result = microscope.load("mycalc.log")      # any supported format
mol = result.molecule                        # the final geometry
print(mol.formula(), result.scf_energies[-1])

for v in result.vibrations[:3]:
    print(v.frequency, v.ir_intensity)

lewis = microscope.perceive(mol)             # bond orders, charges, lone pairs
print(lewis.total_charge, lewis.matches_file)

microscope.save_molecule("next_job.gjf", mol)
```

The names documented here are the supported surface, listed in
`microscope.__all__`. Anything else is internal and may move between releases.

## Reading and writing

::: microscope.load

::: microscope.save_molecule

::: microscope.save_lewis

## What a file gives you

::: microscope.ParseResult

::: microscope.Molecule
    options:
      members:
        - natoms
        - atomic_numbers
        - perceive_bonds
        - formula
        - bounding_sphere

::: microscope.Vibration

::: microscope.ExcitedState

::: microscope.NMRShielding

::: microscope.VolumeData

## Lewis structures

::: microscope.perceive

::: microscope.LewisStructure
    options:
      members:
        - total_charge
        - bracketed
        - matches_file
        - heavy_neighbors
        - order_of

## Styles

::: microscope.Style

::: microscope.load_style

::: microscope.save_style

## Errors

::: microscope.FileFormatError

::: microscope.UnsupportedFormatError

::: microscope.StyleError
