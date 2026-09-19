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

wavefunction = microscope.load("mycalc.fchk").orbitals
homo = wavefunction.volume("homo")          # the HOMO on a grid, as a VolumeData
microscope.save_cube("homo.cube", homo, mol)
```

The names documented here are the supported surface, listed in
`microscope.__all__`. Anything else is internal and may move between releases.

## Reading and writing

::: microscope.load

::: microscope.save_molecule

::: microscope.save_lewis

::: microscope.save_cube

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

## Orbitals

::: microscope.Orbitals
    options:
      members:
        - find
        - volume
        - name
        - describe
        - restricted
        - spin
        - orthonormality_error

::: microscope.OrbitalSet
    options:
      members:
        - nmo
        - homo

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
