# File formats

Every reader is written in-house. Program outputs are recognized by their
contents, not their extension — any `.log` or `.out` is sniffed and handed to
the right parser.

| Format | Read | Write |
|---|---|---|
| Gaussian output (`.log`, `.out`) | geometries, SCF energies, frequencies with IR and normal modes, TD-DFT, NMR shieldings | — |
| Gaussian input (`.gjf`, `.com`) | Cartesian or Z-matrix, with named variables | ✓ |
| Gaussian formatted checkpoint (`.fchk`) | geometry | — |
| ORCA output (`.out`) | geometries, energies, frequencies with IR and normal modes, TD-DFT, NMR shieldings | — |
| ORCA input (`.inp`) | — | ✓ |
| Q-Chem output (`.out`) | geometries, energies, frequencies with IR and normal modes, TD-DFT | — |
| Q-Chem input (`.in`) | — | ✓ |
| Molden (`.molden`) | geometries, frequencies and normal modes | — |
| Cube (`.cube`, `.cub`) | geometry and volumetric grids, including several orbitals in one file | — |
| XYZ (`.xyz`) | single or multi-frame | ✓ |
| PDB (`.pdb`) | the first model, and per-atom formal charges | ✓ |
| MDL molfile (`.mol`, `.sdf`) | — | ✓ 3-D, or the Lewis drawing with bond orders |
| ChemDraw (`.cdxml`) | — | ✓ the Lewis drawing, at the angle on screen |

## Details worth knowing

**Charges.** Gaussian, ORCA and Q-Chem outputs and Gaussian inputs state a
molecule's charge and multiplicity. XYZ, cube and Molden files do not, and PDB
only has a per-atom charge that is almost always blank. When a file does not say,
the viewer does not pretend it said zero — which matters for the
[Lewis structure](../guide/lewis.md), where a charge is only checked against the
file when the file stated one.

**Gaussian.** Older versions label a Z-matrix job's geometry "Z-Matrix
orientation"; that is read. Where a number overflows its column and Gaussian
prints asterisks, that one value is recorded as missing rather than the whole
file being refused.

**Q-Chem.** A fragment or counterpoise job prints each piece after the whole
system; only the whole system is kept, and ghost atoms (`GH`) are skipped.

**Cube files** are read in single precision. The format carries five
significant digits, and a large grid takes half the memory.

**Big structures.** A PDB entry of tens of thousands of atoms opens in a second
or two: bonds are found through a grid of cells rather than by comparing every
pair of atoms.

## Other programs

ADF, GAMESS, NWChem, Psi4, Molpro, Turbomole, Jaguar, DALTON and the rest are not
read natively. If [cclib](https://cclib.github.io) is installed, microscope uses
it as a fallback for them; otherwise it says clearly that it cannot read the file.
