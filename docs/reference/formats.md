# File formats

Every reader is written in-house. Program outputs are recognized by their
contents, not their extension — any `.log` or `.out` is sniffed and handed to
the right parser.

| Format | Read | Write |
|---|---|---|
| Gaussian output (`.log`, `.out`) | geometries, SCF energies, frequencies with IR and normal modes, TD-DFT, NMR shieldings | — |
| Gaussian input (`.gjf`, `.com`) | Cartesian or Z-matrix, with named variables | ✓ |
| Gaussian formatted checkpoint (`.fchk`) | geometry, and the wavefunction: basis set and orbitals | — |
| ORCA output (`.out`) | geometries, energies, frequencies with IR and normal modes, TD-DFT, NMR shieldings | — |
| ORCA input (`.inp`) | — | ✓ |
| Q-Chem output (`.out`) | geometries, energies, frequencies with IR and normal modes, TD-DFT | — |
| Q-Chem input (`.in`) | — | ✓ |
| Molden (`.molden`, ORCA's `.molden.input`) | geometries, frequencies and normal modes, and the wavefunction | — |
| Cube (`.cube`, `.cub`) | geometry and volumetric grids, including several orbitals in one file | ✓ an orbital, or any grid on show |
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

**Wavefunctions.** From an fchk or Molden file microscope reads the basis set
(s to g Cartesian, pure functions to any angular momentum, SP shells) and the
orbitals with their energies and occupations, restricted or unrestricted, so it
can [draw any orbital](../guide/orbitals.md) without a cube file. Molden files
are written differently by different programs; each reading is checked against
the file's own orbitals, which must come out orthonormal, and the one that does
is kept. Ghost atoms keep their basis functions but are not drawn. A Molden
atom is named by its symbol, not its nuclear charge: with a pseudopotential,
ORCA writes the core charge there, 17 for rhodium. A file that writes pure
functions without saying so (no `[5D]`) is read by counting its coefficients,
and one whose orbitals cannot be read at all still opens as a structure.

**Cube files** are read in single precision. The format carries five
significant digits, and a large grid takes half the memory.

**Big structures.** A PDB entry of tens of thousands of atoms opens in a second
or two: bonds are found through a grid of cells rather than by comparing every
pair of atoms.

## Other programs

ADF, GAMESS, NWChem, Psi4, Molpro, Turbomole, Jaguar, DALTON and the rest are not
read natively. If [cclib](https://cclib.github.io) is installed, microscope uses
it as a fallback for them; otherwise it says clearly that it cannot read the file.
