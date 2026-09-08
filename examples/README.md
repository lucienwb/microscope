# examples

`api_demo.py` is the tracked one: a short tour of the Python API that parses a
Gaussian output, reads its spectra and writes the geometry back out in three
formats. Run it from the repo root — it drops its files in `output/`.

```bash
python examples/api_demo.py
```

## examples/data — the validation corpus (not in the repo)

`data/` holds ~600 real structures used to check the parsers, the bond
perception and the Lewis drawing against things people actually calculate. It
is **gitignored**: most of it is downloaded, some of it is unpublished work,
and none of it belongs in a source tree. Fetch or drop in your own; nothing in
the test suite depends on it.

| folder | what is in it | where it came from |
|---|---|---|
| `catalysis/` | Mn PNP pincer and Rh hydroformylation intermediates | the maintainer's own calculations |
| `programs/<Program>/` | outputs from 19 quantum-chemistry codes | [cclib-data](https://github.com/cclib/cclib-data) and [cclib](https://github.com/cclib/cclib) test suites |
| `formats/<ext>/` | one folder per format microscope reads natively | [iodata](https://github.com/theochem/iodata) test data, cclib |
| `biomolecules/` | 72 PDB entries: enzymes, metalloproteins, photosystems, nucleic acids, whole assemblies | [RCSB](https://www.rcsb.org/) |
| `metal-complexes/` | 56 transition-metal complexes, two per metal across 28 metals | [tmQM](https://github.com/uiocompcat/tmQM) (from the CSD) |
| `photocatalysts/` | organic photocatalysts, dyes and chromophores | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) 3D conformers |
| `unsupported/` | formats microscope deliberately refuses (`.in`, `.inp`, `.sdf`, `.mol2`, `.gro`, `.wfn`) | the same sources |

Each source has its own licence and citation requirements — check them before
using any of this for anything but testing.
