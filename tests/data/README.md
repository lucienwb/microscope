# Test data

Parser test fixtures.

- `dvb_*.out/.log/.gjf/.fchk`, `orca_dvb_*.out`, `qchem_dvb_*.out`,
  `water_mp2.log`, `trp.log`, `g09_*.log` (Gaussian 09: unrestricted SP,
  relaxed scan, Mo4OCl4²⁻ single point), `orca5_dvb_gopt.out` (ORCA 5.0),
  `qchem6_water_smd.out` (Q-Chem 6.0) and `dvb_un_sp.fchk` are taken from
  the test suite of the [cclib](https://github.com/cclib/cclib) project
  (BSD-3-Clause license) and are used here solely as parser test inputs.
  We thank the cclib developers for maintaining this collection of real
  program outputs.
- `water.molden` is a hand-written minimal Molden file created for this
  project.
