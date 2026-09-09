"""Parser tests against real sample outputs (from the cclib test suite)."""

from pathlib import Path

import numpy as np
import pytest

import microscope.io as mio
from microscope.core import geometry
from microscope.io import fchk, gaussian, orca, qchem

DATA = Path(__file__).parent / "data"

needs = pytest.mark.skipif(not DATA.is_dir(), reason="tests/data missing")


@needs
def test_gaussian_opt_trajectory():
    result = gaussian.read_log(DATA / "dvb_gopt.out")
    assert result.program == "Gaussian"
    assert result.normal_termination
    mol = result.molecule
    assert mol.formula() == "C10H10"          # divinylbenzene
    assert result.nframes > 1                 # optimization steps
    assert len(result.scf_energies) >= result.nframes - 1
    mol.perceive_bonds()
    assert 20 <= len(mol.bonds) <= 22         # 21 bonds in DVB


@needs
def test_gaussian_frequencies():
    result = gaussian.read_log(DATA / "dvb_ir.out")
    mol = result.molecule
    nmodes_expected = 3 * mol.natoms - 6
    assert len(result.vibrations) == nmodes_expected
    v = result.vibrations[0]
    assert v.ir_intensity is not None
    assert v.displacements is not None
    assert v.displacements.shape == (mol.natoms, 3)
    freqs = [v.frequency for v in result.vibrations]
    assert freqs == sorted(freqs)


@needs
def test_gaussian_tddft():
    result = gaussian.read_log(DATA / "dvb_td.out")
    assert len(result.excited_states) >= 5
    s1 = result.excited_states[0]
    assert 0 < s1.energy_ev < 10
    assert s1.wavelength_nm > 100
    assert s1.osc_strength >= 0


@needs
def test_gaussian_nmr():
    result = gaussian.read_log(DATA / "dvb_nmr.log")
    assert len(result.nmr_shieldings) == result.molecule.natoms
    carbons = [s for s in result.nmr_shieldings if s.symbol == "C"]
    assert len(carbons) == 10


@needs
def test_gaussian_input_roundtrip(tmp_path):
    mol = gaussian.read_log(DATA / "dvb_ir.out").molecule
    out = tmp_path / "roundtrip.gjf"
    gaussian.write_gjf(out, mol)
    again = gaussian.read_gjf(out).molecule
    assert again.symbols == mol.symbols
    np.testing.assert_allclose(again.coords, mol.coords, atol=1e-6)


@needs
def test_gaussian_input_allcheck_rejected():
    from microscope.io.errors import FileFormatError
    with pytest.raises(FileFormatError, match="checkpoint"):
        gaussian.read_gjf(DATA / "dvb_ir.gjf")   # geom=allcheck: no coordinates


@needs
def test_fchk():
    result = fchk.read(DATA / "dvb_ir.fchk")
    assert result.molecule.formula() == "C10H10"


@needs
def test_xyz_roundtrip(tmp_path):
    mol = gaussian.read_log(DATA / "dvb_ir.out").molecule
    out = tmp_path / "mol.xyz"
    mio.save_molecule(out, mol)
    again = mio.load(out).molecule
    assert again.symbols == mol.symbols
    np.testing.assert_allclose(again.coords, mol.coords, atol=1e-6)


@needs
def test_pdb_roundtrip(tmp_path):
    mol = gaussian.read_log(DATA / "dvb_ir.out").molecule
    out = tmp_path / "mol.pdb"
    mio.save_molecule(out, mol)
    again = mio.load(out).molecule
    assert again.symbols == mol.symbols
    np.testing.assert_allclose(again.coords, mol.coords, atol=1e-3)  # pdb: 3 decimals


@needs
def test_water():
    result = gaussian.read_log(DATA / "water_mp2.log")
    assert result.molecule.formula() == "H2O"


# ---------------------------------------------------------------- ORCA

@needs
def test_orca_opt():
    result = mio.load(DATA / "orca_dvb_gopt.out")   # sniffed by content
    assert result.program == "ORCA"
    assert result.normal_termination
    assert result.molecule.formula() == "C10H10"
    assert result.nframes > 1
    assert result.molecule.charge == 0 and result.molecule.multiplicity == 1
    assert len(result.scf_energies) >= 3


@needs
def test_orca_frequencies():
    result = orca.read_log(DATA / "orca_dvb_ir.out")
    mol = result.molecule
    assert len(result.vibrations) == 3 * mol.natoms - 6
    strongest = max(result.vibrations, key=lambda v: v.ir_intensity or 0)
    assert strongest.ir_intensity > 10
    v = result.vibrations[0]
    assert v.displacements is not None
    assert v.displacements.shape == (mol.natoms, 3)


@needs
def test_orca_tddft():
    result = orca.read_log(DATA / "orca_dvb_td.out")
    assert len(result.excited_states) >= 5
    assert all(s.energy_ev > 0 for s in result.excited_states)
    assert any(s.osc_strength > 0.01 for s in result.excited_states)


@needs
def test_orca_nmr():
    result = orca.read_log(DATA / "orca_dvb_nmr.out")
    assert len(result.nmr_shieldings) == result.molecule.natoms
    carbons = [s for s in result.nmr_shieldings if s.symbol == "C"]
    assert len(carbons) == 10
    assert 100 < carbons[0].isotropic < 130


# ---------------------------------------------------------------- Q-Chem

@needs
def test_qchem_opt():
    result = mio.load(DATA / "qchem_dvb_gopt.out")
    assert result.program == "Q-Chem"
    assert result.molecule.formula() == "C10H10"
    assert result.nframes > 1
    assert len(result.scf_energies) >= 3


@needs
def test_qchem_frequencies():
    result = qchem.read_log(DATA / "qchem_dvb_ir.out")
    mol = result.molecule
    assert len(result.vibrations) == 3 * mol.natoms - 6
    v = result.vibrations[0]
    assert v.ir_intensity is not None
    assert v.displacements is not None
    assert v.displacements.shape == (mol.natoms, 3)


@needs
def test_qchem_tddft():
    result = qchem.read_log(DATA / "qchem_dvb_td.out")
    assert len(result.excited_states) >= 5
    assert any(s.osc_strength > 0.01 for s in result.excited_states)
    assert any("Singlet" in s.label or "Triplet" in s.label for s in result.excited_states)


# ---------------------------------------------------------------- Molden

@needs
def test_molden():
    result = mio.load(DATA / "water.molden")
    mol = result.molecule
    assert mol.formula() == "H2O"
    assert len(result.vibrations) == 3
    assert abs(result.vibrations[0].frequency - 1595.0) < 1e-6
    assert result.vibrations[0].ir_intensity == 65.0
    assert result.vibrations[0].displacements.shape == (3, 3)


# ---------------------------------------------------------------------------
# additional program versions / job types (fixtures also from cclib)


@needs
def test_gaussian09_unrestricted_cation():
    result = gaussian.read_log(DATA / "g09_dvb_un_sp.log")
    assert result.normal_termination
    mol = result.molecule
    assert mol.formula() == "C10H10"
    assert (mol.charge, mol.multiplicity) == (1, 2)   # DVB radical cation
    assert abs(result.scf_energies[-1] - (-382.081395204)) < 1e-8


@needs
def test_gaussian09_relaxed_scan_trajectory():
    result = gaussian.read_log(DATA / "g09_dvb_scan.log")
    assert result.normal_termination
    assert result.nframes == 61                       # every scan opt step
    assert len(result.scf_energies) == 61
    assert abs(result.scf_energies[-1] - (-382.308259560)) < 1e-8
    # geometries actually change along the scan
    d0 = result.frames[0].coords
    dN = result.frames[-1].coords
    assert np.abs(d0 - dN).max() > 0.1


@needs
def test_gaussian09_transition_metal():
    result = gaussian.read_log(DATA / "g09_mo4ocl4_sp.log")
    mol = result.molecule
    assert mol.formula() == "Cl4MoO"
    assert (mol.charge, mol.multiplicity) == (-2, 1)
    assert "Mo" in mol.symbols
    assert abs(result.scf_energies[-1] - (-202.713622575)) < 1e-8


@needs
def test_orca5_opt():
    result = orca.read_log(DATA / "orca5_dvb_gopt.out")
    assert result.program == "ORCA"
    assert result.normal_termination
    assert result.molecule.formula() == "C10H10"
    assert result.nframes == 4
    assert abs(result.scf_energies[-1] - (-382.055133372796)) < 1e-9


@needs
def test_qchem6_solvent_sp():
    result = qchem.read_log(DATA / "qchem6_water_smd.out")
    assert result.program == "Q-Chem"
    assert result.molecule.formula() == "H2O"
    assert abs(result.scf_energies[-1] - (-74.9672988069)) < 1e-8


@needs
def test_fchk_unrestricted():
    result = fchk.read(DATA / "dvb_un_sp.fchk")
    mol = result.molecule
    assert mol.formula() == "C10H10"
    assert (mol.charge, mol.multiplicity) == (1, 2)
    mol.perceive_bonds()
    assert 20 <= len(mol.bonds) <= 22


def test_qchem_skips_ghost_atoms_and_fragment_blocks(tmp_path):
    """A counterpoise job prints the fragments, with ghosts, after the whole.

    Taking the last block left the viewer holding one sodium atom, and the
    basis-only centres came through as element X.
    """
    def orientation(rows):
        return ("             Standard Nuclear Orientation (Angstroms)\n"
                "    I     Atom           X                Y                Z\n"
                " ----------------------------------------------------------\n"
                + rows +
                " ----------------------------------------------------------\n")

    out = tmp_path / "cp.out"
    out.write_text(
        "$molecule\n1 2\n$end\n"
        + orientation("    1      C      -1.4476    0.0000    0.0000\n"
                      "    2      H      -1.5627    0.3304   -1.0238\n"
                      "    3      Na      1.2156    0.0000    0.0000\n")
        + orientation("    1      C      -1.4476    0.0000    0.0000\n"
                      "    2      H      -1.5627    0.3304   -1.0238\n"
                      "    3      GH      1.2156    0.0000    0.0000\n")
        + orientation("    1      Na      1.2156    0.0000    0.0000\n"))

    mol = qchem.read_log(out).molecule
    assert mol.symbols == ["C", "H", "Na"]
    assert "X" not in mol.symbols
    assert mol.charge == 1 and mol.charge_known


def test_formats_without_a_charge_field_say_so(tmp_path):
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\ntitle\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n")
    assert not mio.load(xyz).molecule.charge_known

def test_gaussian_z_matrix_input(tmp_path):
    """Internal coordinates, the way Gaussian test decks are written."""
    job = tmp_path / "water.com"
    job.write_text("#p hf/sto-3g\n\nwater\n\n0 1\nO\nH 1 0.956\n"
                   "H 1 0.956 2 104.5\n\n")
    mol = mio.load(job).molecule
    assert mol.symbols == ["O", "H", "H"]
    assert geometry.distance(mol.coords[0], mol.coords[1]) == pytest.approx(0.956)
    assert geometry.angle(mol.coords[1], mol.coords[0],
                          mol.coords[2]) == pytest.approx(104.5)


def test_gaussian_z_matrix_with_named_variables(tmp_path):
    job = tmp_path / "water.com"
    job.write_text("#n hf/sto-3g\n\nwater\n\n0 1\nO1\nH2  1  r2\n"
                   "H3  1  r3  2  a3\n\nr2=0.9732\nr3=0.9641\na3=105.9\n\n")
    mol = mio.load(job).molecule
    assert geometry.distance(mol.coords[0], mol.coords[1]) == pytest.approx(0.9732)
    assert geometry.angle(mol.coords[1], mol.coords[0],
                          mol.coords[2]) == pytest.approx(105.9)


def test_a_dihedral_in_a_z_matrix_comes_back_out(tmp_path):
    job = tmp_path / "butane.com"
    job.write_text("#n hf\n\nb\n\n0 1\nC\nC 1 1.53\nC 2 1.53 1 111.0\n"
                   "C 3 1.53 2 111.0 1 60.0\n\n")
    mol = mio.load(job).molecule
    assert geometry.dihedral(*mol.coords[::-1]) == pytest.approx(60.0)


def test_pdb_reads_the_formal_charge_columns(tmp_path):
    """Columns 79-80 are usually blank, but they are the only charge PDB has."""
    plain = tmp_path / "plain.pdb"
    plain.write_text("HETATM    1  N   UNL     1       0.000   0.000   0.000"
                     "  1.00  0.00           N\n")
    assert not mio.load(plain).molecule.charge_known

    charged = tmp_path / "charged.pdb"
    charged.write_text("HETATM    1  N   UNL     1       0.000   0.000   0.000"
                       "  1.00  0.00           N1+\n")
    mol = mio.load(charged).molecule
    assert mol.charge == 1 and mol.charge_known

def test_a_value_gaussian_could_not_print_is_missing_not_fatal(tmp_path):
    """Gaussian fills a field with asterisks when the number will not fit.

    A whole calculation used to be refused over a Raman activity nobody had
    asked for.
    """
    job = tmp_path / "freq.log"
    job.write_text(
        " Entering Gaussian System\n #p freq\n\n Charge =  0 Multiplicity = 1\n"
        "                         Standard orientation:\n"
        " ---------------------------------------------------------\n"
        " Center     Atomic     Atomic              Coordinates\n"
        " Number     Number      Type              X        Y        Z\n"
        " ---------------------------------------------------------\n"
        "    1          1             0        0.000000  0.000000  0.000000\n"
        "    2          1             0        0.740000  0.000000  0.000000\n"
        " ---------------------------------------------------------\n"
        " Frequencies --  4400.0000\n"
        " Red. masses --     1.0000\n"
        " IR Inten    --***********\n"
        " Raman Activ --***********\n"
        " Atom AN      X      Y      Z\n"
        "   1   1     0.00   0.00   0.71\n"
        "   2   1     0.00   0.00  -0.71\n"
        " Normal termination of Gaussian\n")
    result = mio.load(job)
    assert len(result.vibrations) == 1
    mode = result.vibrations[0]
    assert mode.frequency == pytest.approx(4400.0)
    assert mode.ir_intensity is None and mode.raman_activity is None
