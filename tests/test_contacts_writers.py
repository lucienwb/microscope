"""Hydrogen-bond detection and ORCA/Q-Chem input writer tests."""

from pathlib import Path

import numpy as np

import microscp.io as mio
from microscp.core.contacts import find_hbonds
from microscp.core.molecule import Molecule
from microscp.io import orca, qchem
from microscp.render.scene import build_scene
from microscp.render.styles import Style

DATA = Path(__file__).parent / "data"


def test_water_dimer_hbond():
    mol = mio.load(DATA / "water_dimer.xyz").molecule
    hbonds = find_hbonds(mol)
    assert len(hbonds) == 1
    h, acc, dist = hbonds[0]
    assert (h, acc) == (1, 3)                 # donor H -> acceptor O
    assert abs(dist - 1.89) < 0.01


def test_no_hbond_in_methane():
    mol = Molecule(
        ["C", "H", "H", "H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
                  [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]]))
    assert find_hbonds(mol) == []


def test_hbond_dashes_in_scene():
    mol = mio.load(DATA / "water_dimer.xyz").molecule
    mol.perceive_bonds()
    with_hb = build_scene(mol, Style(show_hbonds=True))
    without = build_scene(mol, Style(show_hbonds=False))
    assert len(with_hb.cylinders) > len(without.cylinders)   # dash segments added
    dashes = with_hb.cylinders[len(without.cylinders):]
    assert np.all(dashes[:, 6] < 0.1)                        # thinner than bonds


def test_houk_style_scene():
    from microscp.render.styles import make_style

    mol = Molecule(["C", "H", "O"],
                   np.array([[0.0, 0.0, 0.0], [1.05, 0.0, 0.0], [0.0, 1.25, 0.0]]))
    mol.perceive_bonds()
    houk = make_style("houk")
    buf = build_scene(mol, houk)
    # sphere rows carry the seam-line flag: heavy atoms yes, hydrogen no
    assert buf.spheres.shape[1] == 8
    np.testing.assert_allclose(buf.spheres[:, 7], [1.0, 0.0, 1.0])
    # Houkmol bonds are uniformly black instead of split by atom color
    nb = len(mol.bonds)
    expected = np.tile(houk.bond_color, (nb, 1))
    np.testing.assert_allclose(buf.cylinders[:nb, 7:10], expected, atol=1e-6)
    np.testing.assert_allclose(buf.cylinders[:nb, 10:13], expected, atol=1e-6)
    # the default style keeps split-color bonds and no seam lines
    default = build_scene(mol, Style())
    assert np.all(default.spheres[:, 7] == 0.0)
    assert not np.allclose(default.cylinders[0, 7:10], default.cylinders[0, 10:13])


def test_mixed_representations_scene():
    from microscp.render.scene import LINE_RADIUS, REP_BALL, REP_LINE, REP_STICK

    mol = Molecule(["C", "C", "C"],
                   np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]]))
    mol.perceive_bonds()
    style = Style()
    buf = build_scene(mol, style, reps=np.array([REP_BALL, REP_STICK, REP_LINE]))
    # ball keeps the styled atom radius; stick shrinks to the bond radius;
    # line shrinks to the wire radius
    np.testing.assert_allclose(
        buf.spheres[:, 3],
        [style.atom_radius(6), style.bond_radius, LINE_RADIUS], atol=1e-6)
    # a bond is drawn at the thinnest representation of its two atoms
    np.testing.assert_allclose(buf.cylinders[:2, 6],
                               [style.bond_radius, LINE_RADIUS], atol=1e-6)


def test_orca_input_writer(tmp_path):
    mol = mio.load(DATA / "dvb_ir.out").molecule
    out = tmp_path / "job.inp"
    orca.write_inp(out, mol, keywords="wB97X-D3 def2-TZVP Opt")
    text = out.read_text()
    assert text.startswith("! wB97X-D3 def2-TZVP Opt")
    assert f"* xyz {mol.charge} {mol.multiplicity}" in text
    atom_lines = [ln for ln in text.splitlines()
                  if ln.startswith(" ") and len(ln.split()) == 4]
    assert len(atom_lines) == mol.natoms


def test_qchem_input_writer(tmp_path):
    mol = mio.load(DATA / "dvb_ir.out").molecule
    out = tmp_path / "job.in"
    qchem.write_in(out, mol, rem={"jobtype": "freq"})
    text = out.read_text()
    assert "$molecule" in text and text.count("$end") == 2
    assert f"{mol.charge} {mol.multiplicity}" in text
    assert "JOBTYPE" in text and "freq" in text


def test_save_dispatch_by_extension(tmp_path):
    mol = mio.load(DATA / "water_dimer.xyz").molecule
    for name, marker in (("a.inp", "* xyz"), ("a.in", "$molecule"),
                         ("a.qcin", "$molecule")):
        path = tmp_path / name
        mio.save_molecule(path, mol)
        assert marker in path.read_text()
