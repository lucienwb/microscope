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
