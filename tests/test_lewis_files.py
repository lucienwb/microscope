"""The ChemDraw writers and the --lewis command line."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from molecules import _benzene, _drawing, _formaldehyde, _molecule, _nitrate, _water

from microscope import cli
from microscope.io import cdxml, molfile
from microscope.io.errors import UnsupportedFormatError


def test_cdxml_is_well_formed_and_keeps_the_chemistry(tmp_path):
    import xml.etree.ElementTree as ET

    structure, points, include = _drawing(_nitrate())
    path = tmp_path / "nitrate.cdxml"
    cdxml.write(path, structure, points, include)
    root = ET.parse(path).getroot()
    assert root.tag == "CDXML"
    nodes = root.findall(".//n")
    assert len(nodes) == 4
    assert {n.get("Element") for n in nodes} == {"7", "8"}
    charges = sorted(n.get("Charge") or "0" for n in nodes)
    assert charges == ["-1", "-1", "0", "1"]
    orders = sorted(b.get("Order") for b in root.findall(".//b"))
    assert orders == ["1", "1", "2"]


def test_cdxml_folds_hydrogens_into_a_count(tmp_path):
    import xml.etree.ElementTree as ET

    structure, points, include = _drawing(_water())
    path = tmp_path / "water.cdxml"
    cdxml.write(path, structure, points, include)
    nodes = ET.parse(path).getroot().findall(".//n")
    assert len(nodes) == 1 and nodes[0].get("NumHydrogens") == "2"


def test_cdxml_uses_chemdraws_bond_length(tmp_path):
    structure, points, include = _drawing(_benzene())
    path = tmp_path / "benzene.cdxml"
    cdxml.write(path, structure, points, include)
    positions = {}
    for line in path.read_text().splitlines():
        if "<n " in line:
            key = line.split('id="')[1].split('"')[0]
            positions[key] = [float(v) for v in
                              line.split('p="')[1].split('"')[0].split()]
    bonds = [(1, 2), (2, 3), (3, 4)]
    lengths = [np.linalg.norm(np.array(positions[str(i)])
                              - np.array(positions[str(j)])) for i, j in bonds]
    assert np.median(lengths) == pytest.approx(cdxml.BOND_LENGTH, rel=0.05)


def test_molfile_columns_follow_the_v2000_layout(tmp_path):
    structure, points, include = _drawing(_formaldehyde())
    path = tmp_path / "ch2o.mol"
    molfile.write_2d(path, structure, points, include)
    lines = path.read_text().splitlines()
    counts = lines[3]
    assert counts[:3].strip() == "2"          # C and O; both hydrogens folded in
    assert counts[3:6].strip() == "1"
    assert counts.endswith("V2000")
    atom = lines[4]
    assert len(atom) == 69
    assert atom[31:34].strip() in ("C", "O")
    assert atom[48:51].strip() in ("3", "4")  # valence field: bonds + implicit H
    assert lines[6].startswith("  1  2  2")   # the double bond


def test_molfile_writes_charges_and_radicals(tmp_path):
    structure, points, include = _drawing(_nitrate())
    path = tmp_path / "nitrate.mol"
    molfile.write_2d(path, structure, points, include)
    text = path.read_text()
    assert "M  CHG  3" in text
    assert text.rstrip().endswith("M  END")

    methyl = _molecule(["C", "H", "H", "H"],
                       [[0, 0, 0], [1.08, 0, 0], [-0.54, 0.94, 0],
                        [-0.54, -0.94, 0]], multiplicity=2)
    radical = tmp_path / "methyl.mol"
    molfile.write_2d(radical, *_drawing(methyl))
    assert "M  RAD  1   1   2" in radical.read_text()   # one doublet on atom 1


def test_molfile_y_axis_points_up(tmp_path):
    # screen y grows downwards, molfiles the other way
    mol = _molecule(["C", "O", "H", "H"],
                    [[0, 0, 0], [0, 1.21, 0], [-0.94, -0.57, 0],
                     [0.94, -0.57, 0]])
    structure, points, include = _drawing(mol)
    path = tmp_path / "up.mol"
    molfile.write_2d(path, structure, points, include)
    lines = path.read_text().splitlines()
    carbon_y, oxygen_y = (float(lines[4 + k][10:20]) for k in (0, 1))
    assert (oxygen_y > carbon_y) == (points[1][1] < points[0][1])


def test_three_dimensional_molfile_keeps_the_geometry(tmp_path):
    import microscope.io as mio

    mol = _water()
    path = tmp_path / "water.mol"
    mio.save_molecule(path, mol)
    lines = path.read_text().splitlines()
    assert lines[3][:3].strip() == "3"
    written = np.array([[float(line[i:i + 10]) for i in (0, 10, 20)]
                        for line in lines[4:7]])
    assert written == pytest.approx(mol.coords, abs=1e-4)


def test_save_lewis_dispatches_on_the_extension(tmp_path):
    import microscope.io as mio

    structure, points, include = _drawing(_water())
    for name in ("a.cdxml", "b.mol", "c.sdf"):
        mio.save_lewis(tmp_path / name, structure, points, include)
        assert (tmp_path / name).stat().st_size > 0
    with pytest.raises(UnsupportedFormatError):
        mio.save_lewis(tmp_path / "d.png", structure, points, include)


# ------------------------------------------------------------------ command line

def test_lewis_output_accepts_pictures_and_structures():
    for good in ("f.png", "f.tif", "f.svg", "f.pdf", "f.cdxml", "f.mol"):
        assert cli.check_output(good, "lewis") == good
    with pytest.raises(cli.CliError):
        cli.check_output("f.jpg", "lewis")
    with pytest.raises(cli.CliError):
        cli.check_output("f.eps", "lewis")       # that one is for spectra


def test_lewis_default_output_and_size():
    args = cli.build_parser().parse_args(["-s", "mol.log", "--lewis"])
    assert cli.default_output("mol.log", "lewis") == "mol_lewis.png"
    assert cli.resolve_size(args, "lewis") == (1200, 900)   # not the wide frame


def test_lewis_option_flags_reach_the_drawing():
    args = cli.build_parser().parse_args(
        ["-s", "mol.log", "--lewis", "--lone-pairs", "--show-hydrogens"])
    assert cli.lewis_option_dict(args) == {
        "show_hydrogens": True, "carbon_labels": False,
        "lone_pairs": True, "color_atoms": False}
    options = cli.build_lewis_options(args)
    assert options.lone_pairs and options.show_hydrogens
    assert not options.color_atoms


@pytest.mark.parametrize("argv", [
    ["-s", "mol.log", "--lewis", "--style", "houk"],
    ["-s", "mol.log", "--lewis", "--rep", "stick"],
    ["-s", "mol.log", "--lewis", "--measure", "1,2"],
    ["-s", "mol.log", "--lewis", "--iso", "0.02"],
])
def test_renderer_flags_are_refused_for_a_drawing(argv):
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        cli.check_lewis_flags(parser.parse_args(argv), parser)


@pytest.mark.parametrize("flag", ["--lone-pairs", "--show-hydrogens",
                                  "--carbon-labels", "--color-atoms"])
def test_drawing_options_need_lewis(flag):
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        cli.check_lewis_flags(parser.parse_args(["-s", "mol.log", flag]), parser)


def test_lewis_and_a_spectrum_are_different_modes():
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        cli.check_mode_flags(
            parser.parse_args(["-s", "f.log", "--ir", "--lewis"]), parser)


def test_spin_turns_the_picture_in_the_page():
    mol = _benzene()
    camera = cli.build_camera(mol, None, None, None, 1.0, spin=90.0)
    assert camera.rotation @ np.array([1.0, 0.0, 0.0]) == pytest.approx([0, 1, 0])


def test_writing_a_chemdraw_file_needs_no_graphics(tmp_path):
    """--lewis -o x.cdxml is pure text, so it must not so much as import Qt.

    Run in a child process and asked afterwards what it loaded: a test that
    only checked our own "start Qt" helper would keep passing if the import
    moved somewhere else.
    """
    out = tmp_path / "trp.cdxml"
    script = (
        "import sys; from microscope.cli import main;"
        f"code = main(['-s', 'tests/data/trp.log', '--lewis', '-o', {str(out)!r}]);"
        "loaded = [m for m in sys.modules if m.split('.')[0] in "
        "('PySide6', 'OpenGL', 'matplotlib')];"
        "print(code, loaded)")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
    proc = subprocess.run([sys.executable, "-c", script], capture_output=True,
                          text=True, env=env, cwd=str(Path(__file__).parent.parent))
    assert proc.returncode == 0, proc.stderr
    verdict = proc.stdout.strip().splitlines()[-1]   # main() prints the path first
    assert verdict == "0 []", proc.stdout
    assert "<CDXML" in out.read_text()
