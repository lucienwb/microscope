"""The flat layout: what is drawn where, through approximate font metrics.

Nothing here builds a QGuiApplication, so the drawing is exercised with
default_measure rather than a painter's real metrics."""

import numpy as np
import pytest
from molecules import H, W, _benzene, _camera, _drawing, _formaldehyde, _molecule, _nitrate, _water

from microscope.core import lewis
from microscope.io import molfile
from microscope.render import lewis2d


def test_hydrogens_fold_into_their_neighbour_by_default():
    structure = lewis.perceive(_water())
    plan = lewis2d.layout(structure, _camera(_water()), W, H)
    assert [a.index for a in plan.atoms] == [0]
    assert plan.atoms[0].text == "OH2"
    assert plan.atoms[0].runs[-1] == ("2", True)            # drawn as a subscript


def test_show_hydrogens_draws_them_all():
    mol = _water()
    structure = lewis.perceive(mol)
    plan = lewis2d.layout(structure, _camera(mol), W, H,
                          lewis2d.LewisOptions(show_hydrogens=True))
    assert len(plan.atoms) == 3
    assert [a.text for a in plan.atoms] == ["O", "H", "H"]


def test_carbon_is_a_bare_vertex_unless_asked_for():
    mol = _benzene()
    structure = lewis.perceive(mol)
    plan = lewis2d.layout(structure, _camera(mol), W, H)
    assert not any(atom.labelled for atom in plan.atoms)
    labelled = lewis2d.layout(structure, _camera(mol), W, H,
                              lewis2d.LewisOptions(carbon_labels=True))
    assert all(atom.text.startswith("C") for atom in labelled.atoms)


def test_a_carbon_with_nothing_drawn_beside_it_is_still_labelled():
    methane = _molecule(["C", "H", "H", "H", "H"],
                        [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
                         [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]])
    plan = lewis2d.layout(lewis.perceive(methane), _camera(methane), W, H)
    assert [a.text for a in plan.atoms] == ["CH4"]


def test_a_bridging_hydrogen_stays_a_vertex_of_its_own():
    # diborane: the two bridging H cannot be folded into a label, and must not
    # be counted into one either (that would draw BH3 twice over)
    diborane = _molecule(
        ["B", "B", "H", "H", "H", "H", "H", "H"],
        [[0, 0, 0], [1.77, 0, 0], [0.885, 0.97, 0], [0.885, -0.97, 0],
         [-0.6, 0, 1.0], [-0.6, 0, -1.0], [2.37, 0, 1.0], [2.37, 0, -1.0]])
    plan = lewis2d.layout(lewis.perceive(diborane), _camera(diborane), W, H)
    assert [a.text for a in plan.atoms] == ["BH2", "BH2", "H", "H"]


def test_hydrogens_lean_away_from_the_bond():
    mol = _molecule(["O", "H", "C", "H", "H", "H"],
                    [[0, 0, 0], [-0.6, 0.8, 0], [1.43, 0, 0],
                     [1.8, 1.03, 0], [1.8, -0.5, 0.9], [1.8, -0.5, -0.9]])
    structure = lewis.perceive(mol)
    plan = lewis2d.layout(structure, _camera(mol), W, H)
    oxygen = next(a for a in plan.atoms if a.index == 0)
    assert oxygen.text == "HO"          # the bond leaves to the right


def test_multiple_bonds_are_drawn_as_parallel_lines():
    mol = _formaldehyde()
    plan = lewis2d.layout(lewis.perceive(mol), _camera(mol), W, H)
    double = next(b for b in plan.bonds if b.order == 2)
    assert len(double.lines) == 2
    acetylene = _molecule(["C", "C", "H", "H"],
                          [[0, 0, 0], [1.20, 0, 0], [-1.06, 0, 0], [2.26, 0, 0]])
    triple = lewis2d.layout(lewis.perceive(acetylene), _camera(acetylene), W, H)
    assert len(triple.bonds[0].lines) == 3


def test_a_ring_double_bond_puts_its_second_line_inside():
    mol = _benzene()
    plan = lewis2d.layout(lewis.perceive(mol), _camera(mol), W, H)
    centre = np.mean([a.pos for a in plan.atoms], axis=0)
    double = next(b for b in plan.bonds if b.order == 2)
    (first, _), (second, _) = double.lines
    assert np.linalg.norm(second - centre) < np.linalg.norm(first - centre)


def test_a_terminal_double_bond_stays_symmetric():
    mol = _formaldehyde()
    plan = lewis2d.layout(lewis.perceive(mol), _camera(mol), W, H)
    double = next(b for b in plan.bonds if b.order == 2)
    axis = plan.atoms[1].pos - plan.atoms[0].pos       # C=O direction
    axis = axis / np.linalg.norm(axis)
    perp = np.array([-axis[1], axis[0]])
    offsets = [float(np.dot(start - plan.atoms[0].pos, perp))
               for start, _ in double.lines]
    assert offsets[0] == pytest.approx(-offsets[1], abs=1e-6)


def test_bonds_stop_short_of_a_label():
    mol = _water()
    structure = lewis.perceive(mol)
    plan = lewis2d.layout(structure, _camera(mol), W, H,
                          lewis2d.LewisOptions(show_hydrogens=True))
    oxygen = next(a for a in plan.atoms if a.index == 0)
    start, _ = plan.bonds[0].lines[0]
    assert np.linalg.norm(start - oxygen.pos) > oxygen.half_width * 0.9


def test_lettering_does_not_shrink_when_the_structure_turns():
    mol = _benzene()
    structure = lewis.perceive(mol)
    camera = _camera(mol)
    face_on = lewis2d.layout(structure, camera, W, H)
    camera.rotate_drag(300, 0)                     # swing the ring edge-on
    edge_on = lewis2d.layout(structure, camera, W, H)
    assert edge_on.font_px == pytest.approx(face_on.font_px)
    assert edge_on.line_width == pytest.approx(face_on.line_width)


def test_stacked_atoms_are_reported_so_the_view_can_be_turned():
    mol = _benzene()
    structure = lewis.perceive(mol)
    camera = _camera(mol)
    assert not lewis2d.layout(structure, camera, W, H).warnings
    camera.rotation = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    edge_on = lewis2d.layout(structure, camera, W, H)
    assert edge_on.overlaps and "turn the structure" in edge_on.warnings[0]


def test_a_disputed_charge_is_carried_into_the_drawing():
    salt = _molecule(["Mo", "Cl"], [[0.0, 0.0, 0.0], [2.40, 0.0, 0.0]], charge=-2)
    plan = lewis2d.layout(lewis.perceive(salt), _camera(salt), W, H)
    assert any("does not match" in w for w in plan.warnings)


def test_a_metal_says_its_charges_are_only_bookkeeping():
    # true whatever the charge works out to: the note has to appear even for a
    # format that never declares one, where there is no charge to disagree with
    salt = _molecule(["Mo", "Cl"], [[0.0, 0.0, 0.0], [2.40, 0.0, 0.0]])
    salt.charge_known = False
    plan = lewis2d.layout(lewis.perceive(salt), _camera(salt), W, H)
    assert any("bookkeeping" in w for w in plan.warnings)
    assert not any("does not match" in w for w in plan.warnings)


def test_the_brackets_enclose_the_whole_drawing():
    ion = _benzene()
    ion.charge, ion.multiplicity = 1, 2
    plan = lewis2d.layout(lewis.perceive(ion), _camera(ion), W, H)
    box = plan.brackets
    assert box is not None
    for atom in plan.atoms:
        assert box.left < atom.pos[0] < box.right
        assert box.top < atom.pos[1] < box.bottom
    assert box.label_pos[0] > box.right          # the charge hangs outside
    assert "".join(t for t, _ in box.runs) == "•+"


def test_a_neutral_molecule_gets_no_brackets():
    assert lewis2d.layout(lewis.perceive(_benzene()), _camera(_benzene()),
                          W, H).brackets is None


def test_a_molfile_records_a_charge_it_cannot_place(tmp_path):
    # V2000 has per-atom charges only, so the species charge goes in the title
    ion = _benzene()
    ion.charge, ion.multiplicity = 1, 2
    structure, points, include = _drawing(ion)
    out = tmp_path / "cation.mol"
    molfile.write_2d(out, structure, points, include)
    first = out.read_text().split("\n")[0]
    assert "charge +1" in first and "unpaired electron" in first


def test_charges_and_lone_pairs_are_placed_only_when_wanted():
    mol = _nitrate()
    structure = lewis.perceive(mol)
    plain = lewis2d.layout(structure, _camera(mol), W, H)
    assert not any(a.dots for a in plain.atoms)
    assert [a.charge_text for a in plain.atoms if a.charge_text].count("−") == 2
    dotted = lewis2d.layout(structure, _camera(mol), W, H,
                            lewis2d.LewisOptions(lone_pairs=True))
    oxygen = next(a for a in dotted.atoms if a.index == 1)
    assert len(oxygen.dots) == 2 * oxygen.lone_pairs
