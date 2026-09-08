"""Lewis perception: bond orders, charges, lone pairs, radicals.

Checked against textbook answers - those are the regressions that
matter, not the pixels."""

import numpy as np
import pytest
from molecules import _acene, _benzene, _dicyanobenzene, _molecule, _nitrate, _orders, _water

from microscope.core import lewis


def test_benzene_comes_out_as_a_kekule_structure():
    # the trap: plain greedy matching strands two ring atoms with no partner
    structure = lewis.perceive(_benzene())
    assert _orders(structure)[2] == 3
    assert not structure.charges.any()
    assert structure.matches_file


@pytest.mark.parametrize("rings", [2, 3, 4, 5, 6])
def test_fused_rings_are_kekulized_all_the_way_along(rings):
    # ordering the choices alone carries the small ones, but the longer rows
    # leave greed stranding two carbons in the middle with no partner left
    structure = lewis.perceive(_acene(rings))
    carbons = 4 * rings + 2
    assert _orders(structure)[2] == carbons // 2
    assert not structure.charges.any()
    assert structure.matches_file


def test_the_adjacency_is_two_arrays_not_a_list_per_atom():
    structure = lewis.perceive(_benzene())
    assert list(structure.neighbors[0]) == [1, 5, 6]
    assert structure.neighbors.flat.dtype == np.int32
    assert structure.neighbors.nbytes < 200          # 12 atoms, 12 bonds
    assert structure.orders.dtype == np.int8
    assert structure.charges.dtype == np.int8


def test_carbon_dioxide_is_two_double_bonds_not_a_triple():
    structure = lewis.perceive(
        _molecule(["O", "C", "O"], [[-1.16, 0, 0], [0, 0, 0], [1.16, 0, 0]]))
    assert _orders(structure)[2] == 2
    assert not structure.charges.any()


@pytest.mark.parametrize("symbols, coords", [
    (["N", "N"], [[0, 0, 0], [1.098, 0, 0]]),
    (["C", "C", "H", "H"], [[0, 0, 0], [1.20, 0, 0], [-1.06, 0, 0], [2.26, 0, 0]]),
])
def test_short_bonds_become_triple(symbols, coords):
    assert _orders(lewis.perceive(_molecule(symbols, coords)))[3] == 1


def test_water_keeps_two_lone_pairs_and_no_charge():
    structure = lewis.perceive(_water())
    assert list(structure.lone_pairs) == [2, 0, 0]
    assert not structure.charges.any()
    assert list(structure.hydrogens) == [2, 0, 0]


def test_nitrate_is_drawn_charge_separated():
    structure = lewis.perceive(_nitrate())
    assert _orders(structure)[2] == 1          # nitrogen cannot take a second
    assert structure.charges[0] == 1
    assert sorted(structure.charges[1:]) == [-1, -1, 0]
    assert structure.total_charge == -1
    assert structure.matches_file


def test_sulfur_may_expand_its_octet():
    corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    corners = corners / np.linalg.norm(corners, axis=1)[:, None] * 1.47
    structure = lewis.perceive(
        _molecule(["S"] + ["O"] * 4, np.vstack([[0, 0, 0], corners]), charge=-2))
    assert _orders(structure)[2] == 2
    assert structure.total_charge == -2 and structure.matches_file


def test_carbon_monoxide_carries_the_textbook_charges():
    structure = lewis.perceive(_molecule(["C", "O"], [[0, 0, 0], [1.128, 0, 0]]))
    assert _orders(structure)[3] == 1
    assert list(structure.charges) == [-1, 1]


def test_a_doublet_becomes_a_radical_rather_than_an_anion():
    methyl = _molecule(["C", "H", "H", "H"],
                       [[0, 0, 0], [1.08, 0, 0], [-0.54, 0.94, 0],
                        [-0.54, -0.94, 0]], multiplicity=2)
    structure = lewis.perceive(methyl)
    assert structure.radicals[0] == 1
    assert not structure.charges.any()
    assert structure.matches_file


def test_a_bare_ion_takes_the_charge_the_file_declares():
    # geometry cannot tell a chloride from a chlorine atom, but the file can,
    # and there is only one atom the charge could possibly be on
    lonely = _molecule(["Cl"], [[0.0, 0.0, 0.0]], charge=-1)
    structure = lewis.perceive(lonely)
    assert int(structure.charges[0]) == -1
    assert int(structure.lone_pairs[0]) == 4
    assert not structure.radicals.any()
    assert structure.matches_file and not structure.bracketed


def test_a_charge_with_nowhere_to_sit_goes_on_the_brackets():
    # a delocalized radical cation: no atom carries the charge, so the whole
    # species does — which is how a chemist draws it too
    ion = _benzene()
    ion.charge, ion.multiplicity = 1, 2
    structure = lewis.perceive(ion)
    assert structure.bracketed
    assert structure.net_charge == 1 and structure.net_radicals == 1
    assert not structure.charges.any()
    assert structure.matches_file


def test_a_format_that_states_no_charge_is_never_contradicted():
    mol = _benzene()
    mol.charge_known = False
    mol.charge = 7                      # a default nobody should believe
    assert lewis.perceive(mol).matches_file


def test_a_perception_that_disagrees_with_the_file_says_so():
    # a metal complex: the charge cannot be reconciled and is not guessed at
    salt = _molecule(["Mo", "Cl"], [[0.0, 0.0, 0.0], [2.40, 0.0, 0.0]], charge=-2)
    structure = lewis.perceive(salt)
    assert structure.total_charge == 0
    assert not structure.matches_file


def test_a_trigonal_boronate_ester_is_neutral():
    # B-O runs 1.36 (trigonal ester) to 1.50 (borate); calling the short end a
    # double bond turned every pinacol boronate into a borate anion
    angles = np.arange(3) * 2.0 * np.pi / 3.0
    oxygens = np.column_stack([1.37 * np.cos(angles), 1.37 * np.sin(angles),
                               np.zeros(3)])
    mol = _molecule(["B", "O", "O", "O"] + ["H"] * 3,
                    np.vstack([[0, 0, 0], oxygens, oxygens * (1 + 0.96 / 1.37)]))
    structure = lewis.perceive(mol)
    assert _orders(structure)[2] == 0
    assert not structure.charges.any()
    assert structure.matches_file


def test_a_four_coordinate_boron_is_still_an_anion():
    mol = _molecule(["B", "H", "H", "H", "N", "H", "H", "H"],
                    [[0, 0, 0], [0.68, 0.68, -0.48], [-0.93, 0.15, -0.48],
                     [0.25, -0.83, -0.48], [0, 0, 1.60],
                     [0.94, 0, 1.93], [-0.47, 0.81, 1.93], [-0.47, -0.81, 1.93]])
    structure = lewis.perceive(mol)
    assert int(structure.charges[0]) == -1        # B
    assert int(structure.charges[4]) == +1        # N
    assert structure.matches_file


def test_a_metal_carbonyl_keeps_its_triple_bond():
    # back-donation stretches a coordinated CO to 1.13-1.18; at 1.17 it used to
    # read as a ketone and leave a carbanion on the carbon
    mol = _molecule(["Mn", "C", "O"],
                    [[0.0, 0, 0], [1.79, 0, 0], [2.96, 0, 0]])
    structure = lewis.perceive(mol)
    assert structure.order_of(1, 2) == 3
    assert int(structure.charges[1]) == 0
    assert int(structure.charges[2]) == +1        # the metal takes no charge


def test_a_file_without_hydrogens_claims_no_charges():
    # a crystal structure rarely resolves hydrogens, and counting the bonds
    # they would have made as formal charges buries a protein in carbanions
    ring = np.arange(6) * np.pi / 3.0
    coords = np.column_stack([1.54 * np.cos(ring), 1.54 * np.sin(ring), np.zeros(6)])
    structure = lewis.perceive(_molecule(["C"] * 6, coords))
    assert structure.hydrogens_missing
    assert not structure.charges.any()
    assert not structure.lone_pairs.any()


def test_a_molecule_that_really_has_no_hydrogens_keeps_its_charges():
    # nitrate has none either, and its charge separation is the whole point
    structure = lewis.perceive(_nitrate())
    assert not structure.hydrogens_missing
    assert structure.total_charge == -1


def test_metals_are_left_out_of_the_bookkeeping():
    complexed = _molecule(["Mo", "Cl"], [[0, 0, 0], [2.4, 0, 0]])
    structure = lewis.perceive(complexed)
    assert structure.charges[0] == 0 and structure.lone_pairs[0] == 0
    assert int(structure.orders[0]) == 1

def test_a_nitrile_keeps_its_triple_next_to_a_ring():
    """The double pass runs first and can spend the carbon's last valence.

    Terephthalonitrile, 4-cyanopyridine, tetracyanoethylene and a real
    fullerene dye all came out with C=N and a nitrogen anion.
    """
    structure = lewis.perceive(_dicyanobenzene())
    assert _orders(structure)[3] == 2          # both C#N
    assert _orders(structure)[2] == 3          # and a Kekule ring
    assert not structure.charges.any()
    assert structure.matches_file


def test_carbon_dioxide_still_refuses_the_triple():
    """The counter-case: its carbon wants two triples and can afford one, so
    neither may claim the valence up front."""
    structure = lewis.perceive(
        _molecule(["O", "C", "O"], [[-1.16, 0, 0], [0, 0, 0], [1.16, 0, 0]]))
    assert _orders(structure)[2] == 2 and _orders(structure)[3] == 0
