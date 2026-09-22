"""Core geometry and molecule tests (no GUI, no GL)."""

import numpy as np

from microscope.core import elements, geometry
from microscope.core.molecule import Molecule


def test_normalize_symbol():
    assert elements.normalize_symbol("c") == "C"
    assert elements.normalize_symbol("CL") == "Cl"
    assert elements.normalize_symbol("6") == "C"
    assert elements.normalize_symbol("Fe") == "Fe"


def test_measurements():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.5, 0.0, 0.0])
    c = np.array([1.5, 1.5, 0.0])
    d = np.array([1.5, 1.5, 1.5])
    assert geometry.distance(a, b) == 1.5
    assert abs(geometry.angle(a, b, c) - 90.0) < 1e-9
    assert abs(abs(geometry.dihedral(a, b, c, d)) - 90.0) < 1e-9


def test_orientation_helpers():
    from microscope.render.camera import orientation_along, orientation_from_plane

    rot = orientation_along(np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(rot) > 0.99
    # the direction becomes the view z axis
    np.testing.assert_allclose(rot @ (np.array([1.0, 2.0, 3.0]) / np.sqrt(14)),
                               [0, 0, 1], atol=1e-12)

    a, b, c = np.array([0.0, 0, 0]), np.array([1.5, 0, 0]), np.array([0.7, 1.2, 0])
    rot = orientation_from_plane(a, b, c)
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-12)
    # plane normal maps onto the view z axis -> plane lies in the screen
    normal = np.cross(b - a, c - a)
    normal /= np.linalg.norm(normal)
    np.testing.assert_allclose(rot @ normal, [0, 0, 1], atol=1e-12)
    # degenerate cases
    assert orientation_along(np.zeros(3)) is None
    assert orientation_from_plane(a, a, b) is None


def test_water_bonds():
    mol = Molecule(
        ["O", "H", "H"],
        np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]),
    )
    bonds = mol.perceive_bonds()
    assert len(bonds) == 2
    assert mol.formula() == "H2O"


def test_bonds_of_a_big_structure_match_the_pairwise_answer():
    """The cell grid must find exactly what comparing every pair finds.

    Proteins are tens of thousands of atoms, where the pairwise difference
    array alone would be hundreds of gigabytes, so the grid is the only path
    that runs at all — but it has to give the same bonds.
    """
    rng = np.random.default_rng(4)
    lattice = np.stack(np.meshgrid(*[np.arange(12) * 1.4] * 3), -1).reshape(-1, 3)
    coords = lattice + rng.normal(scale=0.10, size=lattice.shape)
    mol = Molecule(["C"] * len(coords), coords)
    bonds = mol.perceive_bonds()

    radii = np.array([elements.covalent_radius(z) for z in mol.atomic_numbers])
    dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    cutoff = radii[:, None] + radii[None, :] + 0.45
    ii, jj = np.where(np.triu((dist <= cutoff) & (dist > 0.4), k=1))
    assert len(bonds) == len(ii) > 4000
    assert (bonds == np.column_stack([ii, jj])).all()


def test_an_atom_alone_in_its_cell_gets_no_bonds():
    mol = Molecule(["He", "He"], np.array([[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]]))
    assert len(mol.perceive_bonds()) == 0


def test_result_arrays_use_narrow_types():
    """A formal charge is a single digit; the default int is eight bytes.

    At protein scale these are megabytes each, so the narrowing is deliberate
    and worth pinning down.
    """
    mol = Molecule(["O", "H", "H"],
                   np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469],
                             [0.0, -0.757, -0.469]]))
    assert mol.perceive_bonds().dtype == np.int32
    assert mol.atomic_numbers.dtype == np.int16


def test_angle_arc_points():
    a = np.array([1.5, 0.0, 0.0])
    b = np.array([0.0, 0.0, 0.0])
    c = np.array([0.0, 2.0, 0.0])
    arc = geometry.angle_arc_points(a, b, c, radius=0.5, segments=16)
    assert arc.shape == (17, 3)
    # all points on the circle of radius 0.5 around the vertex, in the plane
    np.testing.assert_allclose(np.linalg.norm(arc - b, axis=1), 0.5, atol=1e-12)
    np.testing.assert_allclose(arc[:, 2], 0.0, atol=1e-12)
    # endpoints aligned with the two arms
    np.testing.assert_allclose(arc[0], [0.5, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(arc[-1], [0.0, 0.5, 0.0], atol=1e-12)
    # collinear atoms have no angle plane
    assert geometry.angle_arc_points(a, b, -a).shape == (0, 3)


def test_dihedral_arc_points():
    t = np.radians(60.0)
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, 0.0])
    c = np.array([0.0, 0.0, 2.0])
    d = np.array([np.cos(t), np.sin(t), 2.0])
    arc = geometry.dihedral_arc_points(a, b, c, d, radius=0.4, segments=20)
    assert arc.shape == (21, 3)
    # arc sits at the central-bond midpoint, perpendicular to the bond
    np.testing.assert_allclose(arc[:, 2], 1.0, atol=1e-12)
    np.testing.assert_allclose(
        np.linalg.norm(arc - [0.0, 0.0, 1.0], axis=1), 0.4, atol=1e-12)
    # sweeps from the a side to the d side by the dihedral magnitude
    np.testing.assert_allclose(arc[0], [0.4, 0.0, 1.0], atol=1e-12)
    np.testing.assert_allclose(arc[-1][:2], 0.4 * np.array([np.cos(t), np.sin(t)]),
                               atol=1e-12)
    swept = np.degrees(np.arccos(np.clip(
        np.dot(arc[0] - [0, 0, 1.0], arc[-1] - [0, 0, 1.0]) / 0.16, -1, 1)))
    assert abs(swept - abs(geometry.dihedral(a, b, c, d))) < 1e-9
    # degenerate: outer atom on the bond axis / planar dihedral of 0
    on_axis = np.array([0.0, 0.0, -1.0])
    assert geometry.dihedral_arc_points(on_axis, b, c, d).shape == (0, 3)
    assert geometry.dihedral_arc_points(a, b, c, a + [0, 0, 2.0]).shape == (0, 3)


def test_dihedral_arm_points():
    t = np.radians(60.0)
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, 0.0])
    c = np.array([0.0, 0.0, 2.0])
    d = np.array([np.cos(t), np.sin(t), 2.0])
    arms = geometry.dihedral_arm_points(a, b, c, d)
    assert arms.shape == (2, 2, 3)
    mid = np.array([0.0, 0.0, 1.0])
    # both arms start at the bond midpoint and end at the in-plane
    # projections of the outer atoms
    np.testing.assert_allclose(arms[0][0], mid, atol=1e-12)
    np.testing.assert_allclose(arms[1][0], mid, atol=1e-12)
    np.testing.assert_allclose(arms[0][1], [1.0, 0.0, 1.0], atol=1e-12)
    np.testing.assert_allclose(arms[1][1], [np.cos(t), np.sin(t), 1.0], atol=1e-12)
    # the arc endpoints lie on those arms
    arc = geometry.dihedral_arc_points(a, b, c, d, radius=0.4)
    np.testing.assert_allclose(arc[0], mid + 0.4 * (arms[0][1] - mid), atol=1e-12)
    np.testing.assert_allclose(arc[-1], mid + 0.4 * (arms[1][1] - mid), atol=1e-12)
    # degenerate input
    on_axis = np.array([0.0, 0.0, -1.0])
    assert geometry.dihedral_arm_points(on_axis, b, c, d).shape == (0, 2, 3)


def test_close_pairs_finds_exactly_the_brute_force_pairs():
    from microscope.core.neighbors import close_pairs

    rng = np.random.default_rng(3)
    points = rng.uniform(-8, 8, size=(600, 3))
    diff = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    for cutoff in (0.9, 2.5):
        want = [(a, b) for a, b in zip(*np.nonzero(np.triu(diff <= cutoff, k=1)))]
        i, j, d = close_pairs(points, cutoff)
        assert list(zip(i.tolist(), j.tolist())) == want
        assert np.allclose(d, diff[i, j])
    subset = np.arange(0, 600, 3)
    i, j, _ = close_pairs(points, 2.5, among=subset)
    assert set(i) <= set(subset) and set(j) <= set(subset)
    assert len(i) == np.triu(diff[np.ix_(subset, subset)] <= 2.5, k=1).sum()


def test_hbonds_match_the_every_pair_answer_in_a_box_of_water():
    from microscope.core.contacts import find_hbonds

    rng = np.random.default_rng(11)
    symbols, coords = [], []
    for centre in rng.uniform(0, 14, size=(160, 3)):
        turn = geometry.rotation_matrix(rng.normal(size=3), rng.uniform(0, 180))
        for s, p in (("O", (0, 0, 0)), ("H", (0.96, 0, 0)), ("H", (-0.24, 0.93, 0))):
            symbols.append(s)
            coords.append(centre + turn @ np.array(p))
    water = Molecule(symbols, np.array(coords))
    water.perceive_bonds()
    found = find_hbonds(water)
    # the definition, pair by pair
    bonded = {tuple(b) for b in water.bonds.tolist()} | {tuple(b[::-1]) for b in water.bonds.tolist()}
    want = []
    for h in range(len(symbols)):
        if symbols[h] != "H":
            continue
        donors = sorted(o for o in range(len(symbols)) if (h, o) in bonded and symbols[o] == "O")
        if not donors:
            continue
        for a in range(len(symbols)):
            if symbols[a] != "O" or a == donors[0] or (h, a) in bonded:
                continue
            d = np.linalg.norm(water.coords[a] - water.coords[h])
            if 1.2 <= d <= 2.6 and geometry.angle(water.coords[donors[0]], water.coords[h],
                                                  water.coords[a]) >= 120:
                want.append((h, a))
    assert [(h, a) for h, a, _ in found] == want and len(want) > 10


def test_atomic_numbers_follow_the_symbols_list():
    mol = Molecule(["C", "O"], np.zeros((2, 3)) + [[0, 0, 0], [1.2, 0, 0]])
    first = mol.atomic_numbers
    assert first.tolist() == [6, 8] and mol.atomic_numbers is first     # kept, not rebuilt
    mol.symbols = ["N", "O"]                   # what undo does: a new list
    assert mol.atomic_numbers.tolist() == [7, 8]
    assert not first.flags.writeable           # shared, so nobody may edit it
