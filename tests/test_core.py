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
