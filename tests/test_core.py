"""Core geometry and molecule tests (no GUI, no GL)."""

import numpy as np

from microscp.core import elements, geometry
from microscp.core.molecule import Molecule


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
    from microscp.render.camera import orientation_along, orientation_from_plane

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
