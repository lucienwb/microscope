"""Manipulator (gizmo) geometry, drag math and the rigid-transform editors.

All headless: hit-testing and drag amounts are pure screen-space geometry on
top of the orthographic camera, so no GL context or QApplication is needed.
"""

import numpy as np
import pytest

from microscope.core import editing
from microscope.core.molecule import Molecule
from microscope.gui import gizmo
from microscope.render.camera import OrthoCamera

W, H = 800.0, 600.0


def _camera():
    cam = OrthoCamera()
    cam.fit(np.zeros(3), 5.0)
    return cam


def _molecule():
    # ethane-like: two carbons plus one hydrogen on each, so bonds are perceived
    mol = Molecule(["C", "C", "H", "H"],
                   np.array([[0.0, 0.0, 0.0],
                             [1.5, 0.0, 0.0],
                             [-1.0, 0.2, 0.0],
                             [2.5, 0.2, 0.0]]))
    mol.perceive_bonds()
    return mol


# ------------------------------------------------------------------ transforms

def test_translate_atoms_moves_only_the_selection():
    mol = _molecule()
    out = editing.translate_atoms(mol, [1, 3], [0.0, 1.0, 0.0])
    assert np.allclose(out[0], mol.coords[0])
    assert np.allclose(out[1], mol.coords[1] + [0, 1, 0])
    assert np.allclose(out[3], mol.coords[3] + [0, 1, 0])
    assert np.allclose(mol.coords, _molecule().coords)   # input untouched


def test_rotate_atoms_is_rigid_about_the_selection_centroid():
    mol = _molecule()
    sel = [1, 3]
    before = np.linalg.norm(mol.coords[1] - mol.coords[3])
    centroid = mol.coords[sel].mean(axis=0)
    out = editing.rotate_atoms(mol, sel, [0, 0, 1], 90.0)
    assert np.isclose(np.linalg.norm(out[1] - out[3]), before)
    assert np.allclose(out[sel].mean(axis=0), centroid)
    assert np.allclose(out[0], mol.coords[0])            # rest stays put


def test_rotate_atoms_about_a_given_pivot():
    mol = _molecule()
    out = editing.rotate_atoms(mol, [1], [0, 0, 1], 90.0, pivot=np.zeros(3))
    assert np.allclose(out[1], [0.0, 1.5, 0.0], atol=1e-9)


def test_rotate_atoms_rejects_a_degenerate_axis():
    with pytest.raises(editing.EditError):
        editing.rotate_atoms(_molecule(), [1], [0, 0, 0], 30.0)


def test_translate_atoms_rejects_an_empty_selection():
    with pytest.raises(editing.EditError):
        editing.translate_atoms(_molecule(), [], [1, 0, 0])


def test_connected_fragment_grows_to_the_whole_bonded_group():
    mol = _molecule()
    assert editing.connected_fragment(mol.bonds, [2]) == [0, 1, 2, 3]
    far = Molecule(["C", "C"], np.array([[0.0, 0, 0], [8.0, 0, 0]]))
    far.perceive_bonds()
    assert editing.connected_fragment(far.bonds, [0]) == [0]


# ------------------------------------------------------------------ hit testing

def test_center_and_axis_handles_are_picked():
    cam = _camera()
    center = np.zeros(3)
    center2d = cam.project(center, W, H)[0]
    assert gizmo.hit_test(cam, W, H, center, center2d) == gizmo.CENTER

    length = gizmo.world_length(cam, int(H))
    for k in range(3):
        tip = cam.project(center + gizmo.AXIS_VECTORS[k] * length, W, H)[0]
        if np.linalg.norm(tip - center2d) < gizmo.MIN_AXIS_PX:
            continue                       # z points at the viewer by default
        mid = (center2d + tip) / 2.0
        assert gizmo.hit_test(cam, W, H, center, mid) == ("axis", k)


def test_ring_handles_are_picked_on_the_circle():
    cam = _camera()
    center = np.zeros(3)
    length = gizmo.world_length(cam, int(H))
    radius = length * gizmo.RING_RADIUS_FRAC
    # the z ring lies in the screen plane when the camera is unrotated; sample
    # it at 45°, away from the arrows (which deliberately win where they cross)
    ring = cam.project(gizmo.ring_points(center, 2, radius), W, H)
    assert gizmo.hit_test(cam, W, H, center, ring[len(ring) // 8]) == ("ring", 2)
    assert gizmo.hit_test(cam, W, H, center, ring[0]) == ("axis", 0)


def test_empty_space_hits_nothing():
    cam = _camera()
    assert gizmo.hit_test(cam, W, H, np.zeros(3), (5.0, 5.0)) is None


def test_edge_on_axis_is_not_pickable():
    """The z axis points straight at the viewer, so its arrow collapses onto
    the centre; it must not swallow clicks meant for something else."""
    cam = _camera()
    center = np.zeros(3)
    center2d = cam.project(center, W, H)[0]
    hit = gizmo.hit_test(cam, W, H, center, center2d + np.array([10.0, 10.0]))
    assert hit != ("axis", 2)


# ------------------------------------------------------------------ drag math

def test_axis_translation_matches_the_projected_drag():
    cam = _camera()
    center = np.zeros(3)
    per_unit = (cam.project(center + gizmo.AXIS_VECTORS[0], W, H)[0]
                - cam.project(center, W, H)[0])
    # dragging exactly one world unit worth of pixels along x gives 1.0 Å
    t = gizmo.axis_translation(cam, W, H, center, 0, per_unit)
    assert np.isclose(t, 1.0)
    # movement perpendicular to the axis contributes nothing
    perp = np.array([-per_unit[1], per_unit[0]])
    assert np.isclose(gizmo.axis_translation(cam, W, H, center, 0, perp), 0.0,
                      atol=1e-9)


def test_axis_translation_follows_the_screen_direction_of_x():
    cam = _camera()
    t = gizmo.axis_translation(cam, W, H, np.zeros(3), 0, (10.0, 0.0))
    assert t > 0.0                     # world +x projects to screen right


def test_plane_translation_uses_camera_right_and_up():
    cam = _camera()
    delta = gizmo.plane_translation(cam, H, (10.0, 0.0))
    assert delta[0] > 0.0 and np.isclose(delta[1], 0.0) and np.isclose(delta[2], 0.0)
    down = gizmo.plane_translation(cam, H, (0.0, 10.0))
    assert down[1] < 0.0               # screen y grows downward


def test_screen_angle_and_rotation_sign():
    cam = _camera()
    center = np.zeros(3)
    center2d = cam.project(center, W, H)[0]
    right = gizmo.screen_angle(cam, W, H, center, center2d + [30.0, 0.0])
    up = gizmo.screen_angle(cam, W, H, center, center2d + [0.0, -30.0])
    assert np.isclose(right, 0.0, atol=1e-9)
    assert np.isclose(up, np.pi / 2, atol=1e-9)
    # z faces the viewer, so its ring turns counter-clockwise on screen
    assert gizmo.rotation_sign(cam, 2) == 1.0


def test_rotation_sign_flips_when_the_axis_turns_away():
    cam = _camera()
    cam.rotation = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
    assert gizmo.rotation_sign(cam, 2) == 1.0
    cam.rotation = np.array([[1.0, 0, 0], [0, -1.0, 0], [0, 0, -1.0]])
    assert gizmo.rotation_sign(cam, 2) == -1.0


def test_wrap_angle_folds_past_half_a_turn():
    assert np.isclose(gizmo.wrap_angle(0.1), 0.1)
    assert np.isclose(gizmo.wrap_angle(2 * np.pi - 0.1), -0.1)
    assert np.isclose(gizmo.wrap_angle(-2 * np.pi + 0.1), 0.1)


def test_gizmo_keeps_a_constant_screen_size():
    cam = _camera()
    near = gizmo.world_length(cam, int(H))
    cam.zoom(5)
    far = gizmo.world_length(cam, int(H))
    assert far < near                  # zoomed in: fewer angstroms per pixel
    for length, camera in ((near, _camera()), (far, cam)):
        tip = camera.project(gizmo.AXIS_VECTORS[0] * length, W, H)[0]
        center2d = camera.project(np.zeros(3), W, H)[0]
        assert np.isclose(np.linalg.norm(tip - center2d), gizmo.AXIS_LENGTH_PX)


def test_selection_center_is_the_centroid():
    mol = _molecule()
    assert gizmo.selection_center(mol.coords, [0, 1]) == pytest.approx(
        [0.75, 0.0, 0.0])
    assert gizmo.selection_center(mol.coords, []) is None
