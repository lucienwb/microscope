"""Cube file parsing and isosurface extraction tests."""

from pathlib import Path

import numpy as np
import pytest

import microscope.io as mio
from microscope.core.isosurface import isosurface_mesh, marching_tetrahedra
from microscope.core.volume import VolumeData
from microscope.io import cube
from microscope.io.errors import FileFormatError
from microscope.render.scene import build_surface_meshes
from microscope.render.styles import Style

DATA = Path(__file__).parent / "data"
BOHR = 0.52917721092


def test_cube_read():
    result = mio.load(DATA / "water_mo.cube")
    assert result.program == "cube"
    mol = result.molecule
    assert mol.formula() == "H2O"
    assert len(result.volumes) == 1
    vol = result.volumes[0]
    assert vol.shape == (20, 20, 20)
    assert vol.is_signed
    # -4 Bohr origin converted to Angstrom; last grid point at +4 Bohr
    assert np.allclose(vol.origin, -4.0 * BOHR, atol=1e-6)
    corner = vol.grid_to_world(np.array([[19.0, 19.0, 19.0]]))[0]
    assert np.allclose(corner, 4.0 * BOHR, atol=1e-6)
    # O-H distance survives the Bohr -> Angstrom conversion
    d = np.linalg.norm(mol.coords[1] - mol.coords[0])
    assert abs(d - np.hypot(1.43, 1.11) * BOHR) < 1e-6
    # psi = z * exp(-1.4 r) is antisymmetric in z (k -> 19 - k mirrors z)
    assert np.allclose(vol.values[:, :, 5], -vol.values[:, :, 14], atol=1e-9)


def test_cube_multi_mo():
    text = "\n".join([
        "two MOs",
        "MO coefficients",
        "   -1    0.0 0.0 0.0",
        "    2    1.0 0.0 0.0",
        "    2    0.0 1.0 0.0",
        "    2    0.0 0.0 1.0",
        "    1    1.0    0.0 0.0 0.0",
        "    2    3    7",
        # 8 grid points x 2 MOs, MO7 = -MO3, z fastest
        " 0.0 -0.0  1.0 -1.0  2.0 -2.0",
        " 3.0 -3.0  4.0 -4.0  5.0 -5.0",
        " 6.0 -6.0  7.0 -7.0",
    ])
    path = DATA / "_tmp_multi.cube"
    path.write_text(text)
    try:
        result = cube.read(path)
        assert [v.label for v in result.volumes] == ["MO 3", "MO 7"]
        v3, v7 = result.volumes
        assert v3.shape == (2, 2, 2)
        assert np.allclose(v3.values.ravel(), np.arange(8.0))
        assert np.allclose(v7.values, -v3.values)
    finally:
        path.unlink()


def test_cube_truncated_rejected():
    path = DATA / "_tmp_trunc.cube"
    path.write_text("\n".join([
        "t", "t", " 1 0.0 0.0 0.0", " 2 1.0 0.0 0.0", " 2 0.0 1.0 0.0",
        " 2 0.0 0.0 1.0", " 1 1.0 0.0 0.0 0.0", " 0.0 1.0",
    ]))
    try:
        with pytest.raises(FileFormatError, match="truncated"):
            cube.read(path)
    finally:
        path.unlink()


def _sphere_volume(n=24, extent=1.6):
    """values = -|r| so the region above level -R is the ball of radius R."""
    lin = np.linspace(-extent, extent, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    values = -np.sqrt(x**2 + y**2 + z**2)
    step = 2.0 * extent / (n - 1)
    return VolumeData(origin=np.full(3, -extent), axes=np.eye(3) * step,
                      values=values)


def test_isosurface_sphere():
    vol = _sphere_volume()
    verts, normals = isosurface_mesh(vol, -1.0)
    assert len(verts) > 300 and len(verts) % 3 == 0
    radii = np.linalg.norm(verts, axis=1)
    assert np.all(np.abs(radii - 1.0) < 0.08)          # on the sphere, ~voxel tol
    assert np.allclose(np.linalg.norm(normals, axis=1), 1.0, atol=1e-6)
    outward = np.sum(normals * (verts / radii[:, None]), axis=1)
    assert np.mean(outward) > 0.95                     # radially outward


def test_isosurface_empty_and_degenerate():
    vol = _sphere_volume()
    verts, normals = isosurface_mesh(vol, 10.0)        # above the data range
    assert len(verts) == 0 and len(normals) == 0
    assert marching_tetrahedra(np.zeros((1, 5, 5)), 0.5).shape == (0, 3, 3)


def test_surface_meshes_two_lobes():
    vol = mio.load(DATA / "water_mo.cube").volumes[0]
    meshes = build_surface_meshes(vol, 0.02, Style())
    assert len(meshes) == 2
    for verts, normals, rgba in meshes:
        assert len(verts) > 0 and len(verts) % 3 == 0
        assert verts.shape == normals.shape and len(rgba) == 4
    # psi = z * exp(-1.4 r): + lobe above the molecular plane, - lobe below
    assert np.mean(meshes[0][0][:, 2]) > 0.2
    assert np.mean(meshes[1][0][:, 2]) < -0.2
    # a density-like (all positive) grid produces only one surface
    dens = VolumeData(origin=vol.origin, axes=vol.axes, values=np.abs(vol.values))
    assert len(build_surface_meshes(dens, 0.02, Style())) == 1
