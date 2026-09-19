"""Orbitals from wavefunction files: the basis read right, and put on a grid right.

The check that matters is C^T S C = I. A program's orbitals are orthonormal
in its own basis, so they come out that way only if every convention in the
file - component order, normalization, pure or Cartesian, signs - was read
the way the program meant it. The Molden dialects are written here from known
orbitals, the way each program writes them, and must be recognized.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import microscope
from microscope.core import basis as B
from microscope.core.orbitals import evaluate, grid_around
from microscope.io import molden

DATA = Path(__file__).parent / "data"


def _brute_force(basis, coefficients, points):
    """An orbital at some points, straight from the definition."""
    values = np.zeros(len(points))
    for shell, start in zip(basis.shells, basis.offsets()):
        c = np.asarray(coefficients[start:start + shell.nfunc], float) @ shell.transform
        d = points - shell.center
        radial = np.exp(-shell.exponents[None, :] * (d ** 2).sum(1)[:, None]) @ shell.weights
        for (a, b, z), coefficient in zip(B.cartesian_components(shell.l), c):
            values += coefficient * d[:, 0] ** a * d[:, 1] ** b * d[:, 2] ** z * radial
    return values


# -------------------------------------------------------------------- basis

def test_solid_harmonics_have_the_usual_signs():
    assert B.solid_harmonic(1, 1) == {(1, 0, 0): 1.0}
    assert B.solid_harmonic(2, 1) == {(1, 0, 1): 2.0}                    # xz
    assert B.solid_harmonic(2, -2) == {(1, 1, 0): 2.0}                   # xy
    d0 = B.solid_harmonic(2, 0)                                          # 2z^2 - x^2 - y^2
    assert d0[(0, 0, 2)] > 0 and d0[(2, 0, 0)] == d0[(0, 2, 0)] == -d0[(0, 0, 2)] / 2
    assert B.solid_harmonic(3, 3) == {(3, 0, 0): 1.0, (1, 2, 0): -3.0}   # x^3 - 3xy^2
    assert B.solid_harmonic(3, -3) == {(2, 1, 0): 3.0, (0, 3, 0): -1.0}  # 3x^2y - y^3


def test_pure_order_is_zero_then_plus_minus():
    assert B.pure_m_order(2) == [0, 1, -1, 2, -2]


@pytest.mark.parametrize("l", range(6))
@pytest.mark.parametrize("pure", [False, True])
def test_every_basis_function_is_normalized(l, pure):
    rows = (B.pure_rows(l, B.pure_m_order(l)) if pure
            else B.cartesian_rows(l, ["x" * a + "y" * b + "z" * c
                                      for a, b, c in B.cartesian_components(l)]))
    shell = B.make_shell([0.0, 0.1, -0.2], l, [3.1, 0.7, 0.2], [0.2, 0.5, 0.4], rows)
    S = B.overlap_matrix([shell])
    assert np.allclose(np.diag(S), 1.0)
    if pure:        # distinct m on one centre are orthogonal
        assert np.allclose(S, np.eye(len(S)), atol=1e-12)


def test_overlap_of_two_s_gaussians_has_its_closed_form():
    a, b, R = 1.3, 0.4, 1.7
    shells = [B.make_shell([0, 0, 0], 0, [a], [1.0], B.cartesian_rows(0, ["s"])),
              B.make_shell([0, 0, R], 0, [b], [1.0], B.cartesian_rows(0, ["s"]))]
    expected = (2 * np.sqrt(a * b) / (a + b)) ** 1.5 * np.exp(-a * b / (a + b) * R ** 2)
    assert B.overlap_matrix(shells)[0, 1] == pytest.approx(expected, rel=1e-12)


# --------------------------------------------------------------------- fchk

@pytest.mark.parametrize("name", ["dvb_ir.fchk", "dvb_un_sp.fchk", "g16_mo4ocl4.fchk"])
def test_fchk_orbitals_come_out_orthonormal(name):
    """sp shells, an unrestricted wavefunction, and pure d on a metal."""
    orbitals = microscope.load(DATA / name).orbitals
    assert orbitals is not None
    assert orbitals.orthonormality_error(occupied_only=False) < 1e-5


def test_fchk_occupations_follow_the_electron_counts():
    closed = microscope.load(DATA / "dvb_ir.fchk").orbitals
    assert closed.restricted and closed.basis.nbf == 60
    assert closed.alpha.homo == 34 and closed.alpha.occupations[34] == 2
    radical = microscope.load(DATA / "dvb_un_sp.fchk").orbitals      # the cation doublet
    assert not radical.restricted
    assert (radical.alpha.homo, radical.beta.homo) == (34, 33)


def test_a_geometry_only_fchk_has_no_orbitals(tmp_path):
    text = (DATA / "dvb_ir.fchk").read_text()
    cut = text[:text.index("Alpha Orbital Energies")]
    (tmp_path / "geometry.fchk").write_text(cut)
    result = microscope.load(tmp_path / "geometry.fchk")
    assert result.orbitals is None and result.molecule.natoms == 20


# -------------------------------------------------------------------- names

def test_orbitals_are_found_by_name():
    orbitals = microscope.load(DATA / "dvb_ir.fchk").orbitals
    assert orbitals.find("homo") == ("alpha", 34)
    assert orbitals.find("LUMO") == ("alpha", 35)
    assert orbitals.find("homo-2") == ("alpha", 32)
    assert orbitals.find("lumo+1") == ("alpha", 36)
    assert orbitals.find("35") == ("alpha", 34)              # numbered as programs do
    assert orbitals.name("alpha", 32) == "HOMO-2"
    assert orbitals.describe("alpha", 35).startswith("LUMO · MO 36 · ")
    for bad in ("homo+99", "0", "beta:homo", "sideways", "gamma:homo"):
        with pytest.raises(ValueError):
            orbitals.find(bad)


def test_unrestricted_orbitals_name_their_spin():
    orbitals = microscope.load(DATA / "dvb_un_sp.fchk").orbitals
    assert orbitals.find("beta:homo") == ("beta", 33)
    assert orbitals.find("b:lumo") == ("beta", 34)
    assert orbitals.describe("beta", 33).startswith("β HOMO")


# --------------------------------------------------------------------- grid

@pytest.mark.parametrize("name, which", [("dvb_ir.fchk", "homo"),
                                         ("g16_mo4ocl4.fchk", "homo"),
                                         ("g16_mo4ocl4.fchk", "lumo+2")])
def test_grid_values_are_the_orbital(name, which):
    """The fast path (separable factors, screened boxes) against the formula."""
    orbitals = microscope.load(DATA / name).orbitals
    spin, index = orbitals.find(which)
    coefficients = orbitals.spin(spin).coefficients[index]
    origin, step, shape = grid_around(orbitals.centers(), spacing=0.15)
    values = evaluate(orbitals.basis, coefficients, origin, step, shape)
    rng = np.random.default_rng(7)
    picks = np.stack([rng.integers(0, n, 500) for n in shape], axis=1)
    expected = _brute_force(orbitals.basis, coefficients, origin + picks * step)
    got = values[picks[:, 0], picks[:, 1], picks[:, 2]]
    assert np.abs(got - expected).max() < 1e-5
    # and the grid holds all of it: a normalized orbital integrates to one
    assert float((values.astype(float) ** 2).sum() * np.prod(step)) == pytest.approx(1, abs=5e-3)


def test_an_orbital_volume_is_ready_to_draw():
    result = microscope.load(DATA / "dvb_ir.fchk")
    volume = result.orbitals.volume("homo")
    assert volume.values.dtype == np.float32
    assert volume.is_signed and volume.suggest_isovalue() == 0.02
    assert volume.label.startswith("HOMO · MO 35")
    assert volume.meta == {"spin": "alpha", "orbital": 34}
    spacing = np.diag(volume.axes)
    assert np.allclose(np.diag(np.diag(volume.axes)), volume.axes)       # axis-aligned
    assert np.allclose(spacing, 0.12)
    far = volume.origin + spacing * (np.array(volume.shape) - 1)
    coords = result.molecule.coords
    assert (coords.min(0) - volume.origin > 3.9).all() and (far - coords.max(0) > 3.9).all()


def test_a_big_box_gets_a_coarser_grid_not_a_bigger_one():
    centers = np.array([[0, 0, 0], [60, 60, 60]], float)
    _, step, shape = grid_around(centers, max_points=1_000_000)
    assert np.prod(shape) <= 1_050_000 and step[0] > 0.12 / B.BOHR_TO_ANGSTROM


# ------------------------------------------------------------------ Molden

# Two atoms close enough that every function overlaps functions on the other,
# so a wrong sign or scale anywhere spoils orthonormality. Several primitives,
# so the two ways of writing contraction coefficients differ too.
_SHELLS = [(1, "s", [5.0, 1.2, 0.3], [0.3, 0.6, 0.4]),
           (1, "p", [2.0, 0.5], [0.5, 0.6]),
           (1, "d", [1.1, 0.4], [0.6, 0.5]),
           (1, "f", [0.9], [1.0]),
           (2, "s", [3.0, 0.4], [0.4, 0.7]),
           (2, "g", [0.8], [1.0]),
           (2, "d", [0.7], [1.0])]
_CENTRES = {1: np.array([0.0, 0.2, -0.1]), 2: np.array([0.3, -0.4, 2.2])}


def _reference_basis(pure: set[int]):
    raw = [molden._RawShell(a, B.angular_momentum(x), np.array(e), np.array(c))
           for a, x, e, c in _SHELLS]
    return raw, molden._build_basis(raw, _CENTRES, pure, molden.Convention())


def _orthonormal_orbitals(basis, count: int, seed: int = 3) -> np.ndarray:
    """Random orbitals, made orthonormal in the metric of the basis."""
    S = basis.overlap()
    w, v = np.linalg.eigh(S)
    inverse_root = v @ np.diag(w ** -0.5) @ v.T
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(len(S), len(S))))
    return (inverse_root @ q).T[:count]


def _as_written(raw, reference, orbitals, pure, convention):
    """What a program using *convention* writes: its contraction coefficients
    and its orbital coefficients for the same orbitals."""
    coefficients = [np.asarray(s.coefficients) * (
        1.0 if convention.normalized_primitives else B.primitive_norm(s.exponents, s.l))
        for s in raw]
    written_shells = [molden._RawShell(s.atom, s.l, s.exponents, c)
                      for s, c in zip(raw, coefficients)]
    target = molden._build_basis(written_shells, _CENTRES, pure, convention)
    cartesian = B.to_cartesian(reference, orbitals)
    written = np.zeros((len(orbitals), target.nbf))
    column = 0
    for shell, start in zip(target.shells, target.offsets()):
        n = B.ncart(shell.l)
        written[:, start:start + shell.nfunc] = (
            cartesian[:, column:column + n] @ np.linalg.pinv(shell.transform))
        column += n
    return coefficients, written


def _molden_text(pure: set[int], coefficients, orbitals, title="written by a test") -> str:
    lines = ["[Molden Format]", "[Title]", f" {title}", "[Atoms] AU"]
    for number, (element, z) in zip(_CENTRES, (("N", 7), ("Cu", 29))):
        x, y, zz = _CENTRES[number]
        lines.append(f"{element} {number} {z} {x:.10f} {y:.10f} {zz:.10f}")
    lines.append("[GTO]")
    atom = None
    for (number, letter, exponents, _), contraction in zip(_SHELLS, coefficients):
        if number != atom:
            if atom is not None:
                lines.append("")
            lines.append(f"  {number} 0")
            atom = number
        lines.append(f" {letter} {len(exponents)} 1.00")
        lines += [f"  {e:.10E} {c:.10E}" for e, c in zip(exponents, contraction)]
    lines.append("")
    if 2 in pure:
        lines.append("[5D7F]")
    if 4 in pure:
        lines.append("[9G]")
    lines.append("[MO]")
    for k, row in enumerate(orbitals):
        lines += [" Sym= A", f" Ene= {-1.0 + 0.1 * k:.6f}", " Spin= Alpha",
                  f" Occup= {2.0 if k < 3 else 0.0}"]
        lines += [f"{i + 1:5d} {c:.12f}" for i, c in enumerate(row)]
    return "\n".join(lines) + "\n"


_DIALECTS = [
    ("standard", {2, 3, 4}, molden.Convention(), "standard"),
    ("ORCA", {2, 3, 4}, molden.Convention(normalized_primitives=False, orca_signs=True),
     "primitive norms in the coefficients, ORCA's signs for pure f, g, h"),
    ("old Psi4", {2, 3, 4}, molden.Convention(normalized_primitives=False),
     "primitive norms in the coefficients"),
    ("Psi4 1.3", set(), molden.Convention(cartesian="shell"), "Cartesians normalized as x^l"),
    ("CFOUR", set(), molden.Convention(cartesian="raw"), "Cartesians normalized as xy"),
    ("Turbomole", set(), molden.Convention(cartesian="turbomole"),
     "Cartesians scaled by (2l-1)!!"),
]


@pytest.mark.parametrize("program, pure, convention, name", _DIALECTS,
                         ids=[d[0] for d in _DIALECTS])
def test_each_programs_molden_dialect_is_recognized(tmp_path, program, pure, convention,
                                                    name):
    raw, reference = _reference_basis(pure)
    truth = _orthonormal_orbitals(reference, 12)
    coefficients, written = _as_written(raw, reference, truth, pure, convention)
    path = tmp_path / "wavefunction.molden"
    path.write_text(_molden_text(pure, coefficients, written))

    orbitals = microscope.load(path).orbitals
    assert orbitals.convention == name
    assert orbitals.orthonormality_error(occupied_only=False) < 1e-6   # float32 storage
    # the orbital drawn is the orbital written
    points = np.random.default_rng(1).normal(scale=1.5, size=(300, 3)) + _CENTRES[1]
    for k in (0, 5, 11):
        assert np.allclose(_brute_force(orbitals.basis, orbitals.alpha.coefficients[k], points),
                           _brute_force(reference, truth[k], points), atol=1e-5)


def test_orbitals_no_reading_explains_are_flagged(tmp_path):
    raw, reference = _reference_basis({2, 3, 4})
    junk = np.random.default_rng(5).normal(size=(6, reference.nbf))
    path = tmp_path / "junk.molden"
    path.write_text(_molden_text({2, 3, 4}, [s.coefficients for s in raw], junk))
    assert microscope.load(path).orbitals.convention == "unrecognized"


def test_molden_atoms_in_parenthesized_au_are_bohr(tmp_path):
    raw, reference = _reference_basis(set())
    text = _molden_text(set(), [s.coefficients for s in raw],
                        _orthonormal_orbitals(reference, 4))
    path = tmp_path / "au.molden"
    path.write_text(text.replace("[Atoms] AU", "[Atoms] (AU)"))
    coords = microscope.load(path).molecule.coords
    assert np.allclose(coords[1], _CENTRES[2] * B.BOHR_TO_ANGSTROM)


def test_a_ghost_atom_keeps_its_basis_functions_but_is_not_drawn(tmp_path):
    raw, reference = _reference_basis(set())
    text = _molden_text(set(), [s.coefficients for s in raw],
                        _orthonormal_orbitals(reference, 4))
    path = tmp_path / "ghost.molden"
    path.write_text(text.replace("Cu 2 29", "Cu 2 0"))
    result = microscope.load(path)
    assert result.molecule.symbols == ["N"]
    assert len({s.atom for s in result.orbitals.basis.shells}) == 2
    assert result.orbitals.convention == "standard"


def test_orca_molden_input_files_open(tmp_path):
    raw, reference = _reference_basis(set())
    path = tmp_path / "job.molden.input"
    path.write_text(_molden_text(set(), [s.coefficients for s in raw],
                                 _orthonormal_orbitals(reference, 4)))
    assert microscope.load(path).orbitals is not None


# ------------------------------------------------------- what real files do

def _small_molden(tmp_path, name="small.molden", edit=lambda text: text, pure=(2, 3, 4)):
    pure = set(pure)
    raw, reference = _reference_basis(pure)
    text = _molden_text(pure, [s.coefficients for s in raw],
                        _orthonormal_orbitals(reference, 6))
    path = tmp_path / name
    path.write_text(edit(text))
    return path


def test_an_ecp_atom_is_named_by_its_symbol_not_its_core_charge(tmp_path):
    """ORCA writes the charge the calculation used: 17 for a Rh with a 28-electron core."""
    path = _small_molden(tmp_path, edit=lambda t: t.replace("Cu 2 29", "Rh 2 17"))
    assert microscope.load(path).molecule.symbols == ["N", "Rh"]


def test_a_broken_wavefunction_costs_the_orbitals_not_the_structure(tmp_path):
    path = _small_molden(tmp_path, edit=lambda t: t.replace("  2 0\n", "  7 0\n"))
    result = microscope.load(path)
    assert result.orbitals is None and result.molecule.natoms == 2
    assert "atom 7" in result.warnings[0]


def test_pure_functions_are_counted_when_the_flags_are_missing(tmp_path):
    """Five d functions written without [5D] mean five, not six."""
    path = _small_molden(tmp_path, edit=lambda t: t.replace("[5D7F]\n", "").replace("[9G]\n", ""))
    orbitals = microscope.load(path).orbitals
    assert orbitals.convention == "standard"
    assert all(s.pure for s in orbitals.basis.shells if s.l >= 2)


def test_a_coefficient_numbered_zero_is_refused(tmp_path):
    path = _small_molden(tmp_path, edit=lambda t: t.replace("\n    1 ", "\n    0 ", 1))
    result = microscope.load(path)
    assert result.orbitals is None and "counts from 1" in result.warnings[0]


@pytest.mark.parametrize("ending", ["\r\n", "\r\r\n"], ids=["crlf", "crlf-twice"])
def test_windows_line_endings_read_the_same(tmp_path, ending):
    """CRLF, and the CR CR LF a file gets from crossing between Windows and
    Unix twice. Written as bytes: write_text would translate the endings
    again on Windows."""
    plain = microscope.load(_small_molden(tmp_path)).orbitals
    text = (tmp_path / "small.molden").read_bytes().decode().replace("\r\n", "\n")
    path = tmp_path / "windows.molden"
    path.write_bytes(text.replace("\n", ending).encode())
    windows = microscope.load(path).orbitals
    assert windows.convention == plain.convention == "standard"
    assert np.array_equal(windows.alpha.coefficients, plain.alpha.coefficients)


def test_homo_needs_occupations(tmp_path):
    path = _small_molden(tmp_path, edit=lambda t: t.replace("Occup= 2.0", "Occup= 0.0"))
    orbitals = microscope.load(path).orbitals
    with pytest.raises(ValueError, match="no occupations"):
        orbitals.find("homo")
    assert orbitals.find("2") == ("alpha", 1)


def test_a_cut_short_fchk_keeps_its_geometry(tmp_path):
    text = (DATA / "dvb_ir.fchk").read_text()
    start = text.index("Alpha MO coefficients")
    cut = text[:start] + "\n".join(text[start:].splitlines()[:40]) + "\n"
    (tmp_path / "cut.fchk").write_text(cut)
    result = microscope.load(tmp_path / "cut.fchk")
    assert result.orbitals is None and result.molecule.natoms == 20
    assert "cut short" in result.warnings[0]


def test_an_fchk_meaning_something_else_is_flagged(tmp_path):
    """Another program's fchk whose functions are scaled its own way must not
    be drawn as if it were Gaussian's."""
    orbitals = microscope.load(DATA / "dvb_ir.fchk").orbitals
    assert orbitals.convention == "standard"
    orbitals.alpha.coefficients[:, 1::4] *= 1.7          # as if some functions meant more
    assert orbitals.reading_error() > 1e-3


def test_a_diffuse_orbital_gets_room_instead_of_being_sliced(tmp_path):
    rows = B.cartesian_rows(0, ["s"])
    basis = B.BasisSet([B.make_shell([0, 0, 0], 0, [0.02], [1.0], rows)])
    from microscope.core.orbitals import EDGE, Orbitals, OrbitalSet, _at_walls
    rydberg = Orbitals(basis, OrbitalSet(np.zeros(1), np.zeros(1, np.float32),
                                         np.ones((1, 1), np.float32)))
    volume = rydberg.volume("1")
    assert _at_walls(volume.values) < EDGE
    assert volume.origin[0] < -4.0 - 1.0                  # it grew past the 4 A start


def test_the_check_samples_big_files_but_keeps_the_homo():
    orbitals = _orbital_set_with(nmo=2000, homo=731)
    rows = orbitals.sample(limit=600)
    assert len(rows) <= 600
    assert any(np.array_equal(r, orbitals.coefficients[731]) for r in rows)


def _orbital_set_with(nmo: int, homo: int):
    from microscope.core.orbitals import OrbitalSet
    occupations = np.zeros(nmo, np.float32)
    occupations[:homo + 1] = 2
    coefficients = np.arange(nmo, dtype=np.float32)[:, None] * np.ones((1, 3), np.float32)
    return OrbitalSet(np.zeros(nmo), occupations, coefficients)
