"""The `scope` command line: argument parsing, view setup and error reporting.

Headless — the rendering itself needs a GL context, so these tests cover
everything up to the point where pixels are drawn.
"""

import numpy as np
import pytest

from microscope import cli
from microscope.core.molecule import Molecule
from microscope.core.results import ExcitedState, NMRShielding, ParseResult, Vibration
from microscope.render.scene import REP_BALL, REP_LINE, REP_STICK


def _molecule():
    return Molecule(["O", "H", "H"],
                    np.array([[0.0, 0.0, 0.0],
                              [0.96, 0.0, 0.0],
                              [-0.24, 0.93, 0.0]]))


# ------------------------------------------------------------------ small parsers

def test_parse_size():
    assert cli.parse_size("1200x900") == (1200, 900)
    assert cli.parse_size("800X600") == (800, 600)
    for bad in ("1200", "12x", "axb", "8x8"):
        with pytest.raises(cli.CliError):
            cli.parse_size(bad)


def test_parse_indices_is_one_based():
    assert cli.parse_indices("1") == [0]
    assert cli.parse_indices("1-3,7") == [0, 1, 2, 6]
    assert cli.parse_indices(" 2 , 4 ") == [1, 3]


def test_parse_indices_rejects_bad_input():
    for bad in ("0", "-1", "3-1", "x", ""):
        with pytest.raises(cli.CliError):
            cli.parse_indices(bad)
    with pytest.raises(cli.CliError):
        cli.parse_indices("5", natoms=3)


def test_parse_color_accepts_names_and_hex():
    assert cli.parse_color("white") == pytest.approx((1.0, 1.0, 1.0))
    assert cli.parse_color("#ff0000") == pytest.approx((1.0, 0.0, 0.0))
    with pytest.raises(cli.CliError):
        cli.parse_color("not-a-color")


def test_parse_reps_whole_molecule_and_ranges():
    assert list(cli.parse_reps("stick", 4)) == [REP_STICK] * 4
    reps = cli.parse_reps("1-2:ball;3,4:line", 4)
    assert list(reps) == [REP_BALL, REP_BALL, REP_LINE, REP_LINE]
    with pytest.raises(cli.CliError):
        cli.parse_reps("1-2:blob", 4)


def test_parse_rotation():
    assert cli.parse_rotation("30,-15") == (30.0, -15.0)
    for bad in ("30", "a,b"):
        with pytest.raises(cli.CliError):
            cli.parse_rotation(bad)


def test_output_paths():
    assert cli.default_output("/tmp/mycalc.log") == "mycalc.png"
    assert cli.check_output("fig.tif") == "fig.tif"
    for bad in ("fig.jpg", "fig"):
        with pytest.raises(cli.CliError):
            cli.check_output(bad)


# ------------------------------------------------------------------ frames & view

def test_pick_frame():
    frames = [_molecule() for _ in range(3)]
    assert cli.pick_frame(frames, None) is frames[-1]
    assert cli.pick_frame(frames, 1) is frames[0]
    assert cli.pick_frame(frames, -1) is frames[-1]
    for bad in (0, 4, -4):
        with pytest.raises(cli.CliError):
            cli.pick_frame(frames, bad)
    with pytest.raises(cli.CliError):
        cli.pick_frame([], None)


def test_named_views_are_proper_rotations():
    for name, rot in cli.NAMED_VIEWS.items():
        assert np.allclose(rot @ rot.T, np.eye(3)), name
        assert np.isclose(np.linalg.det(rot), 1.0), name


def test_named_view_looks_down_the_expected_axis():
    mol = _molecule()
    # "top" looks down the world y axis: it must map world +y to view +z
    cam = cli.build_camera(mol, "top", None, None, 1.0)
    assert np.allclose(cam.rotation @ np.array([0.0, 1.0, 0.0]), [0, 0, 1])
    front = cli.build_camera(mol, "front", None, None, 1.0)
    assert np.allclose(front.rotation, np.eye(3))


def test_view_angles_turn_the_camera():
    mol = _molecule()
    cam = cli.build_camera(mol, "30,-15", None, None, 1.0)
    assert not np.allclose(cam.rotation, np.eye(3))
    assert np.allclose(cam.rotation @ cam.rotation.T, np.eye(3))


def test_rotate_applies_on_top_of_view():
    mol = _molecule()
    base = cli.build_camera(mol, "top", None, None, 1.0)
    spun = cli.build_camera(mol, "top", None, "45,0", 1.0)
    assert not np.allclose(base.rotation, spun.rotation)


def test_align_two_atoms_looks_down_the_bond():
    mol = _molecule()
    cam = cli.build_camera(mol, None, "1,2", None, 1.0)
    bond = mol.coords[1] - mol.coords[0]
    assert np.allclose(cam.rotation @ bond / np.linalg.norm(bond), [0, 0, 1],
                       atol=1e-9)


def test_align_three_atoms_puts_the_plane_in_the_screen():
    mol = _molecule()
    cam = cli.build_camera(mol, None, "1,2,3", None, 1.0)
    for i in (0, 1, 2):                      # every atom lands at the same depth
        assert np.isclose((cam.rotation @ mol.coords[i])[2],
                          (cam.rotation @ mol.coords[0])[2], atol=1e-9)


def test_align_rejects_wrong_atom_counts_and_collinear_atoms():
    mol = _molecule()
    with pytest.raises(cli.CliError):
        cli.build_camera(mol, None, "1,2,3,1", None, 1.0)
    line = Molecule(["H", "H", "H"], np.array([[0.0, 0, 0], [1.0, 0, 0],
                                               [2.0, 0, 0]]))
    with pytest.raises(cli.CliError):
        cli.build_camera(line, None, "1,2,3", None, 1.0)


def test_zoom_shrinks_the_view_height():
    mol = _molecule()
    plain = cli.build_camera(mol, None, None, None, 1.0)
    close = cli.build_camera(mol, None, None, None, 2.0)
    assert np.isclose(close.half_height, plain.half_height / 2.0)
    with pytest.raises(cli.CliError):
        cli.build_camera(mol, None, None, None, -1.0)


# ------------------------------------------------------------------ style & args

def test_build_style_applies_overrides():
    args = cli.build_parser().parse_args(
        ["-s", "x.log", "--style", "houk", "--background", "black",
         "--bond-color", "#101010", "--no-hbonds",
         "--iso-colors", "blue,red", "--iso-opacity", "0.4"])
    style = cli.build_style(args)
    assert style.name == "houk"
    assert style.background == pytest.approx((0.0, 0.0, 0.0))
    assert style.show_hbonds is False
    assert style.surface_positive == pytest.approx((0.0, 0.0, 1.0))
    assert style.surface_opacity == pytest.approx(0.4)


def test_iso_colors_needs_two_colors():
    args = cli.build_parser().parse_args(["-s", "x.log", "--iso-colors", "blue"])
    with pytest.raises(cli.CliError):
        cli.build_style(args)


def test_gui_mode_refuses_render_only_flags():
    parser = cli.build_parser()
    args = parser.parse_args(["mol.log", "--size", "100x100"])
    with pytest.raises(SystemExit):
        cli.check_gui_flags(args, parser)


def test_gui_mode_allows_the_shared_flags():
    parser = cli.build_parser()
    args = parser.parse_args(["mol.log", "--style", "houk", "--labels",
                              "number", "--axes"])
    cli.check_gui_flags(args, parser)          # must not raise


def test_silent_without_a_file_is_an_error(capsys):
    assert cli.main(["-s"]) == 2
    assert "give a file" in capsys.readouterr().err


def test_missing_file_is_reported_cleanly(capsys):
    assert cli.main(["-s", "definitely-not-here.log"]) == 2
    assert "no such file" in capsys.readouterr().err


# ------------------------------------------------------------------ spectra

def _spectral_result():
    return ParseResult(
        frames=[_molecule()],
        vibrations=[Vibration(1600.0, ir_intensity=40.0),
                    Vibration(3100.0, ir_intensity=10.0)],
        excited_states=[ExcitedState(1, "S1", 4.0, 309.96, 0.5)],
        nmr_shieldings=[NMRShielding(0, "C", 120.0), NMRShielding(1, "H", 25.0)])


def _args(*argv):
    return cli.build_parser().parse_args(list(argv))


def test_spectrum_kind_detection():
    assert cli.spectrum_kind(_args("-s", "x.log")) is None
    assert cli.spectrum_kind(_args("-s", "x.log", "--ir")) == "ir"
    assert cli.spectrum_kind(_args("-s", "x.log", "--uvvis")) == "uv"
    assert cli.spectrum_kind(_args("-s", "x.log", "--nmr")) == "nmr"


def test_only_one_spectrum_at_a_time():
    with pytest.raises(SystemExit):
        _args("-s", "x.log", "--ir", "--nmr")


def test_parse_range():
    assert cli.parse_range("400:1800") == (400.0, 1800.0)
    assert cli.parse_range("2..7") == (2.0, 7.0)
    for bad in ("400", "a:b", "5:5"):
        with pytest.raises(cli.CliError):
            cli.parse_range(bad)


def test_spectra_get_a_wider_default_frame():
    args = _args("-s", "x.log", "--ir")
    assert cli.resolve_size(args, "ir") == (1600, 900)
    assert cli.resolve_size(args) == (1200, 900)
    assert cli.resolve_size(_args("-s", "x.log", "--size", "640x480"), "ir") == (640, 480)


def test_spectrum_output_names_and_formats():
    assert cli.default_output("/tmp/mycalc.log", "ir") == "mycalc_ir.png"
    for good in ("s.pdf", "s.svg", "s.eps", "s.png"):
        assert cli.check_output(good, "ir") == good
    with pytest.raises(cli.CliError):
        cli.check_output("s.jpg", "ir")
    with pytest.raises(cli.CliError):
        cli.check_output("s.pdf")          # structure renders stay PNG/TIFF


def test_build_ir_spectrum_applies_scale_and_fwhm():
    result = _spectral_result()
    spec = cli.build_spectrum("ir", result, _args("-s", "x.log", "--ir"))
    assert spec.invert_x and spec.stick_x[0] == pytest.approx(1600.0)
    scaled = cli.build_spectrum(
        "ir", result, _args("-s", "x.log", "--ir", "--freq-scale", "0.96"))
    assert scaled.stick_x[0] == pytest.approx(1536.0)


def test_default_fwhm_per_spectrum_type():
    result = _spectral_result()
    narrow = cli.build_spectrum("ir", result,
                                _args("-s", "x.log", "--ir", "--fwhm", "2"))
    wide = cli.build_spectrum("ir", result,
                              _args("-s", "x.log", "--ir", "--fwhm", "40"))
    assert wide.y.max() > narrow.y.max()      # broader lines overlap more
    with pytest.raises(cli.CliError):
        cli.build_spectrum("ir", result,
                           _args("-s", "x.log", "--ir", "--fwhm", "0"))


def test_build_uv_spectrum_units():
    result = _spectral_result()
    nm = cli.build_spectrum("uv", result, _args("-s", "x.log", "--uv"))
    ev = cli.build_spectrum("uv", result,
                            _args("-s", "x.log", "--uv", "--unit", "eV"))
    assert "nm" in nm.xlabel and "eV" in ev.xlabel


def test_build_nmr_spectrum_reference_and_nucleus():
    result = _spectral_result()
    raw = cli.build_spectrum("nmr", result, _args("-s", "x.log", "--nmr"))
    assert not raw.invert_x                    # sigma scale by default
    assert raw.meta["element"] == "H"          # H preferred when not asked
    shifted = cli.build_spectrum(
        "nmr", result,
        _args("-s", "x.log", "--nmr", "--nucleus", "C", "--reference", "186.4"))
    assert shifted.invert_x
    assert shifted.stick_x[0] == pytest.approx(66.4)


def test_pick_nucleus():
    sh = [NMRShielding(0, "C", 1.0), NMRShielding(1, "H", 2.0)]
    assert cli.pick_nucleus(sh, None) == "H"
    assert cli.pick_nucleus([NMRShielding(0, "C", 1.0)], None) == "C"
    assert cli.pick_nucleus(sh, "C") == "C"
    with pytest.raises(cli.CliError):
        cli.pick_nucleus(sh, "F")
    with pytest.raises(cli.CliError):
        cli.pick_nucleus([], None)


def test_missing_spectral_data_is_explained():
    empty = ParseResult(frames=[_molecule()])
    for kind, hint in (("ir", "freq"), ("uv", "TD-DFT"), ("nmr", "NMR")):
        with pytest.raises(cli.CliError) as exc:
            cli.build_spectrum(kind, empty, _args("-s", "x.log", f"--{kind}"))
        assert hint in str(exc.value)


def test_structure_and_spectrum_flags_do_not_mix():
    parser = cli.build_parser()
    with pytest.raises(SystemExit):        # --style means nothing for a plot
        cli.check_mode_flags(parser.parse_args(
            ["-s", "x.log", "--ir", "--style", "houk"]), parser)
    with pytest.raises(SystemExit):        # --fwhm means nothing for a render
        cli.check_mode_flags(parser.parse_args(
            ["-s", "x.log", "--fwhm", "10"]), parser)
    cli.check_mode_flags(parser.parse_args(
        ["-s", "x.log", "--ir", "--fwhm", "10", "--size", "800x600"]), parser)
    cli.check_mode_flags(parser.parse_args(
        ["-s", "x.log", "--style", "houk"]), parser)
