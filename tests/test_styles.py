"""Rendering styles, and the files that carry them between the two front ends.

The point of a style file is a group standard: saved once in the viewer, used
by `scope -s --style ours.json` for every figure in a paper. So what matters
is that it round-trips exactly, reads when written by hand, and says clearly
what is wrong when it does not.
"""

import json

import pytest

from microscope import cli
from microscope.render import styles


def test_a_style_round_trips_through_a_file(tmp_path):
    original = styles.houk_style()
    original.palette[6] = (0.10, 0.35, 0.55)
    original.atom_scale = 0.42
    path = tmp_path / "ours.json"
    styles.save_style(path, original)
    assert styles.style_to_dict(styles.load_style(path)) == \
        styles.style_to_dict(original)


def test_a_style_file_is_written_to_be_edited(tmp_path):
    """Elements by symbol and colours as #rrggbb, not indices and float triples."""
    path = tmp_path / "ours.json"
    styles.save_style(path, styles.houk_style())
    data = json.loads(path.read_text())
    assert data["palette"]["C"] == "#e6e6e6"
    assert data["background"] == "#ffffff"
    assert set("0123456789abcdef#").issuperset(data["bond_color"])


def test_a_hand_written_style_needs_only_what_it_changes(tmp_path):
    path = tmp_path / "ours.json"
    path.write_text(json.dumps({"name": "ours", "atom_scale": 0.42,
                                "palette": {"C": "#4d4d4d"}}))
    style = styles.load_style(path)
    assert style.atom_scale == 0.42
    assert style.atom_color(6) == pytest.approx((0.302, 0.302, 0.302), abs=0.002)
    assert style.bond_radius == styles.Style().bond_radius     # untouched


@pytest.mark.parametrize("body, complaint", [
    ({"colour": "#ffffff"}, "not part of a style"),
    ({"palette": {"Xx": "#ffffff"}}, "not an element"),
    ({"background": "red"}, "not #rrggbb"),
    ({"palette": {"C": [0.1, 0.2]}}, "three numbers"),
])
def test_a_broken_style_file_says_what_is_wrong(tmp_path, body, complaint):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(body))
    with pytest.raises(styles.StyleError, match=complaint):
        styles.load_style(path)


def test_invalid_json_is_not_a_traceback(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{oops")
    with pytest.raises(styles.StyleError, match="not valid JSON"):
        styles.load_style(path)


def test_make_style_takes_a_preset_or_a_file(tmp_path):
    assert styles.make_style("houk").name == "houk"
    path = tmp_path / "ours.json"
    styles.save_style(path, styles.Style(name="ours"))
    assert styles.make_style(str(path)).name == "ours"


def test_an_unknown_style_names_the_presets():
    with pytest.raises(styles.StyleError, match="cylview, houk"):
        styles.make_style("whatever")


def test_the_command_line_takes_a_style_file(tmp_path):
    path = tmp_path / "ours.json"
    styles.save_style(path, styles.Style(name="ours", atom_scale=0.5))
    args = cli.build_parser().parse_args(
        ["-s", "mol.log", "--style", str(path)])
    assert cli.build_style(args).atom_scale == 0.5


def test_a_bad_style_on_the_command_line_is_a_clean_error():
    args = cli.build_parser().parse_args(["-s", "mol.log", "--style", "nope"])
    with pytest.raises(cli.CliError, match="no style preset or file"):
        cli.build_style(args)


# ------------------------------------------------------------ depth cueing

def test_depth_cueing_is_off_until_asked_for(tmp_path):
    """Some people like the far side faded and some do not: it is a style
    setting, off in both presets, kept in a style file like the rest."""
    assert styles.Style().depth_cue == 0.0 and styles.houk_style().depth_cue == 0.0
    style = styles.Style(depth_cue=0.6)
    styles.save_style(tmp_path / "cued.json", style)
    assert styles.load_style(tmp_path / "cued.json").depth_cue == 0.6


def test_the_command_line_sets_depth_cueing():
    args = cli.build_parser().parse_args(["-s", "x.log", "--depth-cue", "0.5"])
    assert cli.build_style(args).depth_cue == 0.5
    assert cli.build_style(cli.build_parser().parse_args(["-s", "x.log"])).depth_cue == 0.0
    with pytest.raises(cli.CliError):
        cli.build_style(cli.build_parser().parse_args(["-s", "x.log", "--depth-cue", "2"]))
    parser = cli.build_parser()
    with pytest.raises(SystemExit):                    # a Lewis drawing is flat
        cli.check_mode_flags(parser.parse_args(["-s", "x.log", "--lewis",
                                                "--depth-cue", "0.5"]), parser)


def test_depth_cueing_runs_from_the_front_of_the_molecule_to_its_back():
    import numpy as np

    from microscope.render.glrenderer import MoleculeRenderer

    renderer = MoleculeRenderer()                      # no GL needed for the range
    renderer._atoms = np.array([[0, 0, 2.0, 0.5], [0, 0, -3.0, 0.5]], dtype=np.float32)
    view = np.eye(4)
    view[2, 3] = -10.0                                 # the camera stands 10 A off
    near, far = renderer._depth_range(view)
    assert (near, far) == pytest.approx((-7.5, -13.5))
    assert renderer._depth_range(np.eye(4)) == pytest.approx((2.5, -3.5))


def test_every_style_key_is_described_in_the_style_reference():
    """The style-file reference is generated from the Style dataclass, but its
    descriptions are written by hand: a new key must get one."""
    import importlib.util
    from dataclasses import fields
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "mkdocs_hooks", Path(__file__).parent.parent / "scripts" / "mkdocs_hooks.py")
    hooks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hooks)
    table = hooks.style_reference()
    for f in fields(styles.Style):
        row = next(line for line in table.splitlines() if line.startswith(f"| `{f.name}` |"))
        assert not row.rstrip().endswith("|  |"), f"{f.name} has no description"
