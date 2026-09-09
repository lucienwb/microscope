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
