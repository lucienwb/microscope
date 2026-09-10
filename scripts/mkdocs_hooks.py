"""Build-time hooks for the documentation site.

Two reference pages are generated rather than written, because they describe
things the code already defines and a hand-kept copy would drift the first
time a flag or a style key changed:

    <!-- command-line-reference -->   every `scope` flag, from the argparse parser
    <!-- style-reference -->          every style-file key, from the Style dataclass

MkDocs loads this through `hooks:` in mkdocs.yml, so there is no plugin to
install.
"""

from __future__ import annotations

import sys
from dataclasses import fields
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from microscope.cli import build_parser  # noqa: E402
from microscope.render.styles import COLOR_FIELDS, Style  # noqa: E402


def _cell(text: str) -> str:
    return " ".join(str(text).split()).replace("|", "\\|")


def command_line_reference() -> str:
    """One table per argument group, in the order `scope --help` prints them."""
    parser = build_parser()
    out = []
    for group in parser._action_groups:
        actions = [a for a in group._group_actions
                   if a.option_strings or a.dest == "file"]
        if not actions:
            continue
        title = group.title[:1].upper() + group.title[1:]
        out += [f"### {title}", "", "| Flag | Meaning |", "|---|---|"]
        for action in actions:
            if action.option_strings:
                names = " · ".join(f"`{o}`" for o in action.option_strings)
                if action.metavar or (action.nargs != 0 and action.const is None
                                      and not isinstance(action.default, bool)):
                    metavar = action.metavar or action.dest.upper()
                    names += f" `{metavar}`"
            else:
                names = "`FILE`"
            text = _cell(action.help or "")
            if text and action.default not in (None, False, "==SUPPRESS==") \
                    and "default" not in text:
                text += f" *(default: `{action.default}`)*"
            out.append(f"| {names} | {text} |")
        out.append("")
    return "\n".join(out)


def style_reference() -> str:
    """Every key a style file may hold, with the CYLview default for each."""
    default = Style()
    rows = ["| Key | Default | What it is |", "|---|---|---|"]
    notes = {
        "name": "shown in the viewer's style list",
        "bond_radius": "bond cylinder radius, Å",
        "atom_scale": "atom radius as a fraction of the covalent radius",
        "min_atom_radius": "smallest atom radius, Å — keeps hydrogens visible",
        "palette": "colour per element, keyed by symbol",
        "background": "the page behind the molecule",
        "show_hbonds": "draw hydrogen bonds as dashed lines",
        "hbond_radius": "radius of the dashes, Å",
        "hbond_color": "colour of the dashes",
        "hbond_dash": "dash length, Å",
        "hbond_gap": "gap between dashes, Å",
        "surface_positive": "isosurface colour for the + lobe",
        "surface_negative": "isosurface colour for the − lobe",
        "surface_opacity": "isosurface opacity, 0–1",
        "bond_color": "one colour for every bond; `null` splits each by its atoms",
        "quadrant_color": "seam lines on heavy atoms (the Houk style); `null` for none",
        "quadrant_width": "seam line width, as a fraction of the atom radius",
    }
    for f in fields(Style):
        value = getattr(default, f.name)
        if f.name == "palette":
            shown = f"{len(value)} elements"
        elif f.name in COLOR_FIELDS:
            shown = "`null`" if value is None else "`#" + "".join(
                f"{round(c * 255):02x}" for c in value) + "`"
        elif isinstance(value, bool):
            shown = f"`{str(value).lower()}`"
        elif isinstance(value, str):
            shown = f'`"{value}"`'
        else:
            shown = f"`{value}`"
        rows.append(f"| `{f.name}` | {shown} | {notes.get(f.name, '')} |")
    return "\n".join(rows)


GENERATED = {
    "<!-- command-line-reference -->": command_line_reference,
    "<!-- style-reference -->": style_reference,
}


def on_page_markdown(markdown, page, config, files):
    for marker, build in GENERATED.items():
        if marker in markdown:
            markdown = markdown.replace(marker, build())
    return markdown
