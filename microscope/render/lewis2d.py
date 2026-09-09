"""Laying a Lewis structure out flat, the way ChemDraw draws one.

There is no 2-D coordinate generator here and there deliberately never will
be: the drawing *is* the molecule seen through the same orthographic camera as
the 3-D view, so turning the structure picks the viewing angle the drawing is
made from. Choosing that angle is the whole point of the mode — a chemist
rotates until the ring is face-on and the substituents are clear, then saves.

Everything is plain numpy so the geometry can be tested without a display.
Pixel positions come out of ``camera.project``; sizes (font, line width, the
spacing between the two lines of a double bond) are tied to the *unprojected*
median bond length instead, so they hold still while the structure turns
rather than shrinking as a bond swings edge-on.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core import elements
from ..core.lewis import LewisStructure, is_main_group
from .camera import OrthoCamera
from .styles import CYLVIEW_COLORS

# Ink on a white page, ChemDraw-style; --color-atoms tints the heteroatoms.
BLACK = (0.0, 0.0, 0.0)
_LABEL_COLORS = {**CYLVIEW_COLORS, 1: (0.15, 0.15, 0.15), 6: BLACK}

FONT_FRACTION = 0.42        # label height as a fraction of a bond
LINE_FRACTION = 0.052       # bond line width
DOUBLE_GAP = 0.17           # spacing between the lines of a multiple bond
INNER_INSET = 0.14          # how far the inner line of a ring bond is cut back
LABEL_GAP = 0.20            # blank space between a label and its bonds (of font)
DOT_RADIUS = 0.075          # lone-pair dot radius, as a fraction of a bond
BRACKET_PAD = 0.30          # gap between the drawing and the brackets
BRACKET_ARM = 0.16          # length of a bracket's top and bottom arms
OVERLAP_FRACTION = 0.35     # atoms closer than this (of a bond) are stacked up
                            # 0.30 let a ring seen edge-on collapse to a stick
                            # unreported; above ~0.45 the crowded structures
                            # have no clean view left to turn to


@dataclass
class LewisOptions:
    """What to draw. The defaults are the skeletal ChemDraw look."""

    show_hydrogens: bool = False    # draw every H instead of folding it into a label
    carbon_labels: bool = False     # write C on carbons instead of leaving a vertex
    lone_pairs: bool = False        # non-bonding electrons as dots
    color_atoms: bool = False       # tint labels by element instead of all black


@dataclass
class Atom2D:
    index: int
    pos: np.ndarray                        # (2,) pixels
    runs: tuple[tuple[str, bool], ...] = ()  # (text, is_subscript)
    charge: int = 0
    radicals: int = 0
    lone_pairs: int = 0
    color: tuple = BLACK
    half_width: float = 0.0                # the text box itself
    half_height: float = 0.0
    gap: float = 0.0                       # blank space kept around the text
    charge_text: str = ""
    charge_pos: np.ndarray | None = None
    dots: list = field(default_factory=list)   # (x, y) lone-pair/radical dots

    @property
    def text(self) -> str:
        return "".join(t for t, _ in self.runs)

    @property
    def labelled(self) -> bool:
        return bool(self.runs)

    def edge(self, direction: np.ndarray) -> float:
        """Distance from the centre to the edge of the text along *direction*."""
        if not self.labelled:
            return 0.0
        dx, dy = float(abs(direction[0])), float(abs(direction[1]))
        far = 1e9
        return min(self.half_width / dx if dx > 1e-9 else far,
                   self.half_height / dy if dy > 1e-9 else far)

    def clearance(self, direction: np.ndarray) -> float:
        """Where a bond line has to stop: the text plus its breathing space."""
        return self.edge(direction) + (self.gap if self.labelled else 0.0)


@dataclass
class Bond2D:
    i: int
    j: int
    order: int
    lines: list = field(default_factory=list)   # [(start, end)] pixel pairs


@dataclass
class Brackets:
    """Square brackets round the whole drawing, with the species charge.

    What a chemist reaches for when a charge or an unpaired electron belongs
    to the molecule rather than to any one atom of it — a delocalized radical
    cation being the usual case.
    """

    left: float
    right: float
    top: float
    bottom: float
    arm: float
    runs: tuple                    # the label, as (text, is_superscript) runs
    label_pos: np.ndarray | None = None


@dataclass
class LewisLayout:
    atoms: list[Atom2D]
    bonds: list[Bond2D]
    scale: float                  # pixels per bond length
    font_px: float
    line_width: float
    overlaps: int = 0
    warnings: list[str] = field(default_factory=list)
    brackets: Brackets | None = None


def default_measure(runs, font_px: float) -> tuple[float, float]:
    """Rough text extent, used when no real font metrics are available."""
    width = sum(len(text) * font_px * (0.42 if sub else 0.62) for text, sub in runs)
    return width, font_px * 0.72


def hidden_hydrogens(structure: LewisStructure, options: LewisOptions) -> np.ndarray:
    """Hydrogens folded into a neighbour's label instead of being drawn.

    Only a hydrogen with exactly one heavy neighbour disappears, so H2 and a
    stray proton still show up rather than vanishing from the drawing.
    """
    mol = structure.molecule
    hidden = np.zeros(mol.natoms, dtype=bool)
    if options.show_hydrogens:
        return hidden
    for a in range(mol.natoms):
        if mol.symbols[a] != "H":
            continue
        nbrs = structure.neighbors[a]
        if len(nbrs) == 1 and mol.symbols[nbrs[0]] != "H":
            hidden[a] = True
    return hidden


def _label_runs(structure, atom, hidden, flip) -> tuple:
    """Text of an atom label, split into normal and subscript runs.

    Only the hydrogens actually left out of the drawing are folded in — a
    bridging hydrogen stays a vertex of its own and must not be counted here
    as well.
    """
    symbol = structure.molecule.symbols[atom]
    runs = [(symbol, False)]
    count = sum(1 for j in structure.neighbors[atom]
                if structure.molecule.symbols[j] == "H" and hidden[j])
    if count:
        tail = [("H", False)] + ([(str(count), True)] if count > 1 else [])
        runs = tail + runs if flip else runs + tail
    return tuple(runs)


def _charge_text(charge: int) -> str:
    if not charge:
        return ""
    sign = "+" if charge > 0 else "−"
    return sign if abs(charge) == 1 else f"{abs(charge)}{sign}"


def layout(structure: LewisStructure, camera: OrthoCamera,
           width: float, height: float, options: LewisOptions | None = None,
           measure=None) -> LewisLayout:
    """Project *structure* through *camera* and place every drawing element.

    *measure* returns the pixel extent of a run list at a font size; the Qt
    painter passes real font metrics, the default is an approximation good
    enough for tests.
    """
    options = options or LewisOptions()
    measure = measure or default_measure
    mol = structure.molecule
    pts = camera.project(mol.coords, width, height)
    hidden = hidden_hydrogens(structure, options)
    scale = _bond_scale(structure, camera, height)
    font_px = max(scale * FONT_FRACTION, 6.0)

    visible = [a for a in range(mol.natoms) if not hidden[a]]
    shown = {a: [j for j in structure.neighbors[a] if not hidden[j]] for a in visible}
    numbers = mol.atomic_numbers          # rebuilt from the symbols on each read

    atoms: dict[int, Atom2D] = {}
    for a in visible:
        z = int(numbers[a])
        atom = Atom2D(index=a, pos=pts[a], charge=int(structure.charges[a]),
                      radicals=int(structure.radicals[a]),
                      lone_pairs=int(structure.lone_pairs[a]))
        if _needs_label(structure, a, shown[a], options):
            # a lone substituent writes its hydrogens towards the free side
            flip = _flip_label(pts, a, shown[a])
            atom.runs = _label_runs(structure, a, hidden, flip)
            atom.half_width, atom.half_height = (
                v / 2.0 for v in measure(atom.runs, font_px))
            atom.gap = font_px * LABEL_GAP
        if options.color_atoms:
            atom.color = _LABEL_COLORS.get(z, elements.cpk_color(z))
        atom.charge_text = _charge_text(atom.charge)
        atoms[a] = atom

    bonds = _place_bonds(structure, atoms, hidden, shown, scale)
    _place_dots(structure, atoms, shown, options, scale, font_px)
    _place_charges(atoms, shown, font_px)

    overlaps = _count_overlaps(atoms, scale)
    warnings = []
    if structure.hydrogens_missing:
        warnings.append("no hydrogens in this file — formal charges are not "
                        "shown, since the geometry cannot say")
    if _has_metal(structure):
        warnings.append("a metal is present — its bonds are drawn as plain "
                        "lines, so the charges around it are bookkeeping")
    if not structure.matches_file:
        warnings.append(
            f"perceived charge {structure.total_charge:+d} does not match the "
            f"{int(mol.charge):+d} in the file")
    if overlaps:
        plural = "s overlap" if overlaps > 1 else " overlaps"
        warnings.append(f"{overlaps} atom pair{plural} in this view — "
                        "turn the structure")
    return LewisLayout(atoms=[atoms[a] for a in visible], bonds=bonds,
                       scale=scale, font_px=font_px,
                       line_width=max(scale * LINE_FRACTION, 0.8),
                       overlaps=overlaps, warnings=warnings,
                       brackets=_place_brackets(structure, atoms, visible,
                                                scale, font_px, measure))


def _has_metal(structure: LewisStructure) -> bool:
    """Is anything here outside the main group?

    A dative bond drawn as a plain line puts a formal charge on the donor and
    the balancing one on the metal, which is kept out of the bookkeeping — so
    the charges around a metal say more about the drawing convention than
    about the molecule, and the viewer should be told.
    """
    return any(not is_main_group(int(z))
               for z in structure.molecule.atomic_numbers)


def _bracket_runs(charge: int, radicals: int) -> tuple:
    """The label outside the brackets: a dot per unpaired electron, then the
    charge — the order a radical cation is written in."""
    runs = []
    if radicals:
        runs.append(("•" * min(radicals, 3), True))
    text = _charge_text(charge)
    if text:
        runs.append((text, True))
    return tuple(runs)


def _place_brackets(structure, atoms, visible, scale, font_px, measure):
    """Box the drawing in and hang the species charge off the top right."""
    if not structure.bracketed or not visible:
        return None
    pad = scale * BRACKET_PAD
    xs, ys = [], []
    for a in visible:                       # the text boxes, not just the dots
        atom = atoms[a]
        xs += [atom.pos[0] - atom.half_width, atom.pos[0] + atom.half_width]
        ys += [atom.pos[1] - atom.half_height, atom.pos[1] + atom.half_height]
    runs = _bracket_runs(structure.net_charge, structure.net_radicals)
    box = Brackets(left=min(xs) - pad, right=max(xs) + pad,
                   top=min(ys) - pad, bottom=max(ys) + pad,
                   arm=scale * BRACKET_ARM, runs=runs)
    width, height = measure(runs, font_px)
    box.label_pos = np.array([box.right + width / 2.0 + font_px * 0.18,
                              box.top + height / 2.0])
    return box


def _bond_scale(structure: LewisStructure, camera: OrthoCamera,
                height: float) -> float:
    """Pixels per bond, from the unprojected lengths so that turning the
    structure does not change the size of the lettering."""
    mol = structure.molecule
    px_per_angstrom = float(height) / (2.0 * max(camera.half_height, 1e-6))
    if len(structure.bonds):
        d = np.linalg.norm(mol.coords[structure.bonds[:, 0]]
                           - mol.coords[structure.bonds[:, 1]], axis=1)
        typical = float(np.median(d))
    else:
        typical = 1.45
    return px_per_angstrom * max(typical, 0.5)


def _needs_label(structure, atom: int, shown: list[int],
                 options: LewisOptions) -> bool:
    """Carbon is a bare vertex in a skeletal drawing; everything else is
    written out, and so is a carbon with nothing attached to place it."""
    if structure.molecule.symbols[atom] != "C":
        return True
    return options.carbon_labels or not shown


def _flip_label(pts, atom: int, shown: list[int]) -> bool:
    """True when the hydrogens belong on the left of the symbol (HO- rather
    than OH-), which is where they go when the bond leaves to the right."""
    if len(shown) != 1:
        return False
    return bool(pts[shown[0]][0] > pts[atom][0])


def _place_bonds(structure, atoms, hidden, shown, scale) -> list[Bond2D]:
    gap = scale * DOUBLE_GAP
    bonds: list[Bond2D] = []
    for (i, j), order in zip(structure.bonds, structure.orders):
        i, j, order = int(i), int(j), int(order)
        if hidden[i] or hidden[j]:
            continue
        a, b = atoms[i], atoms[j]
        v = b.pos - a.pos
        length = math.hypot(v[0], v[1])       # a two-vector: no numpy needed
        bond = Bond2D(i=i, j=j, order=order)
        bonds.append(bond)
        if length < 1e-6:                       # atoms exactly on top of each other
            continue
        u = v / length
        start = a.pos + u * a.clearance(u)
        end = b.pos - u * b.clearance(-u)
        if float(np.dot(end - start, u)) <= 1.0:   # the labels already touch
            continue
        perp = np.array([-u[1], u[0]])
        if order == 1:
            bond.lines = [(start, end)]
        elif order == 3:
            bond.lines = [(start, end),
                          (start + perp * gap, end + perp * gap),
                          (start - perp * gap, end - perp * gap)]
        else:
            side = _inner_side(atoms, shown, i, j, perp)
            if side is None:
                half = perp * gap / 2.0
                bond.lines = [(start + half, end + half), (start - half, end - half)]
            else:
                span = end - start
                inset = u * math.hypot(span[0], span[1]) * INNER_INSET
                offset = perp * gap * side
                bond.lines = [(start, end),
                              (start + offset + inset, end + offset - inset)]
    return bonds


def _inner_side(atoms, shown, i: int, j: int, perp: np.ndarray):
    """Which side the second line of a double bond goes on, or None for a
    symmetric pair. Inside a ring both neighbours lean the same way, and that
    is the side ChemDraw puts the inner line on."""
    lean = np.zeros(2)
    for a, other in ((i, j), (j, i)):
        neighbours = [n for n in shown[a] if n != other]
        if not neighbours:
            return None                      # a terminal double bond stays centred
        for n in neighbours:
            d = atoms[n].pos - atoms[a].pos
            norm = math.hypot(d[0], d[1])
            if norm > 1e-9:
                lean += d / norm
    projection = float(np.dot(lean, perp))
    if abs(projection) < 0.35:               # trans substituents cancel out
        return None
    return 1.0 if projection > 0 else -1.0


def _place_dots(structure, atoms, shown, options: LewisOptions,
                scale: float, font_px: float) -> None:
    """Lone pairs and radical electrons, in the emptiest directions left."""
    radius = scale * DOT_RADIUS
    for index, atom in atoms.items():
        pairs = atom.lone_pairs if options.lone_pairs else 0
        if not pairs and not atom.radicals:
            continue
        taken = []
        for n in shown[index]:
            d = atoms[n].pos - atom.pos
            norm = float(np.linalg.norm(d))
            if norm > 1e-9:
                taken.append(float(np.arctan2(d[1], d[0])))
        for slot in range(int(pairs) + int(atom.radicals)):
            angle = _free_angle(taken)
            taken.append(angle)
            direction = np.array([np.cos(angle), np.sin(angle)])
            centre = atom.pos + direction * (atom.edge(direction) + radius * 2.2)
            if slot < pairs:
                perp = np.array([-direction[1], direction[0]]) * radius * 1.6
                atom.dots.extend([centre + perp, centre - perp])
            else:
                atom.dots.append(centre)


def _free_angle(taken: list[float], samples: int = 48) -> float:
    """The direction furthest from everything already pointing out of an atom."""
    if not taken:
        return -np.pi / 2.0            # nothing in the way: straight up
    candidates = np.linspace(-np.pi, np.pi, samples, endpoint=False)
    used = np.asarray(taken)
    separation = np.abs(np.angle(np.exp(1j * (candidates[:, None] - used[None, :]))))
    return float(candidates[int(np.argmax(separation.min(axis=1)))])


def _place_charges(atoms, shown, font_px: float) -> None:
    """Charge sign at the top right of a label, or clear of the bonds on a
    bare vertex."""
    for index, atom in atoms.items():
        if not atom.charge_text:
            continue
        if atom.labelled:
            atom.charge_pos = atom.pos + np.array(
                [atom.half_width + font_px * 0.22, -atom.half_height * 1.1])
        else:
            taken = []
            for n in shown[index]:
                d = atoms[n].pos - atom.pos
                norm = float(np.linalg.norm(d))
                if norm > 1e-9:
                    taken.append(float(np.arctan2(d[1], d[0])))
            angle = _free_angle(taken)
            direction = np.array([np.cos(angle), np.sin(angle)])
            atom.charge_pos = atom.pos + direction * font_px * 0.62


def _count_overlaps(atoms, scale: float) -> int:
    """Pairs of atoms projecting on top of each other — the drawing is
    unreadable there and the fix is to turn the structure.

    Counted through a grid of squares one overlap wide rather than by
    comparing every pair: a protein puts tens of thousands of atoms on the
    canvas, and the pairwise distance matrix would be gigabytes.
    """
    pos = np.array([a.pos for a in atoms.values()])
    reach = scale * OVERLAP_FRACTION
    if len(pos) < 2 or reach <= 0:
        return 0

    cells = np.floor(pos / reach).astype(np.int64)
    buckets: dict[tuple[int, int], list[int]] = {}
    for index, cell in enumerate(map(tuple, cells)):
        buckets.setdefault(cell, []).append(index)

    total = 0
    for (cx, cy), members in buckets.items():
        near = [a for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                for a in buckets.get((cx + dx, cy + dy), ())]
        here, there = np.array(members), np.array(near)
        d = np.linalg.norm(pos[here][:, None, :] - pos[there][None, :, :], axis=2)
        total += int(((d < reach) & (here[:, None] < there[None, :])).sum())
    return total
