"""Drawing a Spectrum onto a matplotlib axes.

Shared by the interactive spectra dock and the `scope -s --ir/--uv/--nmr`
command line, so a batch figure looks exactly like the one on screen. No Qt
and no pyplot here: callers bring their own Figure.
"""

from __future__ import annotations

import numpy as np

from .models import Spectrum

CURVE_COLOR = "#303030"
STICK_COLOR = "#9a9a9a"
ACTIVE_COLOR = "#e08214"


def style_axes(ax, spec: Spectrum) -> None:
    """Axis labels, the reversed axis of IR/δ scales, and a clean frame."""
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    if spec.invert_x:
        ax.set_xlim(float(spec.x.max()), float(spec.x.min()))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _top(*arrays) -> float:
    """Head-room above the tallest thing drawn on an axis."""
    peak = max((float(np.max(a)) for a in arrays if len(a)), default=0.0)
    return peak * 1.06 if peak > 0 else 1.0


def draw_spectrum(ax, spec: Spectrum, sticks: bool = True,
                  stick_axis=None, stick_label: str | None = None) -> None:
    """Broadened curve plus the underlying stick spectrum.

    *stick_axis* puts the sticks on their own y axis (UV-Vis draws oscillator
    strengths that way, since they do not share the absorbance scale).

    Both y limits are set explicitly from the data: pinning the bottom to zero
    switches autoscaling off, so an inferred top would freeze at whatever had
    been drawn first and clip everything taller.
    """
    drawn = [spec.y]
    if sticks and len(spec.stick_x):
        bars = stick_axis if stick_axis is not None else ax
        bars.vlines(spec.stick_x, 0.0, spec.stick_y, color=STICK_COLOR, lw=1.0)
        if stick_axis is None:
            drawn.append(spec.stick_y)
        else:
            stick_axis.set_ylim(0.0, _top(spec.stick_y))
            stick_axis.spines["top"].set_visible(False)
            # clear() on a twin axis resets its ticks and label to the left
            stick_axis.yaxis.tick_right()
            stick_axis.yaxis.set_label_position("right")
            if stick_label:
                stick_axis.set_ylabel(stick_label, color=STICK_COLOR)
                stick_axis.tick_params(axis="y", labelcolor=STICK_COLOR)
    ax.plot(spec.x, spec.y, color=CURVE_COLOR, lw=1.3)
    ax.set_ylim(0.0, _top(*drawn))
    style_axes(ax, spec)


def set_xrange(ax, spec: Spectrum, lo: float, hi: float) -> None:
    """Limit the x axis, keeping the reversed direction where it applies."""
    lo, hi = float(min(lo, hi)), float(max(lo, hi))
    ax.set_xlim((hi, lo) if spec.invert_x else (lo, hi))
    inside = (spec.x >= lo) & (spec.x <= hi)
    if inside.any():                     # rescale y to what is actually shown
        top = float(np.max(spec.y[inside]))
        ax.set_ylim(0.0, top * 1.06 if top > 0 else 1.0)
