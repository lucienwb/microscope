"""`scope -s file.log --ir/--uv/--nmr`: plot a spectrum instead of a molecule.

matplotlib only — no Qt, no OpenGL.
"""

from __future__ import annotations

from pathlib import Path

from .. import io as mio
from .parsing import (
    DEFAULT_FWHM,
    CliError,
    check_output,
    default_output,
    parse_range,
    resolve_size,
    spectrum_kind,
)


def pick_nucleus(shieldings, wanted: str | None) -> str:
    """The element to plot: what was asked for, else H, else C, else any."""
    available = sorted({s.symbol for s in shieldings},
                       key=lambda e: {"H": 0, "C": 1}.get(e, 2))
    if not available:
        raise CliError("no NMR shieldings in this file")
    if wanted is None:
        return available[0]
    if wanted not in available:
        raise CliError(f"no {wanted} shieldings in this file "
                       f"(it has {', '.join(available)})")
    return wanted


def build_spectrum(kind: str, result, args):
    """Turn the parsed result into a Spectrum, with the flags applied."""
    from ..spectra import ir_spectrum, nmr_spectrum, uvvis_spectrum

    fwhm = args.fwhm if args.fwhm is not None else DEFAULT_FWHM[kind]
    if fwhm <= 0:
        raise CliError("--fwhm must be positive")
    try:
        if kind == "ir":
            if not result.vibrations:
                raise CliError("no vibrational frequencies in this file "
                               "(it needs a freq job)")
            return ir_spectrum(result.vibrations, scale=args.freq_scale,
                               fwhm=fwhm)
        if kind == "uv":
            if not result.excited_states:
                raise CliError("no excited states in this file "
                               "(it needs a TD-DFT/CIS job)")
            return uvvis_spectrum(result.excited_states, fwhm_ev=fwhm,
                                  unit=args.unit)
        nucleus = pick_nucleus(result.nmr_shieldings, args.nucleus)
        return nmr_spectrum(result.nmr_shieldings, element=nucleus,
                            reference=args.reference, fwhm=fwhm)
    except ValueError as exc:            # the builders' own complaints
        raise CliError(str(exc)) from None


def plot_spectrum(args) -> str:
    """Silent mode with --ir/--uv/--nmr: write the spectrum figure."""
    kind = spectrum_kind(args)
    if not args.file:
        raise CliError(f"give a file to plot, e.g. scope -s --{kind} mycalc.log")
    if not Path(args.file).is_file():
        raise FileNotFoundError(args.file)
    out = check_output(args.output or default_output(args.file, kind), kind)

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    from ..spectra.plot import draw_spectrum, set_xrange

    result = mio.load(args.file)
    spec = build_spectrum(kind, result, args)

    width, height = resolve_size(args, kind)
    dpi = max(int(args.dpi), 1)
    figure = Figure(figsize=(width / dpi, height / dpi), dpi=dpi,
                    tight_layout=True)
    FigureCanvasAgg(figure)              # savefig needs a canvas to start from
    ax = figure.add_subplot(111)
    # UV-Vis oscillator strengths get their own axis, as in the viewer
    stick_axis = ax.twinx() if kind == "uv" else None
    draw_spectrum(ax, spec, sticks=not args.no_sticks, stick_axis=stick_axis,
                  stick_label="oscillator strength f")
    if args.xrange:
        set_xrange(ax, spec, *parse_range(args.xrange))
    if args.title:
        ax.set_title(args.title)

    figure.savefig(out, transparent=not args.opaque)
    if args.csv:
        import pandas as pd
        pd.DataFrame({spec.xlabel: spec.x, spec.ylabel: spec.y}).to_csv(
            args.csv, index=False)
    return out
