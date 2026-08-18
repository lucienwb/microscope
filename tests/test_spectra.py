"""Spectra engine tests."""

import numpy as np
import pytest

from microscope.core.results import ExcitedState, NMRShielding, Vibration
from microscope.spectra import broaden, ir_spectrum, nmr_spectrum, uvvis_spectrum


def test_broaden_peak_height():
    grid = np.linspace(0, 100, 10001)
    y = broaden([50.0], [2.0], grid, fwhm=4.0, shape="lorentzian")
    assert abs(y.max() - 2.0) < 1e-3          # height preserved at center
    assert abs(grid[np.argmax(y)] - 50.0) < 0.02
    yg = broaden([50.0], [3.0], grid, fwhm=4.0, shape="gaussian")
    assert abs(yg.max() - 3.0) < 1e-3
    # half maximum at +/- fwhm/2
    at_half = np.interp(52.0, grid, yg)
    assert abs(at_half - 1.5) < 0.01


def test_ir_spectrum():
    vibs = [Vibration(frequency=-50.0, ir_intensity=99.0),   # imaginary: excluded
            Vibration(frequency=1000.0, ir_intensity=10.0),
            Vibration(frequency=1700.0, ir_intensity=50.0)]
    spec = ir_spectrum(vibs, scale=0.97, fwhm=8.0)
    assert spec.invert_x
    assert len(spec.stick_x) == 2
    assert abs(spec.stick_x[1] - 1700.0 * 0.97) < 1e-9
    peak_x = spec.x[np.argmax(spec.y)]
    assert abs(peak_x - 1700.0 * 0.97) < 2.0


def test_uvvis_units():
    states = [ExcitedState(1, "S", 4.0, 309.96, 0.5),
              ExcitedState(2, "S", 5.0, 247.97, 0.1)]
    nm = uvvis_spectrum(states, unit="nm")
    ev = uvvis_spectrum(states, unit="eV")
    assert abs(nm.stick_x[0] - 309.96) < 0.01
    assert abs(ev.stick_x[0] - 4.0) < 1e-9
    assert np.all(np.diff(nm.x) >= 0)          # ascending nm axis
    peak_ev = ev.x[np.argmax(ev.y)]
    assert abs(peak_ev - 4.0) < 0.05


def test_nmr_reference_switch():
    sh = [NMRShielding(0, "C", 120.0), NMRShielding(1, "C", 60.0),
          NMRShielding(2, "H", 25.0)]
    raw = nmr_spectrum(sh, element="C", reference=0.0, fwhm=1.0)
    assert not raw.invert_x
    assert set(np.round(raw.stick_x, 6)) == {120.0, 60.0}
    delta = nmr_spectrum(sh, element="C", reference=185.0, fwhm=1.0)
    assert delta.invert_x
    assert set(np.round(delta.stick_x, 6)) == {65.0, 125.0}
    with pytest.raises(ValueError):
        nmr_spectrum(sh, element="N")


# ---------------------------------------------------------------- plotting


def _axes():
    """A bare Agg figure — no display, no pyplot, no Qt."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    figure = Figure()
    FigureCanvasAgg(figure)
    return figure, figure.add_subplot(111)


def test_draw_spectrum_shows_the_whole_curve():
    """Regression: pinning the y bottom to 0 switches autoscaling off, so the
    top has to be set from the data or tall peaks get clipped away."""
    from microscope.spectra.plot import draw_spectrum
    # unit-height sticks whose broadened peaks add up well above 1
    sh = [NMRShielding(i, "C", 100.0 + 0.1 * i) for i in range(5)]
    spec = nmr_spectrum(sh, element="C", reference=0.0, fwhm=1.0)
    _, ax = _axes()
    draw_spectrum(ax, spec)
    assert spec.y.max() > 1.0                      # the case that used to clip
    assert ax.get_ylim()[1] >= spec.y.max()
    assert ax.get_ylim()[0] == 0.0


def test_draw_spectrum_puts_uvvis_sticks_on_a_right_hand_axis():
    from microscope.spectra.plot import draw_spectrum
    states = [ExcitedState(1, "S1", 4.0, 309.96, 0.5)]
    spec = uvvis_spectrum(states)
    _, ax = _axes()
    ax2 = ax.twinx()
    ax2.clear()                                    # as the dock does on redraw
    draw_spectrum(ax, spec, stick_axis=ax2, stick_label="oscillator strength f")
    assert ax2.yaxis.get_label_position() == "right"
    assert ax2.get_ylabel() == "oscillator strength f"
    assert ax2.get_ylim()[1] >= 0.5


def test_draw_spectrum_can_hide_the_sticks():
    from microscope.spectra.plot import draw_spectrum
    vibs = [Vibration(1600.0, ir_intensity=30.0)]
    spec = ir_spectrum(vibs, fwhm=10.0)
    _, ax = _axes()
    draw_spectrum(ax, spec, sticks=False)
    assert not ax.collections                      # vlines would land here


def test_ir_axis_runs_from_high_to_low_wavenumber():
    from microscope.spectra.plot import style_axes
    spec = ir_spectrum([Vibration(1600.0, ir_intensity=30.0)], fwhm=10.0)
    _, ax = _axes()
    style_axes(ax, spec)
    lo, hi = ax.get_xlim()
    assert lo > hi                                 # inverted, as spectra are drawn


def test_set_xrange_keeps_the_inverted_direction_and_rescales_y():
    from microscope.spectra.plot import draw_spectrum, set_xrange
    vibs = [Vibration(600.0, ir_intensity=100.0), Vibration(1600.0, ir_intensity=5.0)]
    spec = ir_spectrum(vibs, fwhm=8.0)
    _, ax = _axes()
    draw_spectrum(ax, spec)
    set_xrange(ax, spec, 1400, 1800)
    lo, hi = ax.get_xlim()
    assert (lo, hi) == (1800.0, 1400.0)            # still high -> low
    inside = (spec.x >= 1400) & (spec.x <= 1800)
    assert ax.get_ylim()[1] < spec.y.max()         # tall 600 cm-1 band excluded
    assert ax.get_ylim()[1] >= spec.y[inside].max()
