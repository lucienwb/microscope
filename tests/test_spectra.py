"""Spectra engine tests."""

import numpy as np
import pytest

from microscp.core.results import ExcitedState, NMRShielding, Vibration
from microscp.spectra import broaden, ir_spectrum, nmr_spectrum, uvvis_spectrum


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
