# Spectra

Open any frequency, TD-DFT or NMR job and the spectra panel appears; ++s++
toggles it.

![IR spectrum with the animated C–H stretch selected](../assets/screenshot_ir_spectrum.png)

- **IR** — a Lorentzian-broadened spectrum over its sticks, with an adjustable
  width and frequency scaling factor. Clicking a band
  [animates that mode](trajectories.md#vibrations).
- **UV-Vis** — Gaussian broadening in the energy domain, on a nm or eV axis,
  with oscillator strengths as sticks on their own right-hand axis.
- **NMR** — one nucleus at a time (¹H, ¹³C, …). Enter a reference shielding
  σ₀ to switch from raw shieldings to the chemical-shift scale.

Every spectrum exports as CSV for the data, or PNG, SVG or PDF for the figure.

## From the command line

`--ir`, `--uv` and `--nmr` plot the spectrum instead of the molecule, with the
same broadening and axes as the panel:

```bash
scope -s freq.log --ir -o ir.pdf --fwhm 12 --freq-scale 0.965
scope -s freq.log --ir --xrange 600:1800 --title 'fingerprint region'
scope -s td.log   --uv --unit eV --xrange 2:7 -o uv.svg
scope -s nmr.log  --nmr --nucleus C --reference 186.4 --csv shifts.csv
```

IR and chemical-shift axes run high to low, as those spectra are conventionally
drawn. Figures go to PNG, PDF, SVG, EPS or TIFF at `--dpi` (300 by default), and
`--csv` writes the plotted curve as data. Plotting a spectrum loads only
matplotlib — no Qt, no OpenGL.
