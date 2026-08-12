# Changelog

## 0.1.0 — 2026-08-12

Initial release.

- **Viewer**: CYLview-style OpenGL renderer (ray-cast impostor spheres and
  cylinders, split-color bonds, orthographic camera), atom picking,
  distance/angle/dihedral measurements with pinnable in-scene labels,
  zoom-adaptive atom labels, view alignment (along bond / plane to screen),
  rotation centering on an atom, high-resolution transparent PNG export
- **File formats**: native parsers for Gaussian (`.log`/`.out`, `.gjf`,
  `.fchk`), ORCA, Q-Chem, Molden, XYZ (multi-frame) and PDB; content-based
  program detection; XYZ/GJF/PDB writers; optional cclib fallback
- **Spectra**: IR (Lorentzian broadening, frequency scaling, click a band to
  animate the vibrational mode in 3D), UV-Vis (Gaussian broadening, nm/eV),
  NMR (per-nucleus, σ→δ referencing); CSV and PNG/SVG/PDF export
- **Editing**: fragment-aware distance/angle/dihedral adjustment with live
  preview, atom deletion, bond re-perception, 100-level undo/redo
- **Trajectories**: frame slider with SCF energies and playback button
