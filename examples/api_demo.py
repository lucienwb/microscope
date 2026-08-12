"""microscp Python API demo: parse Gaussian outputs and use the results.

Run from the repo root:  python examples/api_demo.py
"""

from pathlib import Path

import microscp.io as mio
from microscp.core import geometry

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"

# --- 1. geometry optimization: trajectory + energies --------------------------
result = mio.load(DATA / "dvb_gopt.out")
mol = result.molecule
print(f"[opt]  {mol.formula()}  ·  {result.nframes} frames  ·  "
      f"final E = {result.scf_energies[-1]:.6f} Ha  ·  "
      f"converged = {result.normal_termination}")

# --- 2. measurements on the final structure ----------------------------------
mol.perceive_bonds()
i, j = mol.bonds[0]
print(f"[geom] first bond {mol.symbols[i]}{i + 1}–{mol.symbols[j]}{j + 1}: "
      f"{geometry.distance(mol.coords[i], mol.coords[j]):.3f} Å  ·  "
      f"{len(mol.bonds)} bonds perceived")

# --- 3. IR frequencies --------------------------------------------------------
ir = mio.load(DATA / "dvb_ir.out")
print(f"[IR]   {len(ir.vibrations)} normal modes; 5 strongest bands:")
for v in sorted(ir.vibrations, key=lambda v: -(v.ir_intensity or 0))[:5]:
    print(f"         {v.frequency:9.1f} cm⁻¹   I = {v.ir_intensity:7.1f} km/mol")

# --- 4. TD-DFT excited states -------------------------------------------------
td = mio.load(DATA / "dvb_td.out")
print(f"[UV]   {len(td.excited_states)} excited states; lowest 3:")
for s in td.excited_states[:3]:
    print(f"         S{s.index}: {s.energy_ev:.3f} eV = {s.wavelength_nm:.1f} nm, "
          f"f = {s.osc_strength:.4f}  ({s.label})")

# --- 5. NMR shieldings --------------------------------------------------------
nmr = mio.load(DATA / "dvb_nmr.log")
carbons = [s for s in nmr.nmr_shieldings if s.symbol == "C"]
print(f"[NMR]  {len(nmr.nmr_shieldings)} shieldings; carbon isotropic values:")
print("         " + "  ".join(f"{s.isotropic:6.1f}" for s in carbons[:5]) + "  ... ppm")

# --- 6. convert / write files -------------------------------------------------
out = Path(__file__).parent
mio.save_molecule(out / "dvb_final.xyz", mol)
mio.save_molecule(out / "dvb_new_job.gjf", mol)
mio.save_molecule(out / "dvb.pdb", mol)
print(f"[save] wrote dvb_final.xyz, dvb_new_job.gjf, dvb.pdb to {out.name}/")
