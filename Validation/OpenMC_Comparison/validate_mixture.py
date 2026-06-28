"""PyNeut vs OpenMC validation for a multi-isotope material mixture.

Exercises the per-collision isotope-selection machinery (Region(composition=...))
against an independent code on a B4C-like B10+C12 sphere — a strong absorber
(B10) mixed with a moderator (C12), so getting the leakage right requires both
the aggregate macroscopic cross section and the per-collision isotope choice to
be correct. Number densities are set identically in both codes (explicit
atoms/cm^3), so any difference is transport, not the material spec or the data.

Graded by the same leakage / escape-energy z-score and spectrum chi2/dof as the
rest of the suite.

Run:  python validate_mixture.py [--n 20000]
"""
import argparse
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import (ensure_openmc_data, energy_grid, compare_spectra,  # noqa: E402
                     ratio_mean_uncertainty, N_BATCHES, _pyneut_spectrum)

sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from src.cross_section_read import CrossSectionReader  # noqa: E402
from src.material import Nuclide  # noqa: E402
from src.medium import Region, Sphere  # noqa: E402
from src.vt_calc import VelocitySampler  # noqa: E402
from src.settings import Settings  # noqa: E402
from src.random_number_generator import RNGHandler  # noqa: E402
from src.simulation import simulate_single_particle  # noqa: E402
from src.statistics import escape_statistics  # noqa: E402

ENDF = os.path.abspath(os.path.join(HERE, '..', '..', 'endfb'))

# B4C-like (all-B10), 4:1 B:C atom ratio at rho = 2.52 g/cm^3 -> explicit number
# densities used verbatim by both codes:  (element, awr, atomic_mass, N [at/cm^3])
COMPONENTS = [
    ("B10", 9.926900, 10.0129, 1.0992e23),
    ("C12", 11.896910, 12.0000, 2.7480e22),
]
RADIUS = 5.0       # cm
SOURCE_E = 1.0e6   # eV


def run_pyneut_mixture(n):
    reader = CrossSectionReader(ENDF)
    comp = [Nuclide(el, nd, awr, atomic_mass=am) for el, awr, am, nd in COMPONENTS]
    mediums = [Region(surfaces=[Sphere((0, 0, 0), RADIUS)],
                      name="mix", priority=1, composition=comp)]
    settings = Settings("shielding", n)
    rngs = [RNGHandler(12345 + i) for i in range(n)]
    states = [{
        "x": 0.0, "y": 0.0, "z": 0.0,
        "theta": np.arccos(1 - 2 * r.random()), "phi": r.uniform(0, 2 * np.pi),
        "energy": SOURCE_E, "weight": 1.0, "has_interacted": False,
    } for r in rngs]
    # the global A/N/sampler are ignored when a composition is set
    dummy = VelocitySampler(1.0)
    args = [(s, reader, mediums, 1.0, 1.0, dummy, None, False, r, settings)
            for s, r in zip(states, rngs)]

    with Pool() as pool:
        results = pool.map(simulate_single_particle, args)

    leaked = np.array([r["final_weight"] for r in results], dtype=float)
    esc = np.array([r["final_energy"] for r in results], dtype=float)
    bins = energy_grid(SOURCE_E)
    stats = escape_statistics(leaked, esc, n, N_BATCHES, e_bins=bins)
    stats["spectrum"], stats["spectrum_sem"] = _pyneut_spectrum(leaked, esc, bins, N_BATCHES)
    return stats


def run_openmc_mixture(n, batches=20):
    if not ensure_openmc_data():
        return None
    import openmc

    nds = [c[3] for c in COMPONENTS]
    ntot = sum(nds)
    mat = openmc.Material()
    for (el, _awr, _am, nd) in COMPONENTS:
        mat.add_nuclide(el, nd / ntot, "ao")
    mat.set_density("atom/b-cm", ntot * 1e-24)
    materials = openmc.Materials([mat])

    sph = openmc.Sphere(r=RADIUS, boundary_type="vacuum")
    cell = openmc.Cell(region=-sph, fill=mat)
    geometry = openmc.Geometry([cell])
    surfs = [sph]

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = batches
    settings.particles = max(1, n // batches)
    settings.output = {"summary": False, "tallies": False}
    settings.verbosity = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([SOURCE_E], [1.0]),
    )

    e_bins = energy_grid(SOURCE_E)
    e_mid = 0.5 * (e_bins[:-1] + e_bins[1:])
    t_leak = openmc.Tally(name="leakage")
    t_leak.filters = [openmc.SurfaceFilter(surfs)]
    t_leak.scores = ["current"]
    t_spec = openmc.Tally(name="spectrum")
    t_spec.filters = [openmc.SurfaceFilter(surfs), openmc.EnergyFilter(e_bins)]
    t_spec.scores = ["current"]
    tallies = openmc.Tallies([t_leak, t_spec])

    run_dir = os.path.join(HERE, "openmc_runs", "omc_mix_b4c")
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geometry, materials, settings, tallies)
        sp_path = model.run(output=False)
        if sp_path is None:
            import glob
            sp_path = sorted(glob.glob(os.path.join(run_dir, "statepoint.*.h5")),
                             key=os.path.getmtime)[-1]
        with openmc.StatePoint(sp_path) as sp:
            lt = sp.get_tally(name="leakage")
            leak = float(np.abs(lt.mean.flatten()).sum())
            leak_sem = float(np.sqrt((lt.std_dev.flatten() ** 2).sum()))
            st = sp.get_tally(name="spectrum")
            mean = np.abs(st.mean.flatten()).reshape(len(surfs), len(e_mid))
            sd = st.std_dev.flatten().reshape(len(surfs), len(e_mid))
            bin_w = mean.sum(axis=0)
            bin_sd = np.sqrt((sd ** 2).sum(axis=0))
            avg_e, avg_e_sem = ratio_mean_uncertainty(e_mid, bin_w, bin_sd)
    finally:
        os.chdir(cwd)
    return {"leakage": leak, "leakage_sem": leak_sem,
            "avg_energy": avg_e, "avg_energy_sem": avg_e_sem,
            "spectrum": bin_w, "spectrum_sem": bin_sd}


def _z(a, sa, b, sb):
    comb = np.hypot(sa, sb)
    return abs(a - b) / comb if comb > 0 else float("inf")


def main():
    ap = argparse.ArgumentParser(description="PyNeut vs OpenMC mixture validation")
    ap.add_argument("--n", type=int, default=20000, help="histories")
    args = ap.parse_args()

    iso = ", ".join(f"{el} ({nd:.3e}/cm³)" for el, _, _, nd in COMPONENTS)
    print(f"\nB4C-like mixture sphere (R={RADIUS} cm, {SOURCE_E/1e6:g} MeV source)\n  {iso}\n")

    pyn = run_pyneut_mixture(args.n)
    omc = run_openmc_mixture(args.n)
    if omc is None:
        print(f"PyNeut leakage = {pyn['leakage']:.5f} ± {pyn['leakage_sem']:.5f}  (OpenMC unavailable)")
        return

    z_leak = _z(omc["leakage"], omc["leakage_sem"], pyn["leakage"], pyn["leakage_sem"])
    z_e = _z(omc["avg_energy"], omc["avg_energy_sem"], pyn["avg_energy"], pyn["avg_energy_sem"])
    _, _, chi2_dof = compare_spectra(omc, pyn)

    print(f"{'metric':<16}{'OpenMC':>16}{'PyNeut':>16}{'z':>7}")
    print("-" * 55)
    print(f"{'leakage':<16}{omc['leakage']:>10.5f} ± {omc['leakage_sem']:.5f}"
          f"{'':>0}  {pyn['leakage']:>8.5f} ± {pyn['leakage_sem']:.5f}{z_leak:>7.1f}")
    print(f"{'avg energy (eV)':<16}{omc['avg_energy']:>10.0f}        "
          f"  {pyn['avg_energy']:>8.0f}        {z_e:>7.1f}")
    integral = "OK" if max(z_leak, z_e) <= 2 else ("WARNING" if max(z_leak, z_e) <= 4 else "FAIL")
    spec = "GOOD" if chi2_dof <= 2 else ("FAIR" if chi2_dof <= 5 else "POOR")
    print(f"\nintegral: {integral}   spectrum: {spec} (χ²/dof = {chi2_dof:.2f})\n")


if __name__ == "__main__":
    main()
