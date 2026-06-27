"""PyNeut vs OpenMC k-eigenvalue validation for bare U235 metal spheres.

Mirrors the rest of the suite: the OpenMC reference is generated live from the
same ENDF/B HDF5 data with matched geometry and material, and agreement is graded
by the z-score z = |k_OMC - k_PyN| / sqrt(sigma_OMC^2 + sigma_PyN^2).

Both codes run analog criticality. PyNeut uses a Watt fission spectrum (chi) and
the tabulated total nu-bar; OpenMC uses the full ENDF chi, so a small modelling
difference in the spectrum is expected on top of the statistical error.

Run:  python validate_keff.py [--n 2000] [--gens 120] [--inactive 20]
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from _common import ensure_openmc_data, ENDF_PATH, HAS_OPENMC  # noqa: E402

sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from src.cross_section_read import CrossSectionReader  # noqa: E402
from src.material import Material  # noqa: E402
from src.medium import Region, Sphere  # noqa: E402
from src.vt_calc import VelocitySampler  # noqa: E402
from src.settings import Settings  # noqa: E402
from src.random_number_generator import RNGHandler  # noqa: E402
from src.criticality import run_keff, uniform_sphere_source  # noqa: E402

# Bare U235 metal: density and atomic-weight ratio (ENDF/B-VIII).
U235_RHO, U235_M, U235_A = 18.74, 235.0, 233.0248

# Bare U235 critical radius is ~8.7 cm; bracket it.
RADII = [6.0, 8.7, 11.0]


def pyneut_keff(radius, n, gens, inactive, seed=1):
    reader = CrossSectionReader(ENDF_PATH)
    mat = Material("U235", U235_RHO, U235_M, U235_A)
    sampler = VelocitySampler(mat.kg_mass)
    mediums = [Region(surfaces=[Sphere((0, 0, 0), radius)],
                      name="core", priority=1, element="U235")]
    src = uniform_sphere_source((0, 0, 0), radius, n, RNGHandler(seed))
    settings = Settings("criticality", 1)

    from multiprocessing import Pool
    with Pool() as pool:
        out = run_keff(reader, mediums, src, U235_A, mat.number_density, sampler,
                       settings, generations=gens, inactive=inactive,
                       seed=seed, pool=pool)
    return out["k_eff"], out["k_sem"]


def openmc_keff(radius, n, gens, inactive):
    if not ensure_openmc_data():
        return None, None
    import openmc

    mat = openmc.Material()
    mat.add_nuclide("U235", 1.0)
    mat.set_density("g/cm3", U235_RHO)
    materials = openmc.Materials([mat])

    sph = openmc.Sphere(r=radius, boundary_type="vacuum")
    cell = openmc.Cell(region=-sph, fill=mat)
    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = gens
    settings.inactive = inactive
    settings.particles = n
    settings.output = {"summary": False, "tallies": False}
    settings.verbosity = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-radius, -radius, -radius),
                               (radius, radius, radius)),
        constraints={"fissionable": True},
    )

    run_dir = os.path.join(HERE, "openmc_runs", f"omc_keff_u235_r{radius:g}")
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geometry, materials, settings)
        sp_path = model.run(output=False)
        if sp_path is None:
            import glob
            sps = sorted(glob.glob(os.path.join(run_dir, "statepoint.*.h5")),
                         key=os.path.getmtime)
            sp_path = sps[-1]
        with openmc.StatePoint(sp_path) as sp:
            k = sp.keff
        return float(k.nominal_value), float(k.std_dev)
    finally:
        os.chdir(cwd)


def main():
    ap = argparse.ArgumentParser(description="PyNeut vs OpenMC k-eff (U235 spheres)")
    ap.add_argument("--n", type=int, default=2000, help="neutrons per generation")
    ap.add_argument("--gens", type=int, default=120, help="total generations")
    ap.add_argument("--inactive", type=int, default=20, help="inactive generations")
    ap.add_argument("--radii", type=str, default="", help="comma list of radii (cm)")
    args = ap.parse_args()

    radii = [float(x) for x in args.radii.split(",")] if args.radii else RADII

    print(f"\nPyNeut vs OpenMC k-eigenvalue — bare U235 spheres "
          f"(N={args.n}/gen, {args.gens} gens, {args.inactive} inactive)\n")
    header = f"{'R (cm)':>7} | {'PyNeut k_eff':>20} | {'OpenMC k_eff':>20} | {'z':>5} | status"
    print(header)
    print("-" * len(header))

    worst = 0.0
    for R in radii:
        k_p, s_p = pyneut_keff(R, args.n, args.gens, args.inactive)
        k_o, s_o = openmc_keff(R, args.n, args.gens, args.inactive)
        if k_o is None:
            print(f"{R:7.2f} | {k_p:9.5f} ± {s_p:7.5f} | {'(no openmc)':>20} |   -   | NO_REF")
            continue
        denom = np.sqrt(s_p ** 2 + s_o ** 2)
        z = abs(k_o - k_p) / denom if denom > 0 else float("inf")
        worst = max(worst, z)
        status = "OK" if z <= 2 else ("WARNING" if z <= 4 else "FAIL")
        print(f"{R:7.2f} | {k_p:9.5f} ± {s_p:7.5f} | {k_o:9.5f} ± {s_o:7.5f} | {z:5.1f} | {status}")

    print(f"\nworst z = {worst:.1f}  (OK ≤ 2σ, WARNING ≤ 4σ, FAIL > 4σ)\n")


if __name__ == "__main__":
    main()
