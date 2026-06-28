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
from src.material import Nuclide  # noqa: E402
from src.criticality import run_keff, uniform_sphere_source  # noqa: E402

# Bare U235 metal: density and atomic-weight ratio (ENDF/B-VIII).
U235_RHO, U235_M, U235_A = 18.74, 235.0, 233.0248

# Bare U235 critical radius is ~8.7 cm; bracket it.
RADII = [6.0, 8.7, 11.0]

# ICSBEP bare fast-metal critical benchmarks (k_eff ≡ 1.0000 by construction).
# comp entries are (element, number_density [atoms/b-cm], awr, atomic_mass).
BENCHMARKS = {
    "godiva": {
        "name": "Godiva (HEU-MET-FAST-001)", "radius": 8.7407,
        "source": "U235", "k": (1.0000, 0.0010),
        "comp": [
            ("U234", 4.8271e-04, 232.030, 234.04),
            ("U235", 4.4994e-02, 233.025, 235.04),
            ("U238", 2.4984e-03, 236.006, 238.05),
        ],
    },
    "jezebel23": {
        "name": "Jezebel-23 (U233-MET-FAST-001)", "radius": 5.9838,
        "source": "U233", "k": (1.0000, 0.0010),
        "comp": [
            ("U233", 4.6655e-02, 231.038, 233.04),
            ("U234", 5.9170e-04, 232.030, 234.04),
            ("U235", 1.4250e-04, 233.025, 235.04),
            ("U238", 2.9280e-04, 236.006, 238.05),
        ],
    },
    "jezebel": {
        "name": "Jezebel (PU-MET-FAST-001)", "radius": 6.3849,
        "source": "Pu239", "k": (1.0000, 0.0020),
        "comp": [
            ("Pu239", 3.7047e-02, 236.999, 239.05),
            ("Pu240", 1.7512e-03, 237.992, 240.05),
            ("Pu241", 1.1674e-04, 238.985, 241.06),
            ("Ga69",  8.2663e-04, 68.380, 68.93),
            ("Ga71",  5.4857e-04, 70.360, 70.92),
        ],
    },
}


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


def _composition(spec):
    return [Nuclide(el, nd * 1e24, awr, atomic_mass=am)
            for (el, nd, awr, am) in spec["comp"]]


def benchmark_pyneut(spec, n, gens, inactive, seed=1):
    reader = CrossSectionReader(ENDF_PATH)
    R = spec["radius"]
    mediums = [Region(surfaces=[Sphere((0, 0, 0), R)],
                      name="core", priority=1, composition=_composition(spec))]
    src = uniform_sphere_source((0, 0, 0), R, n, RNGHandler(seed),
                                element=spec["source"])
    settings = Settings("criticality", 1)
    from multiprocessing import Pool
    with Pool() as pool:
        out = run_keff(reader, mediums, src, 1.0, 1.0, VelocitySampler(1.0),
                       settings, generations=gens, inactive=inactive,
                       seed=seed, pool=pool)
    return out["k_eff"], out["k_sem"]


def benchmark_openmc(key, spec, n, gens, inactive):
    if not ensure_openmc_data():
        return None, None
    import openmc

    ntot = sum(c[1] for c in spec["comp"])
    mat = openmc.Material()
    for (el, nd, _awr, _am) in spec["comp"]:
        mat.add_nuclide(el, nd / ntot, "ao")
    mat.set_density("atom/b-cm", ntot)
    materials = openmc.Materials([mat])

    R = spec["radius"]
    sph = openmc.Sphere(r=R, boundary_type="vacuum")
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
        space=openmc.stats.Box((-R, -R, -R), (R, R, R)),
        constraints={"fissionable": True})

    run_dir = os.path.join(HERE, "openmc_runs", f"omc_bench_{key}")
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geometry, materials, settings)
        sp_path = model.run(output=False)
        if sp_path is None:
            import glob
            sp_path = sorted(glob.glob(os.path.join(run_dir, "statepoint.*.h5")),
                             key=os.path.getmtime)[-1]
        with openmc.StatePoint(sp_path) as sp:
            k = sp.keff
        return float(k.nominal_value), float(k.std_dev)
    finally:
        os.chdir(cwd)


def run_benchmarks(keys, n, gens, inactive):
    print(f"\nICSBEP fast-metal criticality benchmarks "
          f"(N={n}/gen, {gens} gens, {inactive} inactive)\n")
    hdr = (f"{'benchmark':<30}{'PyNeut k_eff':>18}{'OpenMC k_eff':>18}"
           f"{'ICSBEP':>10}{'z(OMC)':>8}{'z(bench)':>9}")
    print(hdr)
    print("-" * len(hdr))
    for key in keys:
        spec = BENCHMARKS[key]
        missing = [c[0] for c in spec["comp"]
                   if not os.path.exists(os.path.join(ENDF_PATH, "neutron", f"{c[0]}.h5"))]
        if missing:
            print(f"{spec['name']:<30}  SKIPPED — needs {', '.join(f'{m}.h5' for m in missing)}")
            continue
        kp, sp = benchmark_pyneut(spec, n, gens, inactive)
        ko, so = benchmark_openmc(key, spec, n, gens, inactive)
        kb, sb = spec["k"]
        z_bench = abs(kp - kb) / np.hypot(sp, sb)
        if ko is None:
            print(f"{spec['name']:<30}{kp:>10.4f} ± {sp:.4f}{'(no openmc)':>18}"
                  f"{kb:>10.4f}{'-':>8}{z_bench:>9.1f}")
            continue
        z_omc = abs(kp - ko) / np.hypot(sp, so) if np.hypot(sp, so) > 0 else float("inf")
        print(f"{spec['name']:<30}{kp:>10.4f} ± {sp:.4f}{ko:>10.4f} ± {so:.4f}"
              f"{kb:>10.4f}{z_omc:>8.1f}{z_bench:>9.1f}")
    print("\n(z = |Δk| / combined σ; ICSBEP k_eff ≡ 1.0000 by benchmark construction)\n")


def main():
    ap = argparse.ArgumentParser(description="PyNeut vs OpenMC k-eff (U235 spheres)")
    ap.add_argument("--n", type=int, default=2000, help="neutrons per generation")
    ap.add_argument("--gens", type=int, default=120, help="total generations")
    ap.add_argument("--inactive", type=int, default=20, help="inactive generations")
    ap.add_argument("--radii", type=str, default="", help="comma list of radii (cm)")
    ap.add_argument("--benchmarks", nargs="?", const="all", default=None,
                    help="run ICSBEP benchmarks (godiva,jezebel23,jezebel; default all)")
    ap.add_argument("--godiva", action="store_true", help="alias for --benchmarks godiva")
    args = ap.parse_args()

    if args.godiva:
        keys = ["godiva"]
    elif args.benchmarks is not None:
        keys = (list(BENCHMARKS) if args.benchmarks in ("all", "")
                else args.benchmarks.split(","))
    else:
        keys = None

    if keys:
        run_benchmarks(keys, args.n, args.gens, args.inactive)
        return

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
