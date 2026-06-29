"""Runtime benchmark: PyNeut vs OpenMC on the BEAVRS matched-physics pin cell.

Runs the identical k-infinity problem (free-gas, 294 K, all-reflective lattice)
through both codes with the same particle budget on the same core count, and
reports wall-clock, neutron throughput, and the slowdown ratio. This quantifies
the cost of a pure-Python kernel on a thermal problem -- the worst case, since a
reflective thermal lattice maximises collisions per history.

Usage: python benchmark_pincell_runtime.py [--n 600 --gens 40 --inactive 15]
"""
import argparse
import glob
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))

import validate_pincell as VP  # noqa: E402
from src.cross_section_read import CrossSectionReader  # noqa: E402
from src.vt_calc import VelocitySampler  # noqa: E402
from src.settings import Settings  # noqa: E402
from src.random_number_generator import RNGHandler  # noqa: E402
from src.criticality import run_keff  # noqa: E402


def time_pyneut(endf, comps, n, gens, inactive, ncores):
    reader = CrossSectionReader(endf, temperature="294K")
    mediums = VP.pyneut_mediums(comps)
    src = VP.pellet_source(n, RNGHandler(1))
    settings = Settings("criticality", 1)
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(ncores) as pool:
        out = run_keff(reader, mediums, src, 1.0, 1.0, VelocitySampler(1.0),
                       settings, generations=gens, inactive=inactive, seed=1, pool=pool)
    wall = time.time() - t0
    return out["k_eff"], out["k_sem"], wall


def time_openmc(materials, n, gens, inactive):
    import openmc
    half = VP.PITCH / 2.0
    mats = openmc.Materials(list(materials.values()))
    pellet = openmc.ZCylinder(r=VP.R_PELLET)
    inner = openmc.ZCylinder(r=VP.R_INNER_CLAD)
    outer = openmc.ZCylinder(r=VP.R_OUTER_CLAD)
    xl = openmc.XPlane(-half, boundary_type='reflective'); xr = openmc.XPlane(half, boundary_type='reflective')
    yb = openmc.YPlane(-half, boundary_type='reflective'); yt = openmc.YPlane(half, boundary_type='reflective')
    zb = openmc.ZPlane(-half, boundary_type='reflective'); zt = openmc.ZPlane(half, boundary_type='reflective')
    box = +xl & -xr & +yb & -yt & +zb & -zt
    cells = [openmc.Cell(fill=materials['fuel'], region=-pellet & box),
             openmc.Cell(fill=materials['gap'], region=+pellet & -inner & box),
             openmc.Cell(fill=materials['clad'], region=+inner & -outer & box),
             openmc.Cell(fill=materials['water'], region=+outer & box)]
    geom = openmc.Geometry(cells)
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = gens; settings.inactive = inactive; settings.particles = n
    settings.temperature = {'method': 'nearest', 'default': 294.0}
    settings.output = {'summary': False, 'tallies': False}; settings.verbosity = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-VP.R_PELLET, -VP.R_PELLET, -half), (VP.R_PELLET, VP.R_PELLET, half)),
        constraints={'fissionable': True})
    run_dir = os.path.join(HERE, 'openmc_runs', 'omc_pincell_bench')
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geom, mats, settings)
        t0 = time.time()
        sp_path = model.run(output=False)
        wall = time.time() - t0
        if sp_path is None:
            sp_path = sorted(glob.glob(os.path.join(run_dir, 'statepoint.*.h5')), key=os.path.getmtime)[-1]
        with openmc.StatePoint(sp_path) as sp:
            k = sp.keff
        return float(k.nominal_value), float(k.std_dev), wall
    finally:
        os.chdir(cwd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endf", default=os.environ.get("PYNEUT_PINCELL_ENDF", VP.DEFAULT_ENDF))
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--gens", type=int, default=40)
    ap.add_argument("--inactive", type=int, default=15)
    ap.add_argument("--cores", type=int, default=os.cpu_count())
    args = ap.parse_args()

    os.environ["OMP_NUM_THREADS"] = str(args.cores)
    xs_xml = os.path.join(args.endf, "cross_sections.xml")
    if os.path.exists(xs_xml):
        os.environ["OPENMC_CROSS_SECTIONS"] = xs_xml

    histories = args.n * args.gens
    print(f"\nBEAVRS matched-physics pin cell (free-gas, 294 K, reflective lattice)")
    print(f"budget: {args.n} neutrons/gen x {args.gens} gens = {histories:,} histories "
          f"({args.inactive} inactive)  |  {args.cores} cores\n")

    mats = VP.build_openmc_materials(temperature=294.0, thermal=False)
    comps = {k: VP.pyneut_composition(m, 294.0) for k, m in mats.items()}

    kp, sp, wp = time_pyneut(args.endf, comps, args.n, args.gens, args.inactive, args.cores)
    mats = VP.build_openmc_materials(temperature=294.0, thermal=False)
    ko, so, wo = time_openmc(mats, args.n, args.gens, args.inactive)

    thr_p = histories / wp
    thr_o = histories / wo
    z = abs(ko - kp) / np.hypot(sp, so) if np.hypot(sp, so) > 0 else float("nan")
    print(f"{'':12}{'k_inf':>18}{'wall (s)':>12}{'throughput (n/s)':>20}{'per core (n/s)':>16}")
    print("-" * 78)
    print(f"{'PyNeut':12}{kp:>10.5f} ± {sp:.5f}{wp:>12.1f}{thr_p:>20,.0f}{thr_p/args.cores:>16,.0f}")
    print(f"{'OpenMC':12}{ko:>10.5f} ± {so:.5f}{wo:>12.1f}{thr_o:>20,.0f}{thr_o/args.cores:>16,.0f}")
    print("-" * 78)
    print(f"\nk_inf agreement: z = {z:.2f}")
    print(f"PyNeut is {wp / wo:,.0f}x slower than OpenMC on this thermal pin cell "
          f"(both on {args.cores} cores).")


if __name__ == "__main__":
    main()
