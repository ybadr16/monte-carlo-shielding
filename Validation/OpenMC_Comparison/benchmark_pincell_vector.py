"""BEAVRS HZP pin-cell k-infinity on the vector engine + throughput.

Reuses validate_pincell.py's materials (OpenMC-expanded atom densities fed
unchanged to PyNeut), geometry and pellet source law, but transports each
generation through the event-based vector engine (run_keff_vector). Reports
k_inf (collision estimator), the per-generation throughput, and the
comparison against the committed 5-seed references in
pincell_reference.json (PyNeut scalar 1.2737 +/- 0.0017, OpenMC
1.2890 +/- 0.0004, matched physics).

Run:  python benchmark_pincell_vector.py [--n 3000 --gens 140 --inactive 40]
      python benchmark_pincell_vector.py --quick
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))

from _common import ENDF_PATH  # noqa: E402
from validate_pincell import (  # noqa: E402
    PITCH, R_PELLET, build_openmc_materials, pyneut_composition,
    pyneut_mediums,
)
from src.cross_section_read import CrossSectionReader  # noqa: E402
from src.settings import Settings                      # noqa: E402
from src.physics import watt_params_for               # noqa: E402
from src import vector_transport as vt                # noqa: E402

REF_JSON = os.path.join(HERE, 'pincell_reference.json')


def pellet_source_arrays(n, rng):
    """Vector form of validate_pincell.pellet_source: uniform in the pellet,
    isotropic, Watt energies."""
    half = PITCH / 2.0
    r = R_PELLET * np.sqrt(rng.random(n))
    ang = 2.0 * np.pi * rng.random(n)
    mu = 2.0 * rng.random(n) - 1.0
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(1.0 - mu * mu)
    a_w, b_w = watt_params_for('U235')
    return {
        'x': r * np.cos(ang), 'y': r * np.sin(ang),
        'z': rng.uniform(-half, half, n),
        'u': s * np.cos(phi), 'v': s * np.sin(phi), 'w': mu,
        'energy': vt.sample_watt_many(rng, n, a_w, b_w),
    }


def _compare_to_refs(k, sem, label):
    if not os.path.exists(REF_JSON):
        return
    ref = json.load(open(REF_JSON))
    kp, sp = ref['pyneut_mean'], ref['pyneut_sem']
    ko, so = ref['openmc_mean'], ref['openmc_sem']
    print(f"\n  {label} k_inf = {k:.5f} ± {sem:.5f}")
    print(f"  vs scalar PyNeut 5-seed ref {kp:.5f} ± {sp:.5f} : "
          f"dk = {1e5 * (k - kp):+,.0f} pcm, z = {abs(k-kp)/np.hypot(sem,sp):.2f}")
    print(f"  vs OpenMC 5-seed ref        {ko:.5f} ± {so:.5f} : "
          f"dk = {1e5 * (k - ko):+,.0f} pcm, z = {abs(k-ko)/np.hypot(sem,so):.2f}")


def main():
    ap = argparse.ArgumentParser(description='Vector-engine BEAVRS pin cell')
    ap.add_argument('--n', type=int, default=3000, help='neutrons/generation')
    ap.add_argument('--gens', type=int, default=140, help='total generations')
    ap.add_argument('--inactive', type=int, default=40, help='inactive generations')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--seeds', type=int, default=1,
                    help='number of independent replicas (>1 runs them in '
                         'parallel, one per core, and reports the seed-spread '
                         'error bar)')
    ap.add_argument('--workers', type=int, default=None,
                    help='parallel workers (default: one per replica)')
    ap.add_argument('--quick', action='store_true', help='fast smoke test')
    args = ap.parse_args()
    if args.quick:
        args.n, args.gens, args.inactive = 500, 30, 10

    matched = build_openmc_materials(temperature=294.0, thermal=False)
    comps = {k: pyneut_composition(m, 294.0) for k, m in matched.items()}
    mediums = pyneut_mediums(comps)
    settings = Settings('criticality', 1)

    print(f"\nVector-engine BEAVRS HZP pin cell (matched physics: free-gas H, "
          f"294 K, no S(a,b))\n  N={args.n}/gen, {args.gens} gens, "
          f"{args.inactive} inactive")

    if args.seeds > 1:
        # independent replications, one power iteration per core
        seeds = [args.seed + i for i in range(args.seeds)]
        sources = [pellet_source_arrays(args.n, np.random.default_rng(1000 + s))
                   for s in seeds]
        print(f"  {args.seeds} independent seeds in parallel "
              f"({args.workers or args.seeds} workers)\n")
        t0 = time.perf_counter()
        out = vt.run_keff_parallel(ENDF_PATH, mediums, sources, settings, seeds,
                                   generations=args.gens, inactive=args.inactive,
                                   n_workers=args.workers)
        dt = time.perf_counter() - t0
        n_hist = args.n * args.gens * args.seeds
        ks = out['per_seed_k']
        print("  per-seed k_inf : " + "  ".join(f"{k:.5f}" for k in ks))
        print(f"  wall time      : {dt:.1f} s for {n_hist:,} histories "
              f"across {args.seeds} replicas = {n_hist / dt:,.0f} hist/s")
        _compare_to_refs(out['k_mean'], out['k_spread_sem'],
                         f"{args.seeds}-seed mean (seed-spread SEM):")
        print("  (seed-spread SEM is the honest error bar for this DR~1 lattice)")
        return

    # distinct source seed: run_keff_vector builds its own default_rng(seed) for
    # transport, so sharing the seed here would correlate the first generation's
    # source directions with its flight distances (harmless for a reflective
    # lattice, but avoided for cleanliness).
    src = pellet_source_arrays(args.n, np.random.default_rng(1000 + args.seed))
    print(f"  single seed={args.seed}\n")
    t0 = time.perf_counter()
    out = vt.run_keff_vector(reader_single(), mediums, src, settings,
                             generations=args.gens, inactive=args.inactive,
                             seed=args.seed)
    dt = time.perf_counter() - t0
    n_hist = args.n * args.gens
    print(f"  k_inf (collision) : {out['k_eff']:.5f} ± {out['k_sem']:.5f}")
    print(f"  k_inf (source)    : {out['k_source']:.5f} ± {out['k_source_sem']:.5f}")
    print(f"  wall time         : {dt:.1f} s for {n_hist:,} histories "
          f"= {n_hist / dt:,.0f} hist/s (single core)")
    _compare_to_refs(out['k_eff'], out['k_sem'], "single-seed")
    print("  (single-seed run; the DR~1 lattice inflates single-run variance)")


def reader_single():
    return CrossSectionReader(ENDF_PATH, temperature='294K')


if __name__ == '__main__':
    main()
