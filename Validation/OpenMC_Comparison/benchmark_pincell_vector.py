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


def main():
    ap = argparse.ArgumentParser(description='Vector-engine BEAVRS pin cell')
    ap.add_argument('--n', type=int, default=3000, help='neutrons/generation')
    ap.add_argument('--gens', type=int, default=140, help='total generations')
    ap.add_argument('--inactive', type=int, default=40, help='inactive generations')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--quick', action='store_true', help='fast smoke test')
    args = ap.parse_args()
    if args.quick:
        args.n, args.gens, args.inactive = 500, 30, 10

    print(f"\nVector-engine BEAVRS HZP pin cell (matched physics: free-gas H, "
          f"294 K, no S(a,b))\n  N={args.n}/gen, {args.gens} gens, "
          f"{args.inactive} inactive, seed={args.seed}\n")

    matched = build_openmc_materials(temperature=294.0, thermal=False)
    comps = {k: pyneut_composition(m, 294.0) for k, m in matched.items()}
    mediums = pyneut_mediums(comps)
    reader = CrossSectionReader(ENDF_PATH, temperature='294K')
    settings = Settings('criticality', 1)

    src = pellet_source_arrays(args.n, np.random.default_rng(args.seed))

    t0 = time.perf_counter()
    out = vt.run_keff_vector(reader, mediums, src, settings,
                             generations=args.gens, inactive=args.inactive,
                             seed=args.seed)
    dt = time.perf_counter() - t0
    n_hist = args.n * args.gens

    print(f"  k_inf (collision) : {out['k_eff']:.5f} ± {out['k_sem']:.5f}")
    print(f"  k_inf (source)    : {out['k_source']:.5f} ± {out['k_source_sem']:.5f}")
    print(f"  wall time         : {dt:.1f} s for {n_hist:,} histories "
          f"= {n_hist / dt:,.0f} hist/s (single core)")

    if os.path.exists(REF_JSON):
        ref = json.load(open(REF_JSON))
        kp, sp = ref['pyneut_mean'], ref['pyneut_sem']
        ko, so = ref['openmc_mean'], ref['openmc_sem']
        zs = abs(out['k_eff'] - kp) / np.hypot(out['k_sem'], sp)
        zo = abs(out['k_eff'] - ko) / np.hypot(out['k_sem'], so)
        print(f"\n  vs scalar PyNeut 5-seed ref {kp:.5f} ± {sp:.5f} : "
              f"dk = {1e5 * (out['k_eff'] - kp):+,.0f} pcm, z = {zs:.2f}")
        print(f"  vs OpenMC 5-seed ref        {ko:.5f} ± {so:.5f} : "
              f"dk = {1e5 * (out['k_eff'] - ko):+,.0f} pcm, z = {zo:.2f}")
        print("  (single-seed run; the DR~1 lattice inflates single-run "
              "variance — seed-spread SEMs are the honest error bar)")


if __name__ == '__main__':
    main()
