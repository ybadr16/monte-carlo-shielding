"""Vector-engine validation — all 13 fixed-source cases + speedup table.

Runs every case from the per-isotope validate_*.py modules through the
event-based vector engine (src/vector_transport.py) and grades it by z-score
against the committed OpenMC reference values in validation_results.csv (the
same references the paper uses; OpenMC is NOT re-run). The scalar kernel is
timed on a smaller sample of the same cases to produce a per-case speedup
column. Both engines share one physics, so the z-grading is the same
statistical comparison the scalar suite uses — only the RNG consumption
order differs.

Usage
-----
  cd Validation/OpenMC_Comparison
  python run_all_vector.py [--n 100000] [--scalar-n 600] [--no-scalar-timing]

Writes validation_results_vector.csv next to this script.
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from _common import ENDF_PATH, PROJECT_ROOT, energy_grid, N_BATCHES  # noqa: E402

sys.path.insert(0, PROJECT_ROOT)

from src.cross_section_read import CrossSectionReader          # noqa: E402
from src.material import Material, Nuclide, NEUTRON_MASS_AMU   # noqa: E402
from src.medium import Region, Sphere, Cylinder, Box, Plane    # noqa: E402
from src.random_number_generator import RNGHandler             # noqa: E402
from src.settings import Settings                              # noqa: E402
from src.simulation import simulate_single_particle            # noqa: E402
from src.statistics import escape_statistics                   # noqa: E402
from src.vt_calc import VelocitySampler                        # noqa: E402
from src import vector_transport as vt                         # noqa: E402

import validate_Pb208    # noqa: E402
import validate_Fe56     # noqa: E402
import validate_Be9      # noqa: E402
import validate_Al27     # noqa: E402
import validate_C12      # noqa: E402
import validate_O16      # noqa: E402
import validate_thermal  # noqa: E402

CASES = (validate_Pb208._CASES + validate_Fe56._CASES + validate_Be9._CASES
         + validate_Al27._CASES + validate_C12._CASES + validate_O16._CASES
         + validate_thermal._CASES)

REF_CSV = os.path.join(HERE, 'validation_results.csv')
OUT_CSV = os.path.join(HERE, 'validation_results_vector.csv')


def load_reference():
    """OpenMC + scalar-PyNeut reference rows from the committed CSV."""
    ref = {}
    with open(REF_CSV, newline='') as f:
        for row in csv.DictReader(f):
            ref[row['case']] = {k: (float(v) if v not in ('', None) else None)
                                for k, v in row.items()
                                if k not in ('case', 'isotope', 'geometry',
                                             'spectrum_status', 'openmc_source',
                                             'status')}
    return ref


def build_geometry(case, nuc=None):
    """Same geometry as _common.run_pyneut; composition-based when `nuc` is
    given (vector engine), element-based otherwise (scalar kernel)."""
    kw = ({'composition': [nuc]} if nuc is not None
          else {'element': case['el']})
    geo, dims = case['geo'], case['dims']
    if geo == 'sphere':
        surfaces = [Sphere((0, 0, 0), dims[0])]
    elif geo == 'slab':
        d = dims[0]
        surfaces = [Box(-d, d, -50, 50, -50, 50)]
    elif geo == 'cyl':
        R, H = dims
        surfaces = [Cylinder('z', R, (0, 0, 0)),
                    Plane(0, 0, 1, H), Plane(0, 0, -1, H)]
    else:
        raise ValueError(f"Unknown geometry: {geo!r}")
    return [Region(surfaces=surfaces, operation='intersection',
                   name=case['name'], priority=1, **kw)]


def run_vector_case(reader, case, n):
    """Vector-engine run: returns (stats dict, histories/s)."""
    # case['A'] is the atomic-weight ratio; the molar mass is A * m_neutron
    # (see _common.run_pyneut). Passing A as the mass makes the free-gas target
    # ~0.9% too light and biases a light moderator's thermal spectrum hot.
    amass = case['A'] * NEUTRON_MASS_AMU
    mat = Material(case['el'], case['rho'], amass, case['A'])
    nuc = Nuclide(case['el'], mat.number_density, case['A'], atomic_mass=amass)
    mediums = build_geometry(case, nuc=nuc)
    settings = Settings('shielding', n)

    def source(m, seed):
        rng = np.random.default_rng(seed)
        mu = rng.uniform(-1, 1, m)
        phi = rng.uniform(0, 2 * np.pi, m)
        s = np.sqrt(1 - mu ** 2)
        return rng, {'x': np.zeros(m), 'y': np.zeros(m), 'z': np.zeros(m),
                     'u': s * np.cos(phi), 'v': s * np.sin(phi), 'w': mu,
                     'energy': np.full(m, case['E'])}

    # warm the fast tables / angular / secondary caches outside the timing
    rng_w, src_w = source(64, 1)
    vt.run_transport(reader, mediums, src_w, rng_w, settings)

    rng, src = source(n, 12345)
    t0 = time.perf_counter()
    res = vt.run_transport(reader, mediums, src, rng, settings)
    dt = time.perf_counter() - t0

    stats = escape_statistics(res['leaked_weight'], res['escape_energy'],
                              n, N_BATCHES, e_bins=energy_grid(case['E']))
    stats['counters'] = res['counters']
    return stats, n / dt


def scalar_rate(reader, case, n):
    """Single-core scalar-kernel throughput on the same case (histories/s)."""
    # Same mass convention as run_vector/_common: A is the atomic-weight ratio,
    # so the molar mass is A * m_neutron (not A).
    mat = Material(case['el'], case['rho'], case['A'] * NEUTRON_MASS_AMU,
                   case['A'])
    mediums = build_geometry(case)
    settings = Settings('shielding', n)
    sampler = VelocitySampler(mat.kg_mass)
    rngs = [RNGHandler(12345 + i) for i in range(n)]
    states = [{'x': 0.0, 'y': 0.0, 'z': 0.0,
               'theta': np.arccos(1 - 2 * r.random()),
               'phi': r.uniform(0, 2 * np.pi),
               'has_interacted': False, 'energy': case['E'], 'weight': 1.0}
              for r in rngs]
    t0 = time.perf_counter()
    for s, r in zip(states, rngs):
        simulate_single_particle((s, reader, mediums, case['A'],
                                  mat.number_density, sampler, None, False,
                                  r, settings))
    return n / (time.perf_counter() - t0)


def zscore(a, sa, b, sb):
    sa = sa if (sa is not None and np.isfinite(sa)) else 0.0
    sb = sb if (sb is not None and np.isfinite(sb)) else 0.0
    comb = np.hypot(sa, sb)
    return abs(a - b) / comb if comb > 0 else None


def main():
    ap = argparse.ArgumentParser(description='Vector-engine validation suite')
    ap.add_argument('--n', type=int, default=100_000,
                    help='vector histories per case (default 100 000)')
    ap.add_argument('--scalar-n', type=int, default=600,
                    help='scalar histories per case for timing (default 600)')
    ap.add_argument('--no-scalar-timing', action='store_true',
                    help='skip the scalar throughput measurement')
    args = ap.parse_args()

    ref = load_reference()
    reader = CrossSectionReader(ENDF_PATH)
    rows = []
    t0 = time.time()

    print(f"\n  Vector-engine suite: N={args.n:,} per case, "
          f"references = committed validation_results.csv (OpenMC not re-run)")
    print(f"\n  {'Case':<24} {'Leak z':>7} {'E z':>6} {'Status':<8} "
          f"{'vec hist/s':>11} {'scal hist/s':>12} {'speedup':>8}")
    print(f"  {'-' * 82}")

    for case in CASES:
        r = ref.get(case['name'])
        stats, vrate = run_vector_case(reader, case, args.n)
        srate = (None if args.no_scalar_timing
                 else scalar_rate(reader, case, args.scalar_n))

        z_leak = z_energy = None
        if r and r.get('openmc_leakage') is not None:
            z_leak = zscore(stats['leakage'], stats['leakage_sem'],
                            r['openmc_leakage'], r['openmc_leakage_sem'])
            z_energy = zscore(stats['avg_energy'], stats['avg_energy_sem'],
                              r['openmc_energy_eV'], r['openmc_energy_sem'])
        zs = [zv for zv in (z_leak, z_energy) if zv is not None]
        zmax = max(zs) if zs else None
        status = ('NO_REF' if zmax is None else
                  'OK' if zmax <= 2 else 'WARNING' if zmax <= 4 else 'FAIL')

        lz = f"{z_leak:.1f}" if z_leak is not None else '-'
        ez = f"{z_energy:.1f}" if z_energy is not None else '-'
        sr = f"{srate:,.0f}" if srate else '-'
        sp = f"{vrate / srate:.0f}x" if srate else '-'
        print(f"  {case['name']:<24} {lz:>7} {ez:>6} {status:<8} "
              f"{vrate:>11,.0f} {sr:>12} {sp:>8}")

        rows.append({
            'case': case['name'], 'isotope': case['el'],
            'geometry': case['geo'], 'energy_MeV': case['E'] / 1e6,
            'openmc_leakage': r.get('openmc_leakage') if r else None,
            'openmc_leakage_sem': r.get('openmc_leakage_sem') if r else None,
            'vector_leakage': round(stats['leakage'], 5),
            'vector_leakage_sem': round(stats['leakage_sem'], 5),
            'leakage_zscore': round(z_leak, 2) if z_leak is not None else None,
            'openmc_energy_eV': r.get('openmc_energy_eV') if r else None,
            'openmc_energy_sem': r.get('openmc_energy_sem') if r else None,
            'vector_energy_eV': round(stats['avg_energy'], 4),
            'vector_energy_sem': round(stats['avg_energy_sem'], 4),
            'energy_zscore': round(z_energy, 2) if z_energy is not None else None,
            'vector_hist_per_s': round(vrate),
            'scalar_hist_per_s': round(srate) if srate else None,
            'speedup': round(vrate / srate, 1) if srate else None,
            'nxn_events': stats['counters']['nxn_events'],
            'status': status,
        })

    statuses = [r['status'] for r in rows]
    print(f"\n  OK={statuses.count('OK')}  WARNING={statuses.count('WARNING')}"
          f"  FAIL={statuses.count('FAIL')}  (z vs committed OpenMC references)")
    print(f"  Total time: {time.time() - t0:.1f}s")

    with open(OUT_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Results written to: {os.path.relpath(OUT_CSV, HERE)}")


if __name__ == '__main__':
    main()
