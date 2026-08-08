"""
Master validation runner — PyNeut vs OpenMC

Runs all per-isotope validation scripts, collects results, and writes
  validation_results.csv
inside this directory.

Usage
-----
  cd Validation/OpenMC_Comparison
  python run_all.py [--n <particles>] [--skip-xs] [--isotopes Pb208,Fe56,...]

Default: N=10 000 particles per case, all isotopes, XS check included.
"""
import argparse
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import validate_Pb208
import validate_Fe56
import validate_Be9
import validate_Al27
import validate_C12
import validate_O16
import validate_thermal
import validate_xs

CSV_PATH = os.path.join(HERE, 'validation_results.csv')

_MODULES = {
    'Pb208':   validate_Pb208,
    'Fe56':    validate_Fe56,
    'Be9':     validate_Be9,
    'Al27':    validate_Al27,
    'C12':     validate_C12,
    'O16':     validate_O16,
    'Thermal': validate_thermal,
}

_CSV_FIELDS = [
    'case', 'isotope', 'geometry', 'energy_MeV',
    'openmc_leakage', 'openmc_leakage_sem',
    'pyneut_leakage', 'pyneut_leakage_sem',
    'leakage_diff_pct', 'leakage_zscore',
    'openmc_energy_eV', 'openmc_energy_sem',
    'pyneut_energy_eV', 'pyneut_energy_sem',
    'energy_diff_pct', 'energy_zscore',
    'spectrum_chi2_dof', 'spectrum_status',
    'nxn_events', 'parent_maxwellian_pct', 'child_maxwellian_pct', 'n3n_events',
    'openmc_source', 'status',
]


def _banner(text):
    bar = '=' * (len(text) + 4)
    print(f'\n{bar}')
    print(f'  {text}')
    print(bar)


def main():
    parser = argparse.ArgumentParser(description='PyNeut full validation suite')
    parser.add_argument('--n', type=int, default=10_000,
                        help='Particles per case (default 10 000)')
    parser.add_argument('--skip-xs', action='store_true',
                        help='Skip the cross-section interpolation check')
    parser.add_argument('--isotopes', type=str, default='',
                        help='Comma-separated subset, e.g. Pb208,Fe56')
    args = parser.parse_args()

    wanted = set(args.isotopes.split(',')) if args.isotopes else set(_MODULES)
    wanted = {k for k in wanted if k}   # drop empty strings

    all_results = []
    t0 = time.time()

    for name, mod in _MODULES.items():
        if name not in wanted:
            continue
        _banner(f'{name} Validation  (N={args.n:,})')
        results = mod.run(n_particles=args.n)
        all_results.extend(results)

    # XS check appends diagnostic rows with a different schema — write separately
    if not args.skip_xs:
        _banner('Cross-Section Interpolation Check')
        validate_xs.run()

    # -------------------------------------------------------------------------
    # Summary table
    # -------------------------------------------------------------------------
    _banner('SUMMARY')
    statuses = [r['status'] for r in all_results]
    counts = {s: statuses.count(s)
              for s in ('OK', 'WARNING', 'FAIL', 'NO_REF', 'NO_SIGMA')}
    print(f"\n  {'Case':<26} {'Integral':<9} "
          f"{'Leak z':>7} {'E z':>6}   {'Spectrum':<8} {'χ²/dof':>8} {'Max%':>6}")
    print(f"  {'-'*76}")
    for r in all_results:
        lz = f"{r['leakage_zscore']:.1f}" if r['leakage_zscore'] is not None else '-'
        ez = f"{r['energy_zscore']:.1f}"  if r['energy_zscore']  is not None else '-'
        x2 = f"{r['spectrum_chi2_dof']:.2f}" if r.get('spectrum_chi2_dof') is not None else '-'
        ss = r.get('spectrum_status') or '-'
        mx = f"{r['parent_maxwellian_pct']:.0f}" if r.get('parent_maxwellian_pct') is not None else '-'
        print(f"  {r['case']:<26} {r['status']:<9} "
              f"{lz:>7} {ez:>6}   {ss:<8} {x2:>8} {mx:>6}")
    print(f"\n  OK={counts['OK']}  WARNING={counts['WARNING']}  FAIL={counts['FAIL']}"
          f"  NO_REF={counts['NO_REF']}  NO_SIGMA={counts['NO_SIGMA']}")
    print("  (z = |OpenMC - PyNeut| / combined standard error; OK ≤ 2σ, WARNING ≤ 4σ)")
    print(f"  Total time: {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------------
    # Write CSV
    # -------------------------------------------------------------------------
    with open(CSV_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=_CSV_FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_results)
    print(f'\n  Results written to: {os.path.relpath(CSV_PATH, HERE)}')

    # Pin the results to the snapshot of the nuclear data they were measured
    # against. Written as a sidecar rather than a CSV comment so that readers
    # of the CSV (run_all_vector.load_reference) are unaffected. Two results
    # with different digests are not comparable at the pcm level -- see the
    # pin-cell data-snapshot episode in the paper.
    import json
    from datetime import date
    from _common import data_manifest
    man = dict(data_manifest(), measured=date.today().isoformat(),
               n_particles=args.n, applies_to=os.path.basename(CSV_PATH))
    man_path = os.path.join(HERE, 'validation_manifest.json')
    with open(man_path, 'w') as f:
        json.dump(man, f, indent=1)
    print(f"  data manifest: sha256={man['sha256']} ({man['n_files']} files)"
          f" -> {os.path.basename(man_path)}")


if __name__ == '__main__':
    main()
