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
import validate_xs

CSV_PATH = os.path.join(HERE, 'validation_results.csv')

_MODULES = {
    'Pb208': validate_Pb208,
    'Fe56':  validate_Fe56,
    'Be9':   validate_Be9,
    'Al27':  validate_Al27,
    'C12':   validate_C12,
    'O16':   validate_O16,
}

_CSV_FIELDS = [
    'case', 'isotope', 'geometry', 'energy_MeV',
    'openmc_leakage', 'pyneut_leakage', 'leakage_diff_pct',
    'openmc_energy_eV', 'pyneut_energy_eV', 'energy_diff_pct',
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
    counts = {s: statuses.count(s) for s in ('OK', 'WARNING', 'FAIL', 'NO_REF')}
    print(f"\n  {'Case':<28} {'Status':<10} {'Leakage diff':>14} {'Energy diff':>14}")
    print(f"  {'-'*70}")
    for r in all_results:
        ld = f"{r['leakage_diff_pct']:.2f}%" if r['leakage_diff_pct'] is not None else 'N/A'
        ed = f"{r['energy_diff_pct']:.2f}%"  if r['energy_diff_pct']  is not None else 'N/A'
        print(f"  {r['case']:<28} {r['status']:<10} {ld:>14} {ed:>14}")
    print(f"\n  OK={counts['OK']}  WARNING={counts['WARNING']}  FAIL={counts['FAIL']}  NO_REF={counts['NO_REF']}")
    print(f"  Total time: {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------------
    # Write CSV
    # -------------------------------------------------------------------------
    with open(CSV_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=_CSV_FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_results)
    print(f'\n  Results written to: {os.path.relpath(CSV_PATH, HERE)}')


if __name__ == '__main__':
    main()
