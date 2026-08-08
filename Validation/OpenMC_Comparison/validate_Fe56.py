"""
PyNeut vs OpenMC — Fe56 validation cases

Cases
-----
  fe56_slab_14mev   : Fe56 slab half-thickness=2 cm, E=14 MeV
  fe56_sphere_25kev : Fe56 sphere R=15 cm, E=25 keV (resonance region)

OpenMC references are computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'fe56_slab_14mev',   'el': 'Fe56', 'rho': 7.87, 'A': 55.4544,
     'E': 14e6, 'geo': 'slab',   'dims': [2.0]},
    {'name': 'fe56_sphere_1mev',  'el': 'Fe56', 'rho': 7.87, 'A': 55.4544,
     'E': 1e6,  'geo': 'sphere', 'dims': [10.0]},   # fast elastic + inelastic
    {'name': 'fe56_sphere_25kev', 'el': 'Fe56', 'rho': 7.87, 'A': 55.4544,
     'E': 25e3, 'geo': 'sphere', 'dims': [15.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='Fe56')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Fe56 Validation ===')
    run()
