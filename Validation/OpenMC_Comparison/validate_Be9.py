"""
PyNeut vs OpenMC — Be9 validation case

Cases
-----
  be9_n2n_14mev : Be9 sphere R=10 cm, E=14 MeV  ((n,2n) multiplier active)

OpenMC reference is computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'be9_n2n_14mev', 'el': 'Be9', 'rho': 1.848, 'A': 8.93478,
     'E': 14e6, 'geo': 'sphere', 'dims': [10.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='Be9')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Be9 Validation ===')
    run()
