"""
PyNeut vs OpenMC — C12 (graphite) validation case

Cases
-----
  c12_cyl_2mev : Graphite cylinder R=10 cm, half-height=10 cm, E=2 MeV

OpenMC reference is computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'c12_cyl_2mev', 'el': 'C12', 'rho': 2.267, 'A': 12.011,
     'E': 2e6, 'geo': 'cyl', 'dims': [10.0, 10.0]},   # [radius, half-height] cm
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='C12')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== C12 (Graphite) Validation ===')
    run()
