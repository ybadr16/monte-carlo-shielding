"""
PyNeut vs OpenMC — Al27 validation case

Cases
-----
  al27_slab_1mev : Al27 thin slab half-thickness=0.5 cm, E=1 MeV

OpenMC reference is computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'al27_slab_1mev', 'el': 'Al27', 'rho': 2.70, 'A': 26.982,
     'E': 1e6, 'geo': 'slab', 'dims': [0.5]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='Al27')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Al27 Validation ===')
    run()
