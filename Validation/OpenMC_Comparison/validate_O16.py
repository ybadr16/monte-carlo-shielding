"""
PyNeut vs OpenMC — O16 validation case

Cases
-----
  o16_sphere_10mev : O16 sphere R=20 cm, density=0.0014 g/cm3, E=10 MeV

OpenMC reference is computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'o16_sphere_10mev', 'el': 'O16', 'rho': 0.0014, 'A': 15.999,
     'E': 10e6, 'geo': 'sphere', 'dims': [20.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='O16')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== O16 Validation ===')
    run()
