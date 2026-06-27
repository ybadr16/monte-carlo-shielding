"""
PyNeut vs OpenMC — Pb208 validation cases

Cases
-----
  pb208_elastic  : Pb208 sphere R=10 cm, E=2 MeV   (elastic-dominant)
  pb208_n2n      : Pb208 sphere R=10 cm, E=14 MeV  ((n,2n) active)

OpenMC reference values are computed live against the local ENDF/B data
library (endfb/cross_sections.xml, built automatically), using identical
geometry and source to PyNeut.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'pb208_elastic', 'el': 'Pb208', 'rho': 11.35, 'A': 207.2,
     'E': 2e6,  'geo': 'sphere', 'dims': [10.0]},
    {'name': 'pb208_5mev',    'el': 'Pb208', 'rho': 11.35, 'A': 207.2,
     'E': 5e6,  'geo': 'sphere', 'dims': [10.0]},   # inelastic onset
    {'name': 'pb208_n2n',     'el': 'Pb208', 'rho': 11.35, 'A': 207.2,
     'E': 14e6, 'geo': 'sphere', 'dims': [10.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='Pb208')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Pb208 Validation ===')
    run()
