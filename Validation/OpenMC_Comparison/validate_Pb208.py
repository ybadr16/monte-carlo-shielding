"""
PyNeut vs OpenMC — Pb208 validation cases

Cases
-----
  pb208_elastic  : Pb208 sphere R=10 cm, E=2 MeV   (elastic-dominant)
  pb208_n2n      : Pb208 sphere R=10 cm, E=14 MeV  ((n,2n) active)

OpenMC reference values come from the Inelastic_Validations notebook
(N=10 000, shielding mode, hardcoded here for reproducibility).
"""
from _common import run_pyneut, make_result, print_result

# ---------------------------------------------------------------------------
# OpenMC reference (notebook run, shielding mode, N=10000)
# ---------------------------------------------------------------------------
_OMC = {
    'pb208_elastic': (0.9998,  1_962_031),
    'pb208_n2n':     (1.5060,  5_275_196),
}

# ---------------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------------
_CASES = [
    {
        'name':  'pb208_elastic',
        'el':    'Pb208',
        'rho':   11.35,
        'A':     207.2,
        'E':     2e6,
        'geo':   'sphere',
        'dims':  [10.0],
    },
    {
        'name':  'pb208_n2n',
        'el':    'Pb208',
        'rho':   11.35,
        'A':     207.2,
        'E':     14e6,
        'geo':   'sphere',
        'dims':  [10.0],
    },
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        pyn_leak, pyn_e = run_pyneut(case, n_particles)
        omc_leak, omc_e = _OMC[case['name']]
        r = make_result(
            case_name    = case['name'],
            isotope      = 'Pb208',
            geometry     = case['geo'],
            energy_mev   = case['E'] / 1e6,
            omc_leak     = omc_leak,
            omc_energy   = omc_e,
            pyn_leak     = pyn_leak,
            pyn_energy   = pyn_e,
            openmc_source= 'notebook',
        )
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Pb208 Validation ===')
    run()
