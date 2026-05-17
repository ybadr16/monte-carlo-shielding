"""
PyNeut vs OpenMC — Be9 validation cases

Cases
-----
  be9_n2n_14mev : Be9 sphere R=10 cm, E=14 MeV  ((n,2n) active)

Reference: hardcoded from Inelastic_Validations notebook.

Current status: FAIL — BUG 1 (energy too soft by ~4.4%).
Root cause: child neutron in (n,2n) falls back to Maxwellian evaporation
when rejection-sample fails 10 times.  The soft Maxwellian spectrum explains
the cold escape energy.
"""
from _common import run_pyneut, make_result, print_result

# ---------------------------------------------------------------------------
# OpenMC reference (notebook, hardcoded)
# ---------------------------------------------------------------------------
_OMC = (1.5848, 5_558_926)   # leakage, avg_energy_eV

# ---------------------------------------------------------------------------
# Case definition
# ---------------------------------------------------------------------------
_CASE = {
    'name': 'be9_n2n_14mev',
    'el':   'Be9',
    'rho':  1.848,
    'A':    9.012,
    'E':    14e6,
    'geo':  'sphere',
    'dims': [10.0],
}


def run(n_particles=10_000):
    print(f"  Running {_CASE['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE, n_particles)
    omc_leak, omc_e = _OMC
    r = make_result(
        case_name    = _CASE['name'],
        isotope      = 'Be9',
        geometry     = 'sphere',
        energy_mev   = 14.0,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= 'notebook',
    )
    print_result(r)
    return [r]


if __name__ == '__main__':
    print('\n=== Be9 Validation ===')
    run()
