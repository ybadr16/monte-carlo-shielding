"""
PyNeut vs OpenMC — C12 (graphite) validation cases

Cases
-----
  c12_cyl_2mev : Graphite cylinder R=10 cm, H=20 cm (half-height=10), E=2 MeV

Reference: hardcoded from Inelastic_Validations notebook.
"""
from _common import run_pyneut, make_result, print_result

# ---------------------------------------------------------------------------
# OpenMC reference (notebook, hardcoded)
# ---------------------------------------------------------------------------
_OMC = (1.0000, 1_003_104)   # leakage, avg_energy_eV

# ---------------------------------------------------------------------------
# Case definition
# ---------------------------------------------------------------------------
_CASE = {
    'name': 'c12_cyl_2mev',
    'el':   'C12',
    'rho':  2.267,
    'A':    12.011,
    'E':    2e6,
    'geo':  'cyl',
    'dims': [10.0, 10.0],   # [radius_cm, half_height_cm]
}


def run(n_particles=10_000):
    print(f"  Running {_CASE['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE, n_particles)
    omc_leak, omc_e = _OMC
    r = make_result(
        case_name    = _CASE['name'],
        isotope      = 'C12',
        geometry     = 'cyl',
        energy_mev   = 2.0,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= 'notebook',
    )
    print_result(r)
    return [r]


if __name__ == '__main__':
    print('\n=== C12 (Graphite) Validation ===')
    run()
