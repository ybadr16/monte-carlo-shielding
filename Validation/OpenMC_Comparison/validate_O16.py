"""
PyNeut vs OpenMC — O16 validation cases

Cases
-----
  o16_sphere_10mev : O16 sphere R=20 cm, density=0.0014 g/cm3, E=10 MeV

Reference: pre-computed statepoint in benchmark_diverse/openmc_O16_sphere/.
Parameters derived from model.xml: 500 particles × 20 batches = N=10 000.
"""
import os
from _common import run_pyneut, make_result, print_result, read_statepoint, STATEPOINT_DIR

# ---------------------------------------------------------------------------
# Statepoint path
# ---------------------------------------------------------------------------
_SP = os.path.join(STATEPOINT_DIR, 'openmc_O16_sphere', 'statepoint.20.h5')

# ---------------------------------------------------------------------------
# Case definition
# ---------------------------------------------------------------------------
_CASE = {
    'name': 'o16_sphere_10mev',
    'el':   'O16',
    'rho':  0.0014,
    'A':    15.999,
    'E':    10e6,
    'geo':  'sphere',
    'dims': [20.0],
}


def run(n_particles=10_000):
    print(f"  Running {_CASE['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE, n_particles)
    omc_leak, omc_e = read_statepoint(_SP)
    src = 'statepoint' if omc_leak is not None else 'N/A'
    r = make_result(
        case_name    = _CASE['name'],
        isotope      = 'O16',
        geometry     = 'sphere',
        energy_mev   = 10.0,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= src,
    )
    print_result(r)
    return [r]


if __name__ == '__main__':
    print('\n=== O16 Validation ===')
    run()
