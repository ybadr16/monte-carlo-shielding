"""
PyNeut vs OpenMC — Fe56 validation cases

Cases
-----
  fe56_slab_14mev  : Fe56 slab half-thickness=2 cm, E=14 MeV
                     Reference: hardcoded from Inelastic_Validations notebook
  fe56_sphere_25kev: Fe56 sphere R=15 cm, E=25 keV (resonance region)
                     Reference: live statepoint in benchmark_diverse/

Both cases currently FAIL — BUG 2 (XS mismatch, likely missed inelastic MTs).
"""
import os
from _common import run_pyneut, make_result, print_result, read_statepoint, STATEPOINT_DIR

# ---------------------------------------------------------------------------
# OpenMC reference — slab case (notebook, hardcoded)
# ---------------------------------------------------------------------------
_OMC_SLAB = (1.1003, 4_732_820)   # leakage, avg_energy_eV

# ---------------------------------------------------------------------------
# Statepoint path — sphere 25 keV (pre-computed, stored on disk)
# ---------------------------------------------------------------------------
_SP_SPHERE = os.path.join(STATEPOINT_DIR, 'openmc_Fe56_sphere', 'statepoint.20.h5')

# ---------------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------------
_CASE_SLAB = {
    'name': 'fe56_slab_14mev',
    'el':   'Fe56',
    'rho':  7.87,
    'A':    55.845,
    'E':    14e6,
    'geo':  'slab',
    'dims': [2.0],
}

_CASE_SPHERE = {
    'name': 'fe56_sphere_25kev',
    'el':   'Fe56',
    'rho':  7.87,
    'A':    55.845,
    'E':    25e3,
    'geo':  'sphere',
    'dims': [15.0],
}


def run(n_particles=10_000):
    results = []

    # --- slab 14 MeV --------------------------------------------------------
    print(f"  Running {_CASE_SLAB['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE_SLAB, n_particles)
    omc_leak, omc_e = _OMC_SLAB
    r = make_result(
        case_name    = _CASE_SLAB['name'],
        isotope      = 'Fe56',
        geometry     = 'slab',
        energy_mev   = 14.0,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= 'notebook',
    )
    print_result(r)
    results.append(r)

    # --- sphere 25 keV (statepoint) -----------------------------------------
    print(f"  Running {_CASE_SPHERE['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE_SPHERE, n_particles)
    omc_leak, omc_e = read_statepoint(_SP_SPHERE)
    src = 'statepoint' if omc_leak is not None else 'N/A'
    r = make_result(
        case_name    = _CASE_SPHERE['name'],
        isotope      = 'Fe56',
        geometry     = 'sphere',
        energy_mev   = 0.025,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= src,
    )
    print_result(r)
    results.append(r)

    return results


if __name__ == '__main__':
    print('\n=== Fe56 Validation ===')
    run()
