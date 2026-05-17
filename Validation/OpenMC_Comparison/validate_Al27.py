"""
PyNeut vs OpenMC — Al27 validation cases

Cases
-----
  al27_slab_1mev : Al27 thin slab half-thickness=0.5 cm, E=1 MeV

Reference: hardcoded from Inelastic_Validations notebook.
"""
from _common import run_pyneut, make_result, print_result

# ---------------------------------------------------------------------------
# OpenMC reference (notebook, hardcoded)
# ---------------------------------------------------------------------------
_OMC = (0.9995, 933_522)   # leakage, avg_energy_eV

# ---------------------------------------------------------------------------
# Case definition
# ---------------------------------------------------------------------------
_CASE = {
    'name': 'al27_slab_1mev',
    'el':   'Al27',
    'rho':  2.70,
    'A':    26.982,
    'E':    1e6,
    'geo':  'slab',
    'dims': [0.5],
}


def run(n_particles=10_000):
    print(f"  Running {_CASE['name']} …")
    pyn_leak, pyn_e = run_pyneut(_CASE, n_particles)
    omc_leak, omc_e = _OMC
    r = make_result(
        case_name    = _CASE['name'],
        isotope      = 'Al27',
        geometry     = 'slab',
        energy_mev   = 1.0,
        omc_leak     = omc_leak,
        omc_energy   = omc_e,
        pyn_leak     = pyn_leak,
        pyn_energy   = pyn_e,
        openmc_source= 'notebook',
    )
    print_result(r)
    return [r]


if __name__ == '__main__':
    print('\n=== Al27 Validation ===')
    run()
