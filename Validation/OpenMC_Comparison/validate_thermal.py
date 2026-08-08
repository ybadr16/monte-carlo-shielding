"""
PyNeut vs OpenMC — thermal / epithermal regime

This is the regime PyNeut is *least* validated in, and these cases probe it
directly.  OpenMC is run with free-gas scattering (no S(α,β) tables loaded), so
the comparison is apples-to-apples against PyNeut's own free-gas model.

Cases
-----
  b10_epithermal_1kev : B10 sphere R=1 cm, E=1 keV   (strong 1/v absorption)
  cd112_thermal       : Cd112 sphere R=2 cm, E=0.0253 eV (thermal capture)
  c12_thermal         : C12 (graphite) sphere R=5 cm, E=0.0253 eV (pure scatterer)

History: these two thermal cases were once the suite's only failures. Below
THERMAL_KINEMATICS_CUTOFF (10 eV) the elastic branch now routes through the
free-gas kernel, so the target nucleus moves, neutrons can up-scatter and the
spectrum equilibrates to 294 K; three further errors in that kernel (flight-XS
re-broadening, the SVT branch probability, and the sampler mass) were closed
afterwards. Both cases now agree with OpenMC on the mean escape energy within
0.4 sigma with a GOOD spectrum shape. The remaining limitation is that thermal
scattering is free-gas only: there is no bound-atom S(alpha,beta) kernel.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'b10_epithermal_1kev', 'el': 'B10',   'rho': 2.34, 'A': 9.927,
     'E': 1.0e3,  'geo': 'sphere', 'dims': [1.0]},
    {'name': 'cd112_thermal',       'el': 'Cd112', 'rho': 8.65, 'A': 110.95,
     'E': 0.0253, 'geo': 'sphere', 'dims': [2.0]},
    {'name': 'c12_thermal',         'el': 'C12',   'rho': 1.70, 'A': 11.893,
     'E': 0.0253, 'geo': 'sphere', 'dims': [5.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles)
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== Thermal / Epithermal Validation ===')
    run()
