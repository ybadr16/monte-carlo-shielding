"""
Cross-section interpolation check — PyNeut vs OpenMC

Compares PyNeut's interpolated microscopic cross sections against the values
OpenMC reads from the very same ENDF/B HDF5 file.  OpenMC data is loaded
directly with ``openmc.data.IncidentNeutron.from_hdf5`` so no cross_sections.xml
is required for this check.

Channels checked (where present in the file):
    MT 2   elastic scattering
    MT 102 radiative capture
    MT 16  (n,2n)

This isolates the reader/interpolation layer from the transport physics: a
mismatch here points at the data reader, not the Monte Carlo kernel.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
ENDF_PATH = os.path.join(PROJECT_ROOT, 'endfb')
sys.path.insert(0, PROJECT_ROOT)

try:
    import openmc.data
    HAS_OPENMC = True
except ImportError:
    HAS_OPENMC = False

from _common import get_reader

# Spot-check energies per isotope [eV]
_CHECKS = {
    'Pb208': [1e6, 2e6, 5e6, 14e6],
    'Fe56':  [1e5, 1e6, 5e6, 14e6],
    'Be9':   [1e6, 5e6, 14e6],
    'Al27':  [1e5, 1e6, 5e6],
    'C12':   [1e6, 2e6, 5e6],
    'O16':   [1e6, 5e6, 10e6],
}

# (MT, label)
_CHANNELS = [(2, 'elastic'), (102, 'capture'), (16, 'n,2n')]

_WARN_PCT = 1.0   # flag rows above this relative difference


def _openmc_xs(omc_data, mt, energy_ev):
    """Microscopic XS [barns] from an OpenMC IncidentNeutron object, or None."""
    if mt not in omc_data.reactions:
        return None
    rxn = omc_data[mt]
    temps = list(rxn.xs.keys())
    temp = '294K' if '294K' in temps else ('293.6K' if '293.6K' in temps else temps[0])
    return float(rxn.xs[temp](energy_ev))


def run():
    if not HAS_OPENMC:
        print('  [skip] openmc not installed — XS check unavailable')
        return []

    reader = get_reader()
    results = []
    all_ok = True

    print(f'\n  {"Isotope":<8} {"Energy":>10}  {"Channel":<8} '
          f'{"OpenMC":>13} {"PyNeut":>13} {"Diff":>9}')
    print(f'  {"-"*66}')

    for isotope, energies in _CHECKS.items():
        fpath = os.path.join(ENDF_PATH, 'neutron', f'{isotope}.h5')
        if not os.path.exists(fpath):
            print(f'  [skip] {isotope}: {isotope}.h5 not found')
            continue
        omc_data = openmc.data.IncidentNeutron.from_hdf5(fpath)

        for E in energies:
            for mt, label in _CHANNELS:
                omc = _openmc_xs(omc_data, mt, E)
                if omc is None:
                    continue   # channel absent for this isotope
                pyn = reader.get_cross_section(isotope, mt, E)   # barns

                if omc > 0:
                    diff = abs(omc - pyn) / omc * 100
                    diff_str = f'{diff:.3f}%'
                    ok = diff <= _WARN_PCT
                else:
                    diff = None
                    diff_str = 'N/A'
                    ok = True
                all_ok = all_ok and ok

                flag = '' if ok else ' !'
                print(f'  {isotope:<8} {E/1e6:>7.3f} MeV  {label:<8} '
                      f'{omc:>13.5e} {pyn:>13.5e} {diff_str:>9}{flag}')

                results.append({
                    'isotope': isotope, 'energy_eV': E, 'channel': label, 'mt': mt,
                    'openmc_xs': omc, 'pyneut_xs': pyn, 'diff_pct': diff,
                })

    print(f'\n  Overall XS interpolation check: '
          f'{"PASS" if all_ok else "FAIL (see ! rows)"}')
    return results


if __name__ == '__main__':
    print('\n=== Cross-Section Interpolation Check ===')
    run()
