"""
Cross-section interpolation check — PyNeut vs h5py direct read

Compares PyNeut's interpolated sigma_t, sigma_el, sigma_in at spot-check
energies against values read straight from the ENDF HDF5 file via h5py.
No OpenMC dependency.

Isotopes checked: Pb208, Fe56, Be9, Al27, C12, O16
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
ENDF_PATH = os.path.join(PROJECT_ROOT, 'endfb')

sys.path.insert(0, PROJECT_ROOT)

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    print('[warn] h5py not installed — XS validation skipped')

from _common import get_reader

# ---------------------------------------------------------------------------
# Spot-check energies per isotope [eV]
# ---------------------------------------------------------------------------
_CHECKS = {
    'Pb208': [1e6, 2e6, 5e6, 14e6],
    'Fe56':  [1e5, 1e6, 5e6, 14e6],
    'Be9':   [1e6, 5e6, 14e6],
    'Al27':  [1e5, 1e6, 5e6],
    'C12':   [1e6, 2e6, 5e6],
    'O16':   [1e6, 5e6, 10e6],
}

# ENDF MT numbers for total, elastic, inelastic
_MT_TOTAL    = 1
_MT_ELASTIC  = 2
_MT_INELASTIC = 4   # sum of discrete + continuum


def _read_h5_xs(isotope, mt, energy_ev):
    """Interpolate a single reaction XS from the HDF5 file directly."""
    fpath = os.path.join(ENDF_PATH, 'neutron', f'{isotope}.h5')
    if not os.path.exists(fpath):
        return None
    with h5py.File(fpath, 'r') as f:
        root = f'neutron/{isotope}'
        rxn_key = f'reaction_{mt:03d}'
        try:
            grp = f[f'{root}/reactions/{rxn_key}']
            xs_grp = grp['xs/294K'] if '294K' in grp.get('xs', {}) else grp['xs']
            energies = xs_grp['energy'][:]
            values   = xs_grp['value'][:]
        except KeyError:
            return None
    return float(np.interp(energy_ev, energies, values))


def _pyneut_xs(reader, isotope, energy_ev):
    """Return (sigma_t, sigma_el, sigma_in) from PyNeut's reader."""
    try:
        sigma_t  = reader.get_total_cross_section(isotope, energy_ev)
        sigma_el = reader.get_elastic_cross_section(isotope, energy_ev)
        sigma_in = reader.get_inelastic_cross_section(isotope, energy_ev)
        return sigma_t, sigma_el, sigma_in
    except Exception as e:
        return None, None, None


def run():
    if not HAS_H5PY:
        print('  Skipping XS validation (h5py missing)')
        return []

    reader = get_reader()
    results = []
    warn_threshold = 1.0   # %

    print(f'\n  {"Isotope":<8} {"Energy":>10}  {"MT":<4} {"h5py":>12} {"PyNeut":>12} {"Diff":>8}')
    print(f'  {"-"*60}')

    all_ok = True
    for isotope, energies in _CHECKS.items():
        for E in energies:
            pyn_t, pyn_el, pyn_in = _pyneut_xs(reader, isotope, E)

            for mt_label, mt, pyn_val in [
                ('tot', _MT_TOTAL,     pyn_t),
                ('el',  _MT_ELASTIC,   pyn_el),
                ('in',  _MT_INELASTIC, pyn_in),
            ]:
                ref = _read_h5_xs(isotope, mt, E)
                if ref is None or pyn_val is None:
                    diff_str = 'N/A'
                    ok = True
                else:
                    diff = abs(ref - pyn_val) / ref * 100 if ref > 0 else 0.0
                    diff_str = f'{diff:.3f}%'
                    ok = diff <= warn_threshold
                    if not ok:
                        all_ok = False

                flag = '' if ok else ' !'
                e_str = f'{E/1e6:.3f} MeV'
                ref_str  = f'{ref:.5e}'  if ref     is not None else 'N/A'
                pyn_str  = f'{pyn_val:.5e}' if pyn_val is not None else 'N/A'
                print(f'  {isotope:<8} {e_str:>10}  {mt_label:<4} {ref_str:>12} {pyn_str:>12} {diff_str:>8}{flag}')

                results.append({
                    'isotope': isotope,
                    'energy_eV': E,
                    'mt': mt_label,
                    'h5py_xs': ref,
                    'pyneut_xs': pyn_val,
                    'diff_pct': None if ref is None or pyn_val is None else
                                abs(ref - pyn_val) / ref * 100 if ref > 0 else 0.0,
                })

    print(f'\n  Overall XS check: {"PASS" if all_ok else "FAIL (see ! rows)"}')
    return results


if __name__ == '__main__':
    print('\n=== Cross-Section Interpolation Check ===')
    run()
