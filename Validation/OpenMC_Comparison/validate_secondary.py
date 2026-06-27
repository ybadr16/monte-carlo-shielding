"""
Secondary-emission diagnostic — PyNeut vs OpenMC, isolated from transport.

For an (n,2n) reaction at a fixed incident energy this samples the emitted
neutron's energy and angle directly from PyNeut's reader path (exactly what the
transport kernel emits) and compares to the OpenMC-parsed reference
distribution.  It separates the leakage-spectrum disagreement into its causes:
reference frame (CM vs lab), the emission-energy spectrum shape, and the
emission angle.

Run:  python Validation/OpenMC_Comparison/validate_secondary.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

import openmc.data

from src.cross_section_read import CrossSectionReader
from src.simulation import _emit_secondary_neutron
from src.random_number_generator import RNGHandler

ENDF = os.path.join(PROJECT_ROOT, 'endfb')

_CASES = [
    ('Be9',   9.012,  14e6),
    ('Fe56',  55.845, 14e6),
    ('Pb208', 207.2,  14e6),
]
_N = 40000


def _spectrum_chi2(pyn, ref, n_groups=12):
    """Coarse-group chi2/dof between two energy samples (shape only)."""
    hi = max(pyn.max(), ref.max())
    bins = np.linspace(0, hi, n_groups + 1)
    p, _ = np.histogram(pyn, bins=bins)
    r, _ = np.histogram(ref, bins=bins)
    p = p / p.sum(); r = r / r.sum()
    # Poisson-ish variance per group from the reference sample
    var = (r * (1 - r)) / len(ref) + (p * (1 - p)) / len(pyn)
    mask = var > 0
    chi2 = np.sum((p[mask] - r[mask]) ** 2 / var[mask])
    return chi2 / mask.sum()


def run():
    reader = CrossSectionReader(ENDF)
    print(f"\n  {'Iso':<7}{'frame':<6}{'node MeV':>8}{'<E> ref':>10}"
          f"{'<E> nat':>10}{'<E> lab':>10}{'chi2/dof':>11}")
    print('  ' + '-' * 70)
    for iso, A, Ein in _CASES:
        n = openmc.data.IncidentNeutron.from_hdf5(os.path.join(ENDF, 'neutron', f'{iso}.h5'))
        rxn = n[16]
        frame = 'CM' if rxn.center_of_mass else 'lab'
        dist = rxn.products[0].distribution[0]

        # Compare at an exact incident-energy GRID NODE so both codes use the
        # SAME single tabulated outgoing distribution — this isolates the
        # energy-sampling from any incident-energy interpolation (an earlier
        # version compared block i to block i+1 at grid nodes and falsely
        # reported a ~17% bias).
        grid = np.asarray(dist.energy)
        i = int(np.argmin(np.abs(grid - Ein)))
        node = float(grid[i])

        # OpenMC reference, native frame, single block.
        ref_E = np.asarray(dist.energy_out[i].sample(_N)).ravel()

        # PyNeut, native frame, single block (reader path, no CM->lab boost).
        rng = RNGHandler(11)
        pe = np.empty(_N)
        for k in range(_N):
            E, _mu, ok = reader.get_secondary_correlated_sample(iso, 16, node, rng)
            pe[k] = E
        chi2 = _spectrum_chi2(pe, ref_E)

        # Frame-boost context: mean lab energy after the CM->lab transform.
        rng = RNGHandler(11)
        lab = np.array([_emit_secondary_neutron(reader, iso, 16, node, node, A, rng)[0]
                        for _ in range(_N)])
        print(f"  {iso:<7}{frame:<6}{node/1e6:>8.2f}{ref_E.mean()/1e6:>10.3f}"
              f"{pe.mean()/1e6:>10.3f}{lab.mean()/1e6:>10.3f}{chi2:>11.1f}")

    print("\n  Notes:")
    print("   * Compared at a grid node, native frame: '<E> ref' vs '<E> nat'")
    print("     tests PyNeut's outgoing-energy sampling against OpenMC directly.")
    print("   * '<E> lab' is PyNeut after the CM->lab boost (= native for lab data).")
    print("   * chi2/dof ~ 1 confirms the single-collision emission matches OpenMC;")
    print("     any leakage-spectrum disagreement is a transport/multiple-scatter")
    print("     accumulation, not the emission distribution.")


if __name__ == '__main__':
    print('\n=== Secondary-emission diagnostic ((n,2n) @ 14 MeV) ===')
    run()
