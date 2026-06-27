"""
Analytic benchmarks for PyNeut — exact closed-form references (no OpenMC).

These validate individual physics pieces against textbook results, independent
of any other Monte Carlo code:

  1. Uncollided exponential attenuation
     A pencil beam through a slab: the fraction that crosses without any
     collision must follow Beer–Lambert,  T_uncollided = exp(-Σ_t · d),
     exactly, for every material (scattering cannot affect the uncollided
     component).  Exercises the data reader (Σ_t), the free-flight sampling
     and the geometry tracking through the full transport kernel.

  2. Elastic energy-loss kinematics
     For isotropic-CM elastic scattering off a nucleus of mass ratio A (no
     thermal motion), the post-collision energy ratio E'/E is uniform on
     [α, 1] with α = ((A-1)/(A+1))², so

         ⟨E'/E⟩ = (1 + α)/2 ,
         ⟨ln(E/E')⟩ = ξ = 1 + α·ln(α)/(1 - α)   (mean lethargy gain).

     Exercises ``physics.elastic_scattering`` directly.

Run:  python Validation/analytic_benchmarks.py
"""
import os
import sys

import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, PROJECT_ROOT)

from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Box
from src.physics import elastic_scattering
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import simulate_single_particle
from src.vt_calc import VelocitySampler

ENDF = os.path.join(PROJECT_ROOT, 'endfb')
_PASS_Z = 3.0     # |z| within 3σ counts as agreement


def _verdict(z):
    return 'PASS' if (z is not None and abs(z) <= _PASS_Z) else 'FAIL'


# ---------------------------------------------------------------------------
# Benchmark 1 — uncollided exponential attenuation
# ---------------------------------------------------------------------------

def benchmark_attenuation(element='Pb208', rho=11.35, amu=208.0, A=207.2,
                          E=2.0e6, n_particles=40000):
    print('\n' + '=' * 70)
    print('  Benchmark 1: uncollided attenuation  T = exp(-Σ_t·d)')
    print('=' * 70)
    reader = CrossSectionReader(ENDF)
    mat = Material(element, rho, amu, A)
    N = mat.number_density
    sampler = VelocitySampler(mat.kg_mass)

    # Total macroscopic cross section at the (deterministic, E>10 eV) energy.
    *_, sigma_t = reader.get_cross_sections(element, E, sampler, N, A)
    mfp = 1.0 / sigma_t
    print(f"  {element} @ {E/1e6:.2f} MeV:  Σ_t = {sigma_t:.5f} cm⁻¹   (mfp = {mfp:.3f} cm)")
    print(f"\n  {'d (mfp)':>8} {'d (cm)':>9} {'T_analytic':>12} {'T_PyNeut':>14} {'z':>7}  {'verdict':>7}")
    print('  ' + '-' * 66)

    settings = Settings('criticality', n_particles)   # analog, weights stay 1
    results = []
    for d_mfp in (0.5, 1.0, 2.0, 3.0):
        d = d_mfp * mfp
        slab = [Region(surfaces=[Box(-1e4, 1e4, -1e4, 1e4, 0.0, d)],
                       name='slab', priority=1, element=element)]
        rngs = [RNGHandler(20250 + i) for i in range(n_particles)]
        states = [{
            'x': 0.0, 'y': 0.0, 'z': 1e-9,
            'theta': 0.0, 'phi': 0.0,          # pencil beam along +z
            'energy': E, 'weight': 1.0, 'has_interacted': False,
        } for _ in rngs]
        args = [(s, reader, slab, A, N, VelocitySampler(mat.kg_mass), None, False, r, settings)
                for s, r in zip(states, rngs)]
        with Pool() as pool:
            res = pool.map(simulate_single_particle, args)

        k = sum(1 for r in res if r['source_uncollided'])
        T = k / n_particles
        T_an = np.exp(-sigma_t * d)
        sem = np.sqrt(max(T * (1 - T), 1e-12) / n_particles)
        z = (T - T_an) / sem if sem > 0 else None
        v = _verdict(z)
        results.append({'benchmark': 'attenuation', 'd_mfp': d_mfp,
                        'analytic': T_an, 'pyneut': T, 'z': z, 'verdict': v})
        print(f"  {d_mfp:>8.1f} {d:>9.3f} {T_an:>12.5f} "
              f"{T:>9.5f}±{sem:.5f} {z:>7.1f}  {v:>7}")
    return results


# ---------------------------------------------------------------------------
# Benchmark 2 — elastic energy-loss kinematics
# ---------------------------------------------------------------------------

def _sample_elastic_ratios(A, E, n):
    sampler = VelocitySampler(1.0)            # mass irrelevant at E > 10 eV
    rng = RNGHandler(777)
    out = np.empty(n)
    for i in range(n):
        Ep, _, _ = elastic_scattering(E, A, sampler, rng)
        out[i] = Ep / E
    return out


def benchmark_elastic_kinematics(n=200000):
    print('\n' + '=' * 70)
    print('  Benchmark 2: elastic kinematics  (isotropic-CM, no thermal motion)')
    print('=' * 70)
    print(f"\n  {'A':>6} {'quantity':<12} {'analytic':>12} {'PyNeut':>16} {'z':>7}  {'verdict':>7}")
    print('  ' + '-' * 66)

    results = []
    E = 1.0e6
    for A in (12.0, 56.0, 207.0):
        alpha = ((A - 1) / (A + 1)) ** 2
        ratios = _sample_elastic_ratios(A, E, n)

        # (a) bounds: every sample within [alpha, 1]
        frac_in = np.mean((ratios >= alpha - 1e-9) & (ratios <= 1.0 + 1e-9))

        # (b) mean energy ratio
        mean_an = (1 + alpha) / 2
        mean_pyn = ratios.mean()
        sem_mean = ratios.std(ddof=1) / np.sqrt(n)
        z_mean = (mean_pyn - mean_an) / sem_mean

        # (c) mean lethargy gain xi
        leth = -np.log(ratios)
        xi_an = 1 + alpha * np.log(alpha) / (1 - alpha)
        xi_pyn = leth.mean()
        sem_xi = leth.std(ddof=1) / np.sqrt(n)
        z_xi = (xi_pyn - xi_an) / sem_xi

        vbound = 'PASS' if frac_in > 0.999 else 'FAIL'
        results.append({'benchmark': 'elastic_bounds', 'A': A,
                        'analytic': 1.0, 'pyneut': frac_in, 'z': None, 'verdict': vbound})
        results.append({'benchmark': 'elastic_mean_ratio', 'A': A,
                        'analytic': mean_an, 'pyneut': mean_pyn, 'z': z_mean,
                        'verdict': _verdict(z_mean)})
        results.append({'benchmark': 'elastic_xi', 'A': A,
                        'analytic': xi_an, 'pyneut': xi_pyn, 'z': z_xi,
                        'verdict': _verdict(z_xi)})

        print(f"  {A:>6.0f} {'in [α,1]':<12} {'100%':>12} {frac_in*100:>15.3f}% {'':>7}  {vbound:>7}")
        print(f"  {'':>6} {'⟨E'+chr(39)+'/E⟩':<12} {mean_an:>12.5f} {mean_pyn:>16.5f} {z_mean:>7.1f}  {_verdict(z_mean):>7}")
        print(f"  {'':>6} {'ξ (lethargy)':<12} {xi_an:>12.5f} {xi_pyn:>16.5f} {z_xi:>7.1f}  {_verdict(z_xi):>7}")
    return results


def main():
    all_res = []
    all_res += benchmark_attenuation()
    all_res += benchmark_elastic_kinematics()

    n_fail = sum(1 for r in all_res if r['verdict'] == 'FAIL')
    print('\n' + '=' * 70)
    print(f"  ANALYTIC BENCHMARKS: {len(all_res)-n_fail}/{len(all_res)} PASS, {n_fail} FAIL")
    print('  (PASS = within 3σ of the exact analytic result)')
    print('=' * 70)
    return all_res


if __name__ == '__main__':
    main()
