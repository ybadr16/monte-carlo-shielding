"""k-eigenvalue (criticality) driver.

Runs the fission-source power iteration on top of the fixed-source transport
kernel. Each generation transports its source bank with `simulate_single_particle`
(analog `criticality` mode), which banks ν prompt neutrons per fission into a
per-history `fission_sites` list. The driver pools those into the next
generation's source and estimates k from the generation-to-generation
multiplication

    k_gen = (fission neutrons produced) / (source neutrons started),

discarding `inactive` settling generations before averaging k over the active
ones with batch-means statistics (`statistics.mean_sem`). The fission bank is
resampled back to a fixed nominal size each cycle so the population neither dies
nor explodes; that resampling controls only the population, not the k estimate
(which is scored before resampling).

A second, lower-variance **collision estimator** is accumulated in parallel:
every collision contributes w·Σᵢ νᵢ·Σ_f,ᵢ / Σ_t (scored in the transport kernel
as `physics_counts["fission_production"]`), so k_collision = (total fission
production) / (source neutrons started). Both estimators are unbiased and agree
in the mean; the collision one is reported as the primary k_eff.
"""
import numpy as np

from .random_number_generator import RNGHandler
from .simulation import simulate_single_particle
from .statistics import mean_sem
from .physics import sample_watt_spectrum, watt_params_for


def uniform_sphere_source(center, radius, n, rng, element="U235", energy=None):
    """`n` isotropic source neutrons spread uniformly through a sphere.

    Energies are sampled from the fissile isotope's Watt spectrum unless a fixed
    `energy` (eV) is given.
    """
    cx, cy, cz = center
    a_w, b_w = watt_params_for(element)
    states = []
    r2 = radius * radius
    for _ in range(n):
        while True:
            x = rng.uniform(-radius, radius)
            y = rng.uniform(-radius, radius)
            z = rng.uniform(-radius, radius)
            if x * x + y * y + z * z <= r2:
                break
        states.append({
            "x": cx + x, "y": cy + y, "z": cz + z,
            "theta": np.arccos(2 * rng.random() - 1),
            "phi": 2 * np.pi * rng.random(),
            "energy": energy if energy is not None else sample_watt_spectrum(rng, a_w, b_w),
            "weight": 1.0, "has_interacted": False,
        })
    return states


def _resample(sites, target, rng):
    """Resample the fission bank to exactly `target` source neutrons."""
    n = len(sites)
    if n == 0:
        return []
    if n == target:
        return [dict(s) for s in sites]
    idx = rng.choice(n, size=target, replace=(n < target))
    return [dict(sites[i]) for i in idx]


def run_keff(reader, mediums, initial_source, A, N, sampler, settings,
             generations=50, inactive=10, seed=1, pool=None, verbose=False):
    """Estimate k_eff by fission-source iteration.

    Parameters mirror the fixed-source call (reader, mediums, A, N, sampler,
    settings) and add the generation schedule. `settings.mode` must be
    "criticality" (analog transport) so fission neutrons are produced. Pass a
    `multiprocessing.Pool` to parallelise the per-generation transport.

    Returns a dict: k_eff/k_sem (collision estimator), k_source/k_source_sem
    (source estimator), active_k (per-active-generation collision k), and
    cycles (list of (generation, k_collision, k_source, is_active)).
    """
    gen_size = len(initial_source)
    if gen_size == 0:
        raise ValueError("initial_source must contain at least one neutron")

    source = initial_source
    cycles = []
    active_src = []   # per-active-generation source (analog) estimator
    active_col = []   # per-active-generation collision estimator
    rng_state = np.random.default_rng(seed)
    seed_counter = int(seed) * 1_000_003

    for gen in range(generations):
        rngs = [RNGHandler(seed_counter + i) for i in range(len(source))]
        seed_counter += len(source) + 1
        args = [(s, reader, mediums, A, N, sampler, None, False, r, settings)
                for s, r in zip(source, rngs)]

        if pool is not None:
            results = pool.map(simulate_single_particle, args)
        else:
            results = [simulate_single_particle(a) for a in args]

        new_sites = []
        production = 0.0
        for res in results:
            new_sites.extend(res["fission_sites"])
            production += res["physics_counts"].get("fission_production", 0.0)

        started = len(source)
        produced = len(new_sites)
        k_src = produced / started
        k_col = production / started
        is_active = gen >= inactive
        cycles.append((gen, k_col, k_src, is_active))
        if is_active:
            active_src.append(k_src)
            active_col.append(k_col)

        if verbose:
            tag = "active" if is_active else "settle"
            print(f"  gen {gen:3d}  k_col={k_col:.5f}  k_src={k_src:.5f}  "
                  f"src={started:5d}  fiss={produced:5d}  [{tag}]")

        if produced == 0:
            break  # chain died out (deeply subcritical)

        source = _resample(new_sites, gen_size, rng_state)

    # the collision estimator is lower-variance, so it is the primary k_eff
    k_eff, k_sem = mean_sem(active_col)
    k_src_eff, k_src_sem = mean_sem(active_src)
    return {
        "k_eff": k_eff,
        "k_sem": k_sem,
        "k_source": k_src_eff,
        "k_source_sem": k_src_sem,
        "active_k": active_col,
        "cycles": cycles,
    }
