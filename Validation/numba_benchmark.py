"""Numba experiment benchmark: how much does JIT compilation help PyNeut?

Compares the pure-NumPy vector engine against Numba-compiled equivalents on a
C12 thermal sphere. Run: python Validation/numba_benchmark.py
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from src.cross_section_read import CrossSectionReader
from src.material import Material, Nuclide
from src.medium import Region, Sphere
from src.settings import Settings
from src.statistics import escape_statistics
from src.vt_calc import VelocitySampler
from src import vector_transport as vt
from src import numba_kernels as nk

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")


def _best(fn, reps=15):
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fn()
        _ = float(np.sum(out[0])) if isinstance(out, tuple) else 0.0
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    if not nk.HAVE_NUMBA:
        print("Numba not available.")
        return
    reader = CrossSectionReader(ENDF)
    el, rho, A, R, E = "C12", 1.7, 11.8969, 5.0, 0.0253
    mat = Material(el, rho, A, A)
    nuc = Nuclide(el, mat.number_density, A, atomic_mass=A)
    samp = nuc.sampler
    reader._build_fast_table(el)
    tbl = reader._fast_tables[el]
    scale = 1e-24 * nuc.number_density
    grid, el_mac, cap_mac = tbl["grid"], tbl["el"] * scale, tbl["cap"] * scale

    print("Compiling (cold)…")
    t0 = time.perf_counter()
    nk.warmup_all()
    nk.transport_sphere_thermal(np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2),
                                np.zeros(2), np.ones(2), np.full(2, E), R, grid,
                                el_mac, cap_mac, A, samp.beta, 10.0, 1e-4, 0.1, 1000)
    print(f"  cold compile: {time.perf_counter() - t0:.1f}s (one-time; cache=True amortizes)\n")

    n = 200000
    vn = np.full(n, 2200.0)
    rng = np.random.default_rng(0)
    mu = rng.uniform(-1, 1, n); phi = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - mu ** 2)
    u, v, w = s * np.cos(phi), s * np.sin(phi), mu
    Ea = np.full(n, E)

    print(f"{'kernel':<38}{'NumPy':>9}{'Numba':>9}{'speedup':>9}")
    print("-" * 65)

    t_np = _best(lambda: vt.sample_velocity_many(samp, vn, np.random.default_rng(1), return_mu=True))
    t_nb = _best(lambda: nk.svt_sample(vn, samp.beta))
    print(f"{'SVT free-gas sampler (rejection)':<38}{t_np*1e3:>7.1f}ms{t_nb*1e3:>7.1f}ms{t_np/t_nb:>8.1f}x")

    def np_chain():
        Ep, mul = vt.free_gas_elastic(Ea, A, samp, np.random.default_rng(2))
        return vt.rotate_directions(u, v, w, mul, np.random.default_rng(3))
    t_np = _best(np_chain)
    t_nb = _best(lambda: nk.free_gas_collision(Ea, u, v, w, A, samp.beta))
    print(f"{'fused free-gas collision':<38}{t_np*1e3:>7.1f}ms{t_nb*1e3:>7.1f}ms{t_np/t_nb:>8.1f}x")

    # full transport
    def source(seed):
        r = np.random.default_rng(seed); m = r.uniform(-1, 1, n); p = r.uniform(0, 2 * np.pi, n); ss = np.sqrt(1 - m ** 2)
        return {"x": np.zeros(n), "y": np.zeros(n), "z": np.zeros(n),
                "u": ss * np.cos(p), "v": ss * np.sin(p), "w": m, "energy": np.full(n, E)}
    region = Region([Sphere((0, 0, 0), R)], composition=[nuc])
    src = source(7)
    vt.run_transport(reader, [region], source(9), np.random.default_rng(9), Settings("shielding", n))
    t0 = time.perf_counter()
    rv = vt.run_transport(reader, [region], src, np.random.default_rng(7), Settings("shielding", n))
    t_np = time.perf_counter() - t0
    sv = escape_statistics(rv["leaked_weight"], rv["escape_energy"], n)
    s2 = source(7)
    t0 = time.perf_counter()
    lw, le = nk.transport_sphere_thermal(s2["x"].copy(), s2["y"].copy(), s2["z"].copy(),
                                         s2["u"].copy(), s2["v"].copy(), s2["w"].copy(),
                                         s2["energy"].copy(), R, grid, el_mac, cap_mac,
                                         A, samp.beta, 10.0, 1e-4, 0.1, 100000)
    t_nb = time.perf_counter() - t0
    sn = escape_statistics(lw, le, n)
    print(f"{'FULL transport (thermal sphere)':<38}{t_np*1e3:>7.0f}ms{t_nb*1e3:>7.0f}ms{t_np/t_nb:>8.1f}x")
    print(f"\n  full-transport physics: NumPy leak {sv['leakage']:.5f} <E> {sv['avg_energy']:.5f}"
          f" | Numba leak {sn['leakage']:.5f} <E> {sn['avg_energy']:.5f}")
    print(f"  throughput: NumPy {n/t_np/1e3:.0f}k hist/s -> Numba {n/t_nb/1e6:.1f}M hist/s (single process)")


if __name__ == "__main__":
    main()
