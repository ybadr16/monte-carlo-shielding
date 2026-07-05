"""Experimental Numba-compiled hot-path kernels for the vector engine.

This is an exploration branch: it re-implements the vector engine's hottest
inner functions with ``@njit`` and measures whether a JIT-compiled loop beats
the vectorized-NumPy version. Numba wins where NumPy is weak (per-element
rejection loops with early exit) and roughly ties where NumPy is already a thin
wrapper over C (dense elementwise array math). The functions here mirror the
physics of their vector_transport counterparts but sample from Numba's own
per-thread RNG, so results match statistically, not stream-for-stream.

Import is guarded: if Numba is unavailable, HAVE_NUMBA is False and callers
fall back to the pure-NumPy path.
"""
import math

import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        def deco(f):
            return f
        return deco if not args else args[0]

    prange = range


_SQRT_PI = math.sqrt(math.pi)


@njit(parallel=True, fastmath=True, cache=True)
def svt_sample(vn, beta):
    """Free-gas (SVT) target speeds and correlated cosines, one lane at a time.

    Per-lane acceptance-rejection with early exit — the shape Numba compiles
    well and masked-array NumPy handles poorly (every rejected lane still costs
    a full-width array op each iteration). ``prange`` also spreads the lanes
    across threads. Reproduces the physics of
    vector_transport.sample_velocity_many including the SVT branch probability
    in the dimensionless neutron speed ``beta*vn``.
    """
    n = vn.size
    v_t = np.empty(n)
    mu = np.empty(n)
    for i in prange(n):
        vni = vn[i]
        p_first = 2.0 / (_SQRT_PI * vni * beta + 2.0)
        while True:
            r0 = np.random.random()
            r1 = np.random.random()
            r2 = np.random.random()
            r3 = np.random.random()
            r4 = np.random.random()
            r5 = np.random.random()
            if r0 < p_first:
                x = math.sqrt(-math.log((1.0 - r1) * (1.0 - r2)))
            else:
                c = math.cos(0.5 * math.pi * r3)
                x = math.sqrt(-math.log(1.0 - r1) - math.log(1.0 - r2) * c * c)
            vt_i = x / beta
            mu_i = 2.0 * r4 - 1.0
            acc = math.sqrt(vni * vni + vt_i * vt_i
                            - 2.0 * vni * vt_i * mu_i) / (vni + vt_i)
            if r5 < acc:
                v_t[i] = vt_i
                mu[i] = mu_i
                break
    return v_t, mu


def warmup():
    """Trigger compilation so the first real call isn't paying for it."""
    if HAVE_NUMBA:
        svt_sample(np.full(4, 2200.0), 1.5e-3)


_M_N = 1.674927471e-27
_EV = 1.60217663e-19


@njit(parallel=True, fastmath=True, cache=True)
def free_gas_collision(E, u, v, w, A, beta):
    """Whole free-gas elastic collision per lane, fused: sample the target
    (SVT), boost to CM, isotropic-CM exit, boost back, and rotate the incident
    direction by the resulting lab cosine. One compiled loop replaces
    sample_velocity_many + _isotropic_cm_exit + rotate_directions and their
    intermediate arrays. Returns (E', u', v', w')."""
    n = E.size
    Eo = np.empty(n)
    uo = np.empty(n)
    vo = np.empty(n)
    wo = np.empty(n)
    Ap1 = A + 1.0
    for i in prange(n):
        v_n = math.sqrt(2.0 * E[i] * _EV / _M_N)
        # --- SVT target ---
        p_first = 2.0 / (_SQRT_PI * v_n * beta + 2.0)
        while True:
            r0 = np.random.random(); r1 = np.random.random(); r2 = np.random.random()
            r3 = np.random.random(); r4 = np.random.random(); r5 = np.random.random()
            if r0 < p_first:
                x = math.sqrt(-math.log((1.0 - r1) * (1.0 - r2)))
            else:
                c = math.cos(0.5 * math.pi * r3)
                x = math.sqrt(-math.log(1.0 - r1) - math.log(1.0 - r2) * c * c)
            v_t = x / beta
            mu_t = 2.0 * r4 - 1.0
            acc = math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t)
            if r5 < acc:
                break
        # --- target vector, CM boost ---
        phi_t = 2.0 * math.pi * np.random.random()
        sin_t = math.sqrt(max(0.0, 1.0 - mu_t * mu_t)) if v_t > 0.0 else 0.0
        vtx = v_t * sin_t * math.cos(phi_t) if v_t > 0.0 else 0.0
        vty = v_t * sin_t * math.sin(phi_t) if v_t > 0.0 else 0.0
        vtz = v_t * mu_t if v_t > 0.0 else 0.0
        vcm_x = A * vtx / Ap1
        vcm_y = A * vty / Ap1
        vcm_z = (v_n + A * vtz) / Ap1
        Vx = -vcm_x; Vy = -vcm_y; Vz = v_n - vcm_z
        Vmag = math.sqrt(Vx * Vx + Vy * Vy + Vz * Vz)
        # --- isotropic CM exit, boost back ---
        mu_cm = 2.0 * np.random.random() - 1.0
        phi_cm = 2.0 * math.pi * np.random.random()
        sin_cm = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
        vpx = Vmag * sin_cm * math.cos(phi_cm) + vcm_x
        vpy = Vmag * sin_cm * math.sin(phi_cm) + vcm_y
        vpz = Vmag * mu_cm + vcm_z
        sp2 = vpx * vpx + vpy * vpy + vpz * vpz
        Eo[i] = max(1e-5, 0.5 * _M_N * sp2 / _EV)
        sp = math.sqrt(sp2)
        mu_lab = vpz / sp if sp > 0.0 else 1.0
        if mu_lab > 1.0:
            mu_lab = 1.0
        elif mu_lab < -1.0:
            mu_lab = -1.0
        # --- rotate incident (u,v,w) by mu_lab ---
        ui = u[i]; vi = v[i]; wi = w[i]
        phi = 2.0 * math.pi * np.random.random()
        cphi = math.cos(phi); sphi = math.sin(phi)
        sin_theta = math.sqrt(max(0.0, 1.0 - mu_lab * mu_lab))
        if abs(wi) >= 0.999999:
            sign = 1.0 if wi > 0.0 else -1.0
            uo[i] = sin_theta * cphi
            vo[i] = sin_theta * sphi
            wo[i] = sign * mu_lab
        else:
            denom = math.sqrt(max(1e-12, 1.0 - wi * wi))
            ratio = sin_theta / denom
            un = mu_lab * ui + ratio * (ui * wi * cphi - vi * sphi)
            vn2 = mu_lab * vi + ratio * (vi * wi * cphi + ui * sphi)
            wn = mu_lab * wi - sin_theta * denom * cphi
            nrm = math.sqrt(un * un + vn2 * vn2 + wn * wn)
            uo[i] = un / nrm; vo[i] = vn2 / nrm; wo[i] = wn / nrm
    return Eo, uo, vo, wo


def warmup_all():
    if HAVE_NUMBA:
        svt_sample(np.full(4, 2200.0), 1.5e-3)
        z = np.zeros(4); o = np.ones(4)
        free_gas_collision(np.full(4, 0.0253), z, z, o, 12.0, 1.5e-3)


@njit(fastmath=True, cache=True)
def _interp1(grid, y, e):
    """Clamped lin-lin interp of y on grid at scalar e (np.interp semantics)."""
    n = grid.size
    if e <= grid[0]:
        return y[0]
    if e >= grid[n - 1]:
        return y[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if grid[mid] <= e:
            lo = mid
        else:
            hi = mid
    t = (e - grid[lo]) / (grid[hi] - grid[lo])
    return y[lo] + t * (y[hi] - y[lo])


@njit(parallel=True, fastmath=True, cache=True)
def transport_sphere_thermal(x0, y0, z0, u0, v0, w0, E0,
                             R, grid, el_mac, cap_mac,
                             A, beta, cutoff, wcut, rsurv, max_coll):
    """Full single-nuclide thermal-sphere transport, one compiled history per
    lane (prange). Implicit capture + roulette; free-gas isotropic-CM elastic
    (valid where energies stay sub-cutoff, i.e. a thermal sphere). Replaces the
    entire Python event loop. Returns per-history (leaked_weight, escape_energy).
    """
    n = E0.size
    lw = np.zeros(n)
    le = np.zeros(n)
    Ap1 = A + 1.0
    for i in prange(n):
        x = x0[i]; y = y0[i]; z = z0[i]
        u = u0[i]; v = v0[i]; w = w0[i]
        E = E0[i]; wt = 1.0
        alive = True
        for _c in range(max_coll):
            sig_el = _interp1(grid, el_mac, E)
            sig_cap = _interp1(grid, cap_mac, E)
            sig_t = sig_el + sig_cap
            if sig_t <= 0.0:
                lw[i] = wt; le[i] = E; alive = False; break
            s = -math.log(1.0 - np.random.random()) / sig_t
            # distance to sphere surface from (x,y,z) along (u,v,w)
            k = x * u + y * v + z * w
            c = x * x + y * y + z * z - R * R
            disc = k * k - c
            if disc < 0.0:
                dist = 1e30
            else:
                sq = math.sqrt(disc)
                d1 = -k - sq; d2 = -k + sq
                dist = d2 if c < 0.0 else (d1 if d1 >= 0.0 else (d2 if d2 >= 0.0 else 1e30))
            if s > dist:
                x += dist * u; y += dist * v; z += dist * w
                lw[i] = wt; le[i] = E; alive = False; break
            x += s * u; y += s * v; z += s * w
            # roulette then implicit capture
            if wt < wcut:
                if np.random.random() < rsurv:
                    wt /= rsurv
                else:
                    alive = False; break
            wt *= sig_el / sig_t   # implicit capture (scatter fraction survives)
            # free-gas isotropic-CM elastic collision (inline)
            v_n = math.sqrt(2.0 * E * _EV / _M_N)
            p_first = 2.0 / (_SQRT_PI * v_n * beta + 2.0)
            while True:
                r1 = np.random.random(); r2 = np.random.random(); r3 = np.random.random()
                r4 = np.random.random(); r5 = np.random.random()
                if np.random.random() < p_first:
                    xx = math.sqrt(-math.log((1.0 - r1) * (1.0 - r2)))
                else:
                    cc = math.cos(0.5 * math.pi * r3)
                    xx = math.sqrt(-math.log(1.0 - r1) - math.log(1.0 - r2) * cc * cc)
                v_t = xx / beta; mu_t = 2.0 * r4 - 1.0
                acc = math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t)
                if r5 < acc:
                    break
            phi_t = 2.0 * math.pi * np.random.random()
            sin_t = math.sqrt(max(0.0, 1.0 - mu_t * mu_t)) if v_t > 0.0 else 0.0
            vtx = v_t * sin_t * math.cos(phi_t) if v_t > 0.0 else 0.0
            vty = v_t * sin_t * math.sin(phi_t) if v_t > 0.0 else 0.0
            vtz = v_t * mu_t if v_t > 0.0 else 0.0
            vcx = A * vtx / Ap1; vcy = A * vty / Ap1; vcz = (v_n + A * vtz) / Ap1
            Vx = -vcx; Vy = -vcy; Vz = v_n - vcz
            Vmag = math.sqrt(Vx * Vx + Vy * Vy + Vz * Vz)
            mu_cm = 2.0 * np.random.random() - 1.0
            phi_cm = 2.0 * math.pi * np.random.random()
            sin_cm = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
            vpx = Vmag * sin_cm * math.cos(phi_cm) + vcx
            vpy = Vmag * sin_cm * math.sin(phi_cm) + vcy
            vpz = Vmag * mu_cm + vcz
            sp2 = vpx * vpx + vpy * vpy + vpz * vpz
            E = max(1e-5, 0.5 * _M_N * sp2 / _EV)
            sp = math.sqrt(sp2)
            mu_lab = vpz / sp if sp > 0.0 else 1.0
            if mu_lab > 1.0: mu_lab = 1.0
            elif mu_lab < -1.0: mu_lab = -1.0
            phi = 2.0 * math.pi * np.random.random()
            cphi = math.cos(phi); sphi = math.sin(phi)
            st = math.sqrt(max(0.0, 1.0 - mu_lab * mu_lab))
            if abs(w) >= 0.999999:
                sign = 1.0 if w > 0.0 else -1.0
                u = st * cphi; v = st * sphi; w = sign * mu_lab
            else:
                dn = math.sqrt(max(1e-12, 1.0 - w * w)); rt = st / dn
                un = mu_lab * u + rt * (u * w * cphi - v * sphi)
                vv = mu_lab * v + rt * (v * w * cphi + u * sphi)
                wn = mu_lab * w - st * dn * cphi
                nn = math.sqrt(un * un + vv * vv + wn * wn)
                u = un / nn; v = vv / nn; w = wn / nn
    return lw, le
