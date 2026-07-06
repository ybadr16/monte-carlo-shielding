"""Stage-1 fully-njit fixed-source transport (single-nuclide, full CSG).

Compiles a single-nuclide fixed-source problem into flat arrays and transports
every history in one @njit(parallel) loop over the compiled geometry
(numba_geometry): free-flight sampling, CSG boundary crossing with reflection,
free-gas (sub-cutoff) and static-tabulated-angular (above-cutoff) elastic
scattering, discrete-level inelastic (MT 51-90), implicit capture and roulette.
No (n,xn)/fission yet -- those stay on the vector engine. Numba's per-thread RNG
is used, so results match the scalar kernel statistically, not stream-for-stream.
The scalar kernel remains the reference/oracle.
"""
import math

import numpy as np

from .numba_geometry import (
    CompiledGeometry, surf_normal, resolve_medium, nearest_crossing,
    contains_medium, HAVE_NUMBA,
)

try:
    from numba import njit, prange
except Exception:  # pragma: no cover
    def njit(*a, **k):
        return (lambda f: f) if not a else a[0]
    prange = range

_M_N = 1.674927471e-27
_EV = 1.60217663e-19
_SQRT_PI = math.sqrt(math.pi)


@njit(fastmath=True, cache=True)
def _interp1(grid, y, e):
    n = grid.size
    if e <= grid[0]:
        return y[0]
    if e >= grid[n - 1]:
        return y[n - 1]
    lo = 0; hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if grid[mid] <= e:
            lo = mid
        else:
            hi = mid
    t = (e - grid[lo]) / (grid[hi] - grid[lo])
    return y[lo] + t * (y[hi] - y[lo])


@njit(fastmath=True, cache=True)
def _sample_mu_tab(E, eg, off, mus, cdfs):
    """CM elastic cosine from the tabulated angular data; statistical
    interpolation between incident-energy blocks then CDF inversion. Mirrors
    AngularDistribution.sample_mu."""
    ne = eg.size
    idx = 0
    if E <= eg[0]:
        idx = 0
    elif E >= eg[ne - 1]:
        idx = ne - 2
    else:
        lo = 0; hi = ne - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if eg[mid] <= E:
                lo = mid
            else:
                hi = mid
        idx = lo
    f = (E - eg[idx]) / (eg[idx + 1] - eg[idx])
    if f < 0.0: f = 0.0
    elif f > 1.0: f = 1.0
    block = idx if np.random.random() > f else idx + 1
    s = off[block]; e = off[block + 1]
    r = np.random.random()
    # invert the block CDF (linear in cdf<->mu), matching np.interp
    if r <= cdfs[s]:
        return mus[s]
    if r >= cdfs[e - 1]:
        return mus[e - 1]
    lo = s; hi = e - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if cdfs[mid] <= r:
            lo = mid
        else:
            hi = mid
    c0 = cdfs[lo]; c1 = cdfs[hi]
    if c1 == c0:
        return mus[lo]
    return mus[lo] + (r - c0) / (c1 - c0) * (mus[hi] - mus[lo])


@njit(parallel=True, cache=True)
def transport_fixed(x0, y0, z0, u0, v0, w0, E0,
                    n_media, prim, prim_off, priority, is_void, bc,
                    tok, tok_off, stype, params,
                    grid, el_mac, in_mac, cap_mac,
                    ang_eg, ang_off, ang_mu, ang_cdf, has_ang,
                    inel_xs, inel_Q, n_inel,
                    A, beta, cutoff, wcut, rsurv, eps, max_events):
    n = x0.size
    lw = np.zeros(n)
    le = np.zeros(n)
    Ap1 = A + 1.0
    for i in prange(n):
        x = x0[i]; y = y0[i]; z = z0[i]
        u = u0[i]; v = v0[i]; w = w0[i]
        E = E0[i]; wt = 1.0
        for _ev in range(max_events):
            m = resolve_medium(n_media, priority, tok, tok_off, stype, params, x, y, z)
            if m < 0:
                lw[i] = wt; le[i] = E; break
            dist, s_surf = nearest_crossing(n_media, prim, prim_off, priority,
                                            tok, tok_off, stype, params,
                                            x, y, z, u, v, w)
            if s_surf < 0:
                lw[i] = wt; le[i] = E; break

            crossed = False
            if is_void[m] == 1:
                crossed = True
                s = 1e30
            else:
                sig_el = _interp1(grid, el_mac, E)
                sig_in = _interp1(grid, in_mac, E)
                sig_cap = _interp1(grid, cap_mac, E)
                sig_t = sig_el + sig_in + sig_cap
                if sig_t <= 0.0:
                    lw[i] = wt; le[i] = E; break
                s = -math.log(1.0 - np.random.random()) / sig_t
                if s > dist:
                    crossed = True

            if crossed:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:  # reflective
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u = u - 2.0 * dot * nx; v = v - 2.0 * dot * ny; w = w - 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue

            # collision
            x += s * u; y += s * v; z += s * w
            if not contains_medium(m, tok, tok_off, stype, params, x, y, z):
                continue
            if wt < wcut:
                if np.random.random() < rsurv:
                    wt /= rsurv
                else:
                    break
            p_scat = (sig_el + sig_in) / sig_t
            wt *= p_scat
            if wt <= 0.0:
                break

            p_inel = sig_in / (sig_el + sig_in) if (sig_el + sig_in) > 0.0 else 0.0
            if n_inel > 0 and np.random.random() < p_inel:
                # discrete inelastic: pick MT by cumulative xs, static two-body w/ Q
                tot = 0.0
                for j in range(n_inel):
                    tot += _interp1(grid, inel_xs[j], E)
                roll = np.random.random() * tot
                acc = 0.0; Q = 0.0
                for j in range(n_inel):
                    acc += _interp1(grid, inel_xs[j], E)
                    if roll <= acc:
                        Q = inel_Q[j]; break
                v_n = math.sqrt(2.0 * E * _EV / _M_N)
                vcz = v_n / Ap1
                Vz = v_n - vcz
                E_cm = 0.5 * _M_N * Vz * Vz / _EV
                tot_cm = E_cm * Ap1 / A
                if tot_cm <= Q:
                    mag = Vz            # cannot pay Q -> elastic-like
                else:
                    E_cm_p = (tot_cm - Q) * A / Ap1
                    mag = math.sqrt(2.0 * E_cm_p * _EV / _M_N)
                mu_cm = 2.0 * np.random.random() - 1.0
                phi_cm = 2.0 * math.pi * np.random.random()
                sin_cm = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
                vpx = mag * sin_cm * math.cos(phi_cm)
                vpy = mag * sin_cm * math.sin(phi_cm)
                vpz = mag * mu_cm + vcz
                sp2 = vpx * vpx + vpy * vpy + vpz * vpz
                E = max(1e-5, 0.5 * _M_N * sp2 / _EV)
                sp = math.sqrt(sp2)
                mu_lab = vpz / sp if sp > 0.0 else 1.0
            else:
                # elastic
                if E < cutoff:
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
                        acc2 = math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t)
                        if r5 < acc2:
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
                else:
                    # static target, tabulated CM angular
                    if has_ang == 1:
                        mu_cm = _sample_mu_tab(E, ang_eg, ang_off, ang_mu, ang_cdf)
                    else:
                        mu_cm = 2.0 * np.random.random() - 1.0
                    term = A * A + 1.0 + 2.0 * A * mu_cm
                    E = max(1e-5, E * term / (Ap1 * Ap1))
                    mu_lab = (1.0 + A * mu_cm) / math.sqrt(term)
            if mu_lab > 1.0: mu_lab = 1.0
            elif mu_lab < -1.0: mu_lab = -1.0
            # rotate (u,v,w) by mu_lab
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


# ---------------------------------------------------------------------------
# Python driver: compile a single-nuclide problem and run
# ---------------------------------------------------------------------------

THERMAL_CUTOFF = 10.0


def _angular_arrays(reader, element):
    """Flat (energy_grid, offsets, mu, cdf, has_ang) for MT=2 elastic angular."""
    from .angular_distribution import AngularDistribution
    emap = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
    ad = AngularDistribution(reader.base_path, emap.get(element, element), mt=2)
    ad.load()
    if not ad.loaded:
        return (np.zeros(1), np.zeros(2, np.int64), np.zeros(1), np.zeros(1), 0)
    eg = np.asarray(ad.energy_grid, dtype=float)
    off = np.asarray(ad.offsets, dtype=np.int64)
    # offsets has one entry per energy; append the total length as the final bound
    off = np.concatenate([off, np.array([ad.mu_data.shape[1]], dtype=np.int64)])
    return eg, off, np.ascontiguousarray(ad.mu_data[0]), np.ascontiguousarray(ad.mu_data[2]), 1


def _inelastic_arrays(reader, element, grid, N):
    """(inel_xs[n,len(grid)] macroscopic, inel_Q[n]) for discrete MT 51-90."""
    reader._scan_inelastic_mts(element)
    mts = [mt for mt in reader._available_inelastic_mts.get(element, []) if 51 <= mt <= 90]
    if not mts:
        return np.zeros((0, grid.size)), np.zeros(0), 0
    xs = np.zeros((len(mts), grid.size))
    Q = np.zeros(len(mts))
    for j, mt in enumerate(mts):
        d = reader.get_cross_section_data(element, mt)
        if d is None:
            continue
        xs[j] = np.interp(grid, d["energy"], d["xs"]) * 1e-24 * N
        Q[j] = abs(d["q_value"])
    return xs, Q, len(mts)


class NumbaProblem:
    """Compiled single-nuclide fixed-source problem for the njit engine."""

    def __init__(self, reader, mediums, element, N, A, beta):
        self.geo = CompiledGeometry(mediums)
        reader._build_fast_table(element)
        tbl = reader._fast_tables[element]
        sc = 1e-24 * N
        self.grid = np.ascontiguousarray(tbl["grid"])
        self.el = np.ascontiguousarray(tbl["el"] * sc)
        self.inl = np.ascontiguousarray(tbl["in"] * sc)
        self.cap = np.ascontiguousarray(tbl["cap"] * sc)
        self.ang = _angular_arrays(reader, element)
        self.inel_xs, self.inel_Q, self.n_inel = _inelastic_arrays(reader, element, self.grid, N)
        self.A = float(A); self.beta = float(beta)

    def run(self, source, wcut=1e-4, rsurv=0.1, eps=1e-6, max_events=100000):
        g = self.geo
        x = np.ascontiguousarray(source["x"], dtype=float)
        y = np.ascontiguousarray(source["y"], dtype=float)
        z = np.ascontiguousarray(source["z"], dtype=float)
        u = np.ascontiguousarray(source["u"], dtype=float)
        v = np.ascontiguousarray(source["v"], dtype=float)
        w = np.ascontiguousarray(source["w"], dtype=float)
        E = np.ascontiguousarray(source["energy"], dtype=float)
        lw, le = transport_fixed(
            x, y, z, u, v, w, E,
            g.n_media, g.prim, g.prim_off, g.priority, g.is_void, g.bc,
            g.tok, g.tok_off, g.stype, g.params,
            self.grid, self.el, self.inl, self.cap,
            self.ang[0], self.ang[1], self.ang[2], self.ang[3], self.ang[4],
            self.inel_xs, self.inel_Q, self.n_inel,
            self.A, self.beta, THERMAL_CUTOFF, wcut, rsurv, eps, max_events)
        return {"leaked_weight": lw, "escape_energy": le}


# ---------------------------------------------------------------------------
# Stage 2: multi-nuclide fixed-source (elastic + capture + isotope selection)
# ---------------------------------------------------------------------------

@njit(fastmath=True, cache=True)
def _sample_mu_tab_k(E, k, eg_off, bo_off, data_off, ang_eg, ang_bo, ang_mu, ang_cdf):
    """Tabulated CM cosine for nuclide k (flat ragged angular arrays)."""
    e0 = eg_off[k]; e1 = eg_off[k + 1]
    ne = e1 - e0
    if ne < 2:
        return 2.0 * np.random.random() - 1.0
    eg = ang_eg[e0:e1]
    if E <= eg[0]:
        idx = 0
    elif E >= eg[ne - 1]:
        idx = ne - 2
    else:
        lo = 0; hi = ne - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if eg[mid] <= E:
                lo = mid
            else:
                hi = mid
        idx = lo
    f = (E - eg[idx]) / (eg[idx + 1] - eg[idx])
    if f < 0.0: f = 0.0
    elif f > 1.0: f = 1.0
    block = idx if np.random.random() > f else idx + 1
    bo = ang_bo[bo_off[k]:bo_off[k + 1]]
    d0 = data_off[k]
    s = d0 + bo[block]; e = d0 + bo[block + 1]
    r = np.random.random()
    if r <= ang_cdf[s]:
        return ang_mu[s]
    if r >= ang_cdf[e - 1]:
        return ang_mu[e - 1]
    lo = s; hi = e - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ang_cdf[mid] <= r:
            lo = mid
        else:
            hi = mid
    c0 = ang_cdf[lo]; c1 = ang_cdf[hi]
    if c1 == c0:
        return ang_mu[lo]
    return ang_mu[lo] + (r - c0) / (c1 - c0) * (ang_mu[hi] - ang_mu[lo])


@njit(parallel=True, cache=True)
def transport_fixed_multi(x0, y0, z0, u0, v0, w0, E0,
                          n_media, prim, prim_off, priority, is_void, bc,
                          tok, tok_off, stype, params,
                          grid, goff, xs_el, xs_cap, nuc_A, nuc_beta,
                          ang_eg, ang_eg_off, ang_bo, ang_bo_off,
                          ang_mu, ang_cdf, ang_data_off, ang_has,
                          med_off, med_nuc, med_N,
                          cutoff, wcut, rsurv, eps, max_events):
    n = x0.size
    lw = np.zeros(n)
    le = np.zeros(n)
    for i in prange(n):
        x = x0[i]; y = y0[i]; z = z0[i]
        u = u0[i]; v = v0[i]; w = w0[i]
        E = E0[i]; wt = 1.0
        el_k = np.empty(64)
        for _ev in range(max_events):
            m = resolve_medium(n_media, priority, tok, tok_off, stype, params, x, y, z)
            if m < 0:
                lw[i] = wt; le[i] = E; break
            dist, s_surf = nearest_crossing(n_media, prim, prim_off, priority,
                                            tok, tok_off, stype, params,
                                            x, y, z, u, v, w)
            if s_surf < 0:
                lw[i] = wt; le[i] = E; break

            crossed = False
            ms = med_off[m]; me = med_off[m + 1]; nk = me - ms
            Sig_el = 0.0; Sig_cap = 0.0
            sig_t = 0.0
            if is_void[m] == 1:
                crossed = True; s = 1e30
            else:
                for j in range(nk):
                    gk = med_nuc[ms + j]; Nk = med_N[ms + j]
                    g0 = goff[gk]; g1 = goff[gk + 1]
                    ej = _interp1(grid[g0:g1], xs_el[g0:g1], E) * Nk
                    cj = _interp1(grid[g0:g1], xs_cap[g0:g1], E) * Nk
                    el_k[j] = ej
                    Sig_el += ej; Sig_cap += cj
                sig_t = Sig_el + Sig_cap
                if sig_t <= 0.0:
                    lw[i] = wt; le[i] = E; break
                s = -math.log(1.0 - np.random.random()) / sig_t
                if s > dist:
                    crossed = True

            if crossed:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u = u - 2.0 * dot * nx; v = v - 2.0 * dot * ny; w = w - 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue

            x += s * u; y += s * v; z += s * w
            if not contains_medium(m, tok, tok_off, stype, params, x, y, z):
                continue
            if wt < wcut:
                if np.random.random() < rsurv:
                    wt /= rsurv
                else:
                    break
            wt *= Sig_el / sig_t
            if wt <= 0.0:
                break
            roll = np.random.random() * Sig_el
            acc = 0.0; sel = 0
            for j in range(nk):
                acc += el_k[j]
                if roll <= acc:
                    sel = j; break
            gk = med_nuc[ms + sel]; A = nuc_A[gk]; beta = nuc_beta[gk]; Ap1 = A + 1.0
            if E < cutoff:
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
                    acc2 = math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t)
                    if r5 < acc2:
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
            else:
                if ang_has[gk] == 1:
                    mu_cm = _sample_mu_tab_k(E, gk, ang_eg_off, ang_bo_off,
                                             ang_data_off, ang_eg, ang_bo, ang_mu, ang_cdf)
                else:
                    mu_cm = 2.0 * np.random.random() - 1.0
                term = A * A + 1.0 + 2.0 * A * mu_cm
                E = max(1e-5, E * term / (Ap1 * Ap1))
                mu_lab = (1.0 + A * mu_cm) / math.sqrt(term)
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


class MultiNumbaProblem:
    """Stage-2 multi-nuclide fixed-source problem (elastic + capture)."""

    def __init__(self, reader, mediums):
        from .numba_geometry import CompiledGeometry
        from .numba_materials import compile_media
        self.geo = CompiledGeometry(mediums)
        self.tbl, self.med_off, self.med_nuc, self.med_N = compile_media(reader, mediums)

    def run(self, source, wcut=1e-4, rsurv=0.1, eps=1e-6, max_events=100000):
        g = self.geo; t = self.tbl
        x = np.ascontiguousarray(source["x"], float); y = np.ascontiguousarray(source["y"], float)
        z = np.ascontiguousarray(source["z"], float); u = np.ascontiguousarray(source["u"], float)
        v = np.ascontiguousarray(source["v"], float); w = np.ascontiguousarray(source["w"], float)
        E = np.ascontiguousarray(source["energy"], float)
        lw, le = transport_fixed_multi(
            x, y, z, u, v, w, E,
            g.n_media, g.prim, g.prim_off, g.priority, g.is_void, g.bc,
            g.tok, g.tok_off, g.stype, g.params,
            t.grid, t.goff, t.el, t.cap, t.A, t.beta,
            t.ang_eg, t.ang_eg_off, t.ang_bo, t.ang_bo_off,
            t.ang_mu, t.ang_cdf, t.ang_data_off, t.ang_has,
            self.med_off, self.med_nuc, self.med_N,
            THERMAL_CUTOFF, wcut, rsurv, eps, max_events)
        return {"leaked_weight": lw, "escape_energy": le}


# ---------------------------------------------------------------------------
# Stage 3a: single-nuclide analog criticality (fission + Watt + power iteration)
# ---------------------------------------------------------------------------

@njit(fastmath=True, cache=True)
def _watt_sample(a, b):
    """Watt fission energy (eV), Maxwellian-shift method (mirrors physics)."""
    r1 = np.random.random(); r2 = np.random.random(); r3 = np.random.random()
    w = a * (-math.log(r1) - math.log(r2) * math.cos(math.pi * r3 / 2.0) ** 2)
    a2b = a * a * b
    return w + a2b / 4.0 + (2.0 * np.random.random() - 1.0) * math.sqrt(a2b * w)


@njit(parallel=True, cache=True)
def keff_generation_single(x0, y0, z0, u0, v0, w0, E0,
                           n_media, prim, prim_off, priority, is_void, bc,
                           tok, tok_off, stype, params,
                           grid, el_mac, in_mac, cap_mac, fis_mac,
                           nu_grid, nu_val, watt_a, watt_b,
                           ang_eg, ang_off, ang_mu, ang_cdf, has_ang,
                           inel_xs, inel_Q, n_inel,
                           A, beta, cutoff, eps, max_events, maxf):
    """One fission generation, analog, single fissile nuclide. Returns per-history
    fission-production (collision estimator) and a preallocated fission bank."""
    n = x0.size
    fp = np.zeros(n)
    fb_x = np.zeros((n, maxf)); fb_y = np.zeros((n, maxf)); fb_z = np.zeros((n, maxf))
    fb_u = np.zeros((n, maxf)); fb_v = np.zeros((n, maxf)); fb_w = np.zeros((n, maxf))
    fb_E = np.zeros((n, maxf)); fb_c = np.zeros(n, np.int64)
    Ap1 = A + 1.0
    for i in prange(n):
        x = x0[i]; y = y0[i]; z = z0[i]
        u = u0[i]; v = v0[i]; w = w0[i]; E = E0[i]
        for _ev in range(max_events):
            m = resolve_medium(n_media, priority, tok, tok_off, stype, params, x, y, z)
            if m < 0:
                break
            dist, s_surf = nearest_crossing(n_media, prim, prim_off, priority,
                                            tok, tok_off, stype, params, x, y, z, u, v, w)
            if s_surf < 0:
                break
            if is_void[m] == 1:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u -= 2.0 * dot * nx; v -= 2.0 * dot * ny; w -= 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue
            sig_el = _interp1(grid, el_mac, E); sig_in = _interp1(grid, in_mac, E)
            sig_cap = _interp1(grid, cap_mac, E); sig_fis = _interp1(grid, fis_mac, E)
            sig_t = sig_el + sig_in + sig_cap + sig_fis
            if sig_t <= 0.0:
                break
            s = -math.log(1.0 - np.random.random()) / sig_t
            if s > dist:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u -= 2.0 * dot * nx; v -= 2.0 * dot * ny; w -= 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue
            x += s * u; y += s * v; z += s * w
            if not contains_medium(m, tok, tok_off, stype, params, x, y, z):
                continue
            # collision estimator of fission production
            if sig_fis > 0.0:
                fp[i] += _interp1(nu_grid, nu_val, E) * sig_fis / sig_t
            # analog reaction selection
            roll = np.random.random() * sig_t
            if roll < sig_el + sig_in:
                inelastic = roll >= sig_el
                if inelastic and n_inel > 0:
                    tot = 0.0
                    for j in range(n_inel):
                        tot += _interp1(grid, inel_xs[j], E)
                    rr = np.random.random() * tot; acc = 0.0; Q = 0.0
                    for j in range(n_inel):
                        acc += _interp1(grid, inel_xs[j], E)
                        if rr <= acc:
                            Q = inel_Q[j]; break
                    v_n = math.sqrt(2.0 * E * _EV / _M_N); vcz = v_n / Ap1; Vz = v_n - vcz
                    E_cm = 0.5 * _M_N * Vz * Vz / _EV; tot_cm = E_cm * Ap1 / A
                    mag = Vz if tot_cm <= Q else math.sqrt(2.0 * ((tot_cm - Q) * A / Ap1) * _EV / _M_N)
                    mu_cm = 2.0 * np.random.random() - 1.0
                    phic = 2.0 * math.pi * np.random.random(); sic = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
                    vpx = mag * sic * math.cos(phic); vpy = mag * sic * math.sin(phic); vpz = mag * mu_cm + vcz
                    sp2 = vpx * vpx + vpy * vpy + vpz * vpz
                    E = max(1e-5, 0.5 * _M_N * sp2 / _EV); sp = math.sqrt(sp2)
                    mu_lab = vpz / sp if sp > 0.0 else 1.0
                else:
                    if E < cutoff:
                        v_n = math.sqrt(2.0 * E * _EV / _M_N)
                        p_first = 2.0 / (_SQRT_PI * v_n * beta + 2.0)
                        while True:
                            r1 = np.random.random(); r2 = np.random.random(); r3 = np.random.random()
                            r4 = np.random.random(); r5 = np.random.random()
                            xx = (math.sqrt(-math.log((1.0 - r1) * (1.0 - r2))) if np.random.random() < p_first
                                  else math.sqrt(-math.log(1.0 - r1) - math.log(1.0 - r2) * math.cos(0.5 * math.pi * r3) ** 2))
                            v_t = xx / beta; mu_t = 2.0 * r4 - 1.0
                            if r5 < math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t):
                                break
                        phit = 2.0 * math.pi * np.random.random(); sit = math.sqrt(max(0.0, 1.0 - mu_t * mu_t)) if v_t > 0 else 0.0
                        vtx = v_t * sit * math.cos(phit) if v_t > 0 else 0.0
                        vty = v_t * sit * math.sin(phit) if v_t > 0 else 0.0
                        vtz = v_t * mu_t if v_t > 0 else 0.0
                        vcx = A * vtx / Ap1; vcy = A * vty / Ap1; vcz = (v_n + A * vtz) / Ap1
                        Vmag = math.sqrt(vcx * vcx + vcy * vcy + (v_n - vcz) ** 2)
                        mu_cm = 2.0 * np.random.random() - 1.0; phic = 2.0 * math.pi * np.random.random()
                        sic = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
                        vpx = Vmag * sic * math.cos(phic) + vcx; vpy = Vmag * sic * math.sin(phic) + vcy; vpz = Vmag * mu_cm + vcz
                        sp2 = vpx * vpx + vpy * vpy + vpz * vpz
                        E = max(1e-5, 0.5 * _M_N * sp2 / _EV); sp = math.sqrt(sp2)
                        mu_lab = vpz / sp if sp > 0.0 else 1.0
                    else:
                        mu_cm = _sample_mu_tab(E, ang_eg, ang_off, ang_mu, ang_cdf) if has_ang == 1 else 2.0 * np.random.random() - 1.0
                        term = A * A + 1.0 + 2.0 * A * mu_cm
                        E = max(1e-5, E * term / (Ap1 * Ap1)); mu_lab = (1.0 + A * mu_cm) / math.sqrt(term)
                if mu_lab > 1.0: mu_lab = 1.0
                elif mu_lab < -1.0: mu_lab = -1.0
                phi = 2.0 * math.pi * np.random.random(); cphi = math.cos(phi); sphi = math.sin(phi)
                sti = math.sqrt(max(0.0, 1.0 - mu_lab * mu_lab))
                if abs(w) >= 0.999999:
                    sign = 1.0 if w > 0.0 else -1.0
                    u = sti * cphi; v = sti * sphi; w = sign * mu_lab
                else:
                    dn = math.sqrt(max(1e-12, 1.0 - w * w)); rt = sti / dn
                    un = mu_lab * u + rt * (u * w * cphi - v * sphi)
                    vv = mu_lab * v + rt * (v * w * cphi + u * sphi)
                    wn = mu_lab * w - sti * dn * cphi
                    nn = math.sqrt(un * un + vv * vv + wn * wn)
                    u = un / nn; v = vv / nn; w = wn / nn
            elif roll < sig_el + sig_in + sig_cap:
                break  # capture
            else:
                # fission: bank nu prompt neutrons at (x,y,z), Watt energies, isotropic
                nu = int(_interp1(nu_grid, nu_val, E) + np.random.random())
                if nu > maxf:
                    nu = maxf
                for j in range(nu):
                    mu = 2.0 * np.random.random() - 1.0; ph = 2.0 * math.pi * np.random.random()
                    ss = math.sqrt(max(0.0, 1.0 - mu * mu))
                    fb_x[i, j] = x; fb_y[i, j] = y; fb_z[i, j] = z
                    fb_u[i, j] = ss * math.cos(ph); fb_v[i, j] = ss * math.sin(ph); fb_w[i, j] = mu
                    fb_E[i, j] = _watt_sample(watt_a, watt_b)
                fb_c[i] = nu
                break
    return fp, fb_x, fb_y, fb_z, fb_u, fb_v, fb_w, fb_E, fb_c


class SingleKeffProblem:
    """Stage-3a single-fissile-nuclide analog criticality via the njit generation
    kernel + a Python power-iteration driver. Elastic + discrete inelastic +
    capture + fission (nu-bar + Watt bank). Fast systems with significant
    (n,2n)/(n,3n) or continuum-inelastic multiplication need Stage 4 (not ported);
    thermal lattices do not."""

    def __init__(self, reader, mediums, element, N, A, beta):
        from .physics import watt_params_for
        self.geo = CompiledGeometry(mediums)
        reader._build_fast_table(element)
        tbl = reader._fast_tables[element]
        sc = 1e-24 * N
        self.grid = np.ascontiguousarray(tbl["grid"])
        self.el = np.ascontiguousarray(tbl["el"] * sc)
        self.inl = np.ascontiguousarray(tbl["in"] * sc)
        self.cap = np.ascontiguousarray(tbl["cap"] * sc)
        self.fis = np.ascontiguousarray(tbl["fis"] * sc)
        reader._load_nu(element)
        nue, nuv = reader._nu_cache[element]
        self.nu_e = np.ascontiguousarray(nue, float)
        self.nu_v = np.ascontiguousarray(nuv, float)
        self.wa, self.wb = watt_params_for(element)
        self.ang = _angular_arrays(reader, element)
        self.inel_xs, self.inel_Q, self.n_inel = _inelastic_arrays(reader, element, self.grid, N)
        self.A = float(A); self.beta = float(beta)

    def _generation(self, src, maxf=8, eps=1e-6, max_events=100000):
        g = self.geo
        return keff_generation_single(
            *(np.ascontiguousarray(src[k]) for k in ("x", "y", "z", "u", "v", "w", "energy")),
            g.n_media, g.prim, g.prim_off, g.priority, g.is_void, g.bc,
            g.tok, g.tok_off, g.stype, g.params,
            self.grid, self.el, self.inl, self.cap, self.fis,
            self.nu_e, self.nu_v, self.wa, self.wb,
            self.ang[0], self.ang[1], self.ang[2], self.ang[3], self.ang[4],
            self.inel_xs, self.inel_Q, self.n_inel,
            self.A, self.beta, THERMAL_CUTOFF, eps, max_events, maxf)

    def run_keff(self, initial_source, generations=90, inactive=30, seed=1, maxf=8):
        from .statistics import mean_sem
        r = np.random.default_rng(seed)
        src = {k: np.array(initial_source[k], float) for k in
               ("x", "y", "z", "u", "v", "w", "energy")}
        n = src["x"].size
        active, cyc = [], []
        for gen in range(generations):
            fp, fx, fy, fz, fu, fv, fw, fE, fc = self._generation(src, maxf)
            k = float(fp.sum()) / n
            mask = np.arange(maxf)[None, :] < fc[:, None]
            bx, by, bz, bu, bv, bw, bE = (a[mask] for a in (fx, fy, fz, fu, fv, fw, fE))
            produced = int(fc.sum())
            cyc.append((gen, k, produced / n, gen >= inactive))
            if gen >= inactive:
                active.append(k)
            if produced == 0:
                break
            pick = r.choice(produced, size=n, replace=produced < n)
            src = {"x": bx[pick], "y": by[pick], "z": bz[pick],
                   "u": bu[pick], "v": bv[pick], "w": bw[pick], "energy": bE[pick]}
        keff, ksem = mean_sem(active)
        return {"k_eff": keff, "k_sem": ksem, "active_k": active, "cycles": cyc}


@njit(fastmath=True, cache=True)
def _ss_right(arr, start, nb, r):
    b = 0
    while b < nb and arr[start + b] <= r:
        b += 1
    return b


@njit(fastmath=True, cache=True)
def _urr_sample(E, k, mic_el, mic_cap, mic_fis, r,
                urr_energy, urr_e_off, urr_tab_off, urr_nband, urr_interp, urr_mult,
                urr_cumul, urr_el, urr_fis, urr_cap):
    """Self-shielded (el, cap, fis) micro xs from nuclide k's URR probability
    table at energy E with random r. Mirrors reader._urr_micro_xs."""
    e0 = urr_e_off[k]; e1 = urr_e_off[k + 1]; ne = e1 - e0
    # locate bracket (searchsorted right - 1)
    i = 0
    if E <= urr_energy[e0]:
        i = 0
    elif E >= urr_energy[e1 - 1]:
        i = ne - 2
    else:
        lo = 0; hi = ne - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if urr_energy[e0 + mid] <= E:
                lo = mid
            else:
                hi = mid
        i = lo
    if i < 0: i = 0
    if i > ne - 2: i = ne - 2
    E0 = urr_energy[e0 + i]; E1 = urr_energy[e0 + i + 1]
    log_interp = urr_interp[k] == 5
    if log_interp and E0 > 0.0 and E1 > 0.0 and E > 0.0 and E1 != E0:
        f = math.log(E / E0) / math.log(E1 / E0)
    else:
        f = 0.0 if E1 == E0 else (E - E0) / (E1 - E0)
    if f < 0.0: f = 0.0
    elif f > 1.0: f = 1.0
    nb = urr_nband[k]; toff = urr_tab_off[k]
    row_i = toff + i * nb; row_i1 = toff + (i + 1) * nb
    b0 = _ss_right(urr_cumul, row_i, nb, r)
    if b0 > nb - 1: b0 = nb - 1
    b1 = _ss_right(urr_cumul, row_i1, nb, r)
    if b1 > nb - 1: b1 = nb - 1
    ve0 = urr_el[row_i + b0]; ve1 = urr_el[row_i1 + b1]
    vf0 = urr_fis[row_i + b0]; vf1 = urr_fis[row_i1 + b1]
    vc0 = urr_cap[row_i + b0]; vc1 = urr_cap[row_i1 + b1]
    if log_interp and ve0 > 0.0 and ve1 > 0.0:
        el = math.exp((1 - f) * math.log(ve0) + f * math.log(ve1))
    else:
        el = (1 - f) * ve0 + f * ve1
    if log_interp and vf0 > 0.0 and vf1 > 0.0:
        fis = math.exp((1 - f) * math.log(vf0) + f * math.log(vf1))
    else:
        fis = (1 - f) * vf0 + f * vf1
    if log_interp and vc0 > 0.0 and vc1 > 0.0:
        cap = math.exp((1 - f) * math.log(vc0) + f * math.log(vc1))
    else:
        cap = (1 - f) * vc0 + f * vc1
    if urr_mult[k] == 1:
        el *= mic_el; fis *= mic_fis; cap *= mic_cap
    return el, cap, fis


# ---------------------------------------------------------------------------
# Stage 3b: multi-nuclide analog criticality (isotope selection + URR + fission)
# ---------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def keff_generation_multi(x0, y0, z0, u0, v0, w0, E0,
                          n_media, prim, prim_off, priority, is_void, bc,
                          tok, tok_off, stype, params,
                          grid, goff, xs_el, xs_in, xs_cap, xs_fis, nuc_A, nuc_beta,
                          ang_eg, ang_eg_off, ang_bo, ang_bo_off, ang_mu, ang_cdf, ang_data_off, ang_has,
                          in_xs, in_Q, in_mt_off, in_xs_off,
                          nu_e, nu_v, nu_off, fissile, watt_a, watt_b,
                          urr_has, urr_emin, urr_emax, urr_interp, urr_mult,
                          urr_energy, urr_e_off, urr_tab_off, urr_nband,
                          urr_cumul, urr_el, urr_fis, urr_cap,
                          med_off, med_nuc, med_N, cutoff, eps, max_events, maxf):
    n = x0.size
    fp = np.zeros(n)
    fb_x = np.zeros((n, maxf)); fb_y = np.zeros((n, maxf)); fb_z = np.zeros((n, maxf))
    fb_u = np.zeros((n, maxf)); fb_v = np.zeros((n, maxf)); fb_w = np.zeros((n, maxf))
    fb_E = np.zeros((n, maxf)); fb_c = np.zeros(n, np.int64)
    for i in prange(n):
        x = x0[i]; y = y0[i]; z = z0[i]; u = u0[i]; v = v0[i]; w = w0[i]; E = E0[i]
        el_k = np.empty(64); in_k = np.empty(64); cap_k = np.empty(64)
        fis_k = np.empty(64); tot_k = np.empty(64)
        for _ev in range(max_events):
            m = resolve_medium(n_media, priority, tok, tok_off, stype, params, x, y, z)
            if m < 0:
                break
            dist, s_surf = nearest_crossing(n_media, prim, prim_off, priority,
                                            tok, tok_off, stype, params, x, y, z, u, v, w)
            if s_surf < 0:
                break
            ms = med_off[m]; me = med_off[m + 1]; nk = me - ms
            if is_void[m] == 1:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u -= 2.0 * dot * nx; v -= 2.0 * dot * ny; w -= 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue
            Sig_t = 0.0; fprod = 0.0
            for j in range(nk):
                gk = med_nuc[ms + j]; Nk = med_N[ms + j]
                g0 = goff[gk]; g1 = goff[gk + 1]
                e_ = _interp1(grid[g0:g1], xs_el[g0:g1], E)
                i_ = _interp1(grid[g0:g1], xs_in[g0:g1], E)
                c_ = _interp1(grid[g0:g1], xs_cap[g0:g1], E)
                f_ = _interp1(grid[g0:g1], xs_fis[g0:g1], E)
                if urr_has[gk] == 1 and urr_emin[gk] <= E <= urr_emax[gk]:
                    e_, c_, f_ = _urr_sample(E, gk, e_, c_, f_, np.random.random(),
                                             urr_energy, urr_e_off, urr_tab_off, urr_nband,
                                             urr_interp, urr_mult, urr_cumul, urr_el, urr_fis, urr_cap)
                ee = e_ * Nk; ii = i_ * Nk; cc = c_ * Nk; ff = f_ * Nk
                el_k[j] = ee; in_k[j] = ii; cap_k[j] = cc; fis_k[j] = ff
                tk = ee + ii + cc + ff; tot_k[j] = tk; Sig_t += tk
                if fissile[gk] == 1 and ff > 0.0:
                    fprod += _interp1(nu_e[nu_off[gk]:nu_off[gk + 1]],
                                      nu_v[nu_off[gk]:nu_off[gk + 1]], E) * ff
            if Sig_t <= 0.0:
                break
            s = -math.log(1.0 - np.random.random()) / Sig_t
            if s > dist:
                px = x + dist * u; py = y + dist * v; pz = z + dist * w
                if bc[s_surf] == 1:
                    nx, ny, nz = surf_normal(stype[s_surf], params[s_surf], px, py, pz)
                    dot = u * nx + v * ny + w * nz
                    u -= 2.0 * dot * nx; v -= 2.0 * dot * ny; w -= 2.0 * dot * nz
                    st = eps if (u * nx + v * ny + w * nz) > 0.0 else -eps
                    x = px + st * nx; y = py + st * ny; z = pz + st * nz
                else:
                    x = px + eps * u; y = py + eps * v; z = pz + eps * w
                continue
            x += s * u; y += s * v; z += s * w
            if not contains_medium(m, tok, tok_off, stype, params, x, y, z):
                continue
            fp[i] += fprod / Sig_t
            # isotope selection by total per-nuclide XS
            roll = np.random.random() * Sig_t; acc = 0.0; sel = 0
            for j in range(nk):
                acc += tot_k[j]
                if roll <= acc:
                    sel = j; break
            gk = med_nuc[ms + sel]; A = nuc_A[gk]; beta = nuc_beta[gk]; Ap1 = A + 1.0
            el_ = el_k[sel]; in_ = in_k[sel]; cap_ = cap_k[sel]; fis_ = fis_k[sel]
            tt = el_ + in_ + cap_ + fis_
            roll2 = np.random.random() * tt
            if roll2 < el_ + in_:
                do_inel = roll2 >= el_
                nmt = in_mt_off[gk + 1] - in_mt_off[gk]
                if do_inel and nmt > 0:
                    g0 = goff[gk]; g1 = goff[gk + 1]; base = in_mt_off[gk]
                    tot = 0.0
                    for jm in range(nmt):
                        xo = in_xs_off[base + jm]
                        tot += _interp1(grid[g0:g1], in_xs[xo:xo + (g1 - g0)], E)
                    rr = np.random.random() * tot; accm = 0.0; Q = 0.0
                    for jm in range(nmt):
                        xo = in_xs_off[base + jm]
                        accm += _interp1(grid[g0:g1], in_xs[xo:xo + (g1 - g0)], E)
                        if rr <= accm:
                            Q = in_Q[base + jm]; break
                    v_n = math.sqrt(2.0 * E * _EV / _M_N); vcz = v_n / Ap1; Vz = v_n - vcz
                    E_cm = 0.5 * _M_N * Vz * Vz / _EV; tot_cm = E_cm * Ap1 / A
                    mag = Vz if tot_cm <= Q else math.sqrt(2.0 * ((tot_cm - Q) * A / Ap1) * _EV / _M_N)
                    mu_cm = 2.0 * np.random.random() - 1.0; phic = 2.0 * math.pi * np.random.random()
                    sic = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
                    vpx = mag * sic * math.cos(phic); vpy = mag * sic * math.sin(phic); vpz = mag * mu_cm + vcz
                    sp2 = vpx * vpx + vpy * vpy + vpz * vpz; E = max(1e-5, 0.5 * _M_N * sp2 / _EV)
                    sp = math.sqrt(sp2); mu_lab = vpz / sp if sp > 0.0 else 1.0
                elif E < cutoff:
                    v_n = math.sqrt(2.0 * E * _EV / _M_N); p_first = 2.0 / (_SQRT_PI * v_n * beta + 2.0)
                    while True:
                        r1 = np.random.random(); r2 = np.random.random(); r3 = np.random.random()
                        r4 = np.random.random(); r5 = np.random.random()
                        xx = (math.sqrt(-math.log((1.0 - r1) * (1.0 - r2))) if np.random.random() < p_first
                              else math.sqrt(-math.log(1.0 - r1) - math.log(1.0 - r2) * math.cos(0.5 * math.pi * r3) ** 2))
                        v_t = xx / beta; mu_t = 2.0 * r4 - 1.0
                        if r5 < math.sqrt(v_n * v_n + v_t * v_t - 2.0 * v_n * v_t * mu_t) / (v_n + v_t):
                            break
                    phit = 2.0 * math.pi * np.random.random(); sit = math.sqrt(max(0.0, 1.0 - mu_t * mu_t)) if v_t > 0 else 0.0
                    vtx = v_t * sit * math.cos(phit) if v_t > 0 else 0.0; vty = v_t * sit * math.sin(phit) if v_t > 0 else 0.0
                    vtz = v_t * mu_t if v_t > 0 else 0.0
                    vcx = A * vtx / Ap1; vcy = A * vty / Ap1; vcz = (v_n + A * vtz) / Ap1
                    Vmag = math.sqrt(vcx * vcx + vcy * vcy + (v_n - vcz) ** 2)
                    mu_cm = 2.0 * np.random.random() - 1.0; phic = 2.0 * math.pi * np.random.random()
                    sic = math.sqrt(max(0.0, 1.0 - mu_cm * mu_cm))
                    vpx = Vmag * sic * math.cos(phic) + vcx; vpy = Vmag * sic * math.sin(phic) + vcy; vpz = Vmag * mu_cm + vcz
                    sp2 = vpx * vpx + vpy * vpy + vpz * vpz; E = max(1e-5, 0.5 * _M_N * sp2 / _EV)
                    sp = math.sqrt(sp2); mu_lab = vpz / sp if sp > 0.0 else 1.0
                else:
                    mu_cm = (_sample_mu_tab_k(E, gk, ang_eg_off, ang_bo_off, ang_data_off, ang_eg, ang_bo, ang_mu, ang_cdf)
                             if ang_has[gk] == 1 else 2.0 * np.random.random() - 1.0)
                    term = A * A + 1.0 + 2.0 * A * mu_cm
                    E = max(1e-5, E * term / (Ap1 * Ap1)); mu_lab = (1.0 + A * mu_cm) / math.sqrt(term)
                if mu_lab > 1.0: mu_lab = 1.0
                elif mu_lab < -1.0: mu_lab = -1.0
                phi = 2.0 * math.pi * np.random.random(); cphi = math.cos(phi); sphi = math.sin(phi)
                sti = math.sqrt(max(0.0, 1.0 - mu_lab * mu_lab))
                if abs(w) >= 0.999999:
                    sign = 1.0 if w > 0.0 else -1.0; u = sti * cphi; v = sti * sphi; w = sign * mu_lab
                else:
                    dn = math.sqrt(max(1e-12, 1.0 - w * w)); rt = sti / dn
                    un = mu_lab * u + rt * (u * w * cphi - v * sphi); vv = mu_lab * v + rt * (v * w * cphi + u * sphi)
                    wn = mu_lab * w - sti * dn * cphi; nn = math.sqrt(un * un + vv * vv + wn * wn)
                    u = un / nn; v = vv / nn; w = wn / nn
            elif roll2 < el_ + in_ + cap_:
                break
            else:
                nu = int(_interp1(nu_e[nu_off[gk]:nu_off[gk + 1]], nu_v[nu_off[gk]:nu_off[gk + 1]], E) + np.random.random())
                if nu > maxf: nu = maxf
                for jf in range(nu):
                    mu = 2.0 * np.random.random() - 1.0; ph = 2.0 * math.pi * np.random.random()
                    ss = math.sqrt(max(0.0, 1.0 - mu * mu))
                    fb_x[i, jf] = x; fb_y[i, jf] = y; fb_z[i, jf] = z
                    fb_u[i, jf] = ss * math.cos(ph); fb_v[i, jf] = ss * math.sin(ph); fb_w[i, jf] = mu
                    fb_E[i, jf] = _watt_sample(watt_a[gk], watt_b[gk])
                fb_c[i] = nu
                break
    return fp, fb_x, fb_y, fb_z, fb_u, fb_v, fb_w, fb_E, fb_c


class MultiKeffProblem:
    """Stage-3b multi-nuclide analog criticality: isotope selection + URR +
    per-nuclide elastic/discrete-inelastic/capture/fission over compiled CSG
    (reflective BCs supported). Thermal lattices; fast (n,2n) is Stage 4."""

    def __init__(self, reader, mediums):
        from .numba_geometry import CompiledGeometry
        from .numba_materials import compile_media
        self.geo = CompiledGeometry(mediums)
        self.t, self.med_off, self.med_nuc, self.med_N = compile_media(reader, mediums)

    def _generation(self, src, maxf=8, eps=1e-6, max_events=100000):
        g = self.geo; t = self.t
        return keff_generation_multi(
            *(np.ascontiguousarray(src[k]) for k in ("x", "y", "z", "u", "v", "w", "energy")),
            g.n_media, g.prim, g.prim_off, g.priority, g.is_void, g.bc,
            g.tok, g.tok_off, g.stype, g.params,
            t.grid, t.goff, t.el, t.inl, t.cap, t.fis, t.A, t.beta,
            t.ang_eg, t.ang_eg_off, t.ang_bo, t.ang_bo_off, t.ang_mu, t.ang_cdf, t.ang_data_off, t.ang_has,
            t.in_xs, t.in_Q, t.in_mt_off, t.in_xs_off,
            t.nu_e, t.nu_v, t.nu_off, t.fissile, t.watt_a, t.watt_b,
            t.urr_has, t.urr_emin, t.urr_emax, t.urr_interp, t.urr_mult,
            t.urr_energy, t.urr_e_off, t.urr_tab_off, t.urr_nband,
            t.urr_cumul, t.urr_el, t.urr_fis, t.urr_cap,
            self.med_off, self.med_nuc, self.med_N, THERMAL_CUTOFF, eps, max_events, maxf)

    def run_keff(self, initial_source, generations=90, inactive=30, seed=1, maxf=8):
        from .statistics import mean_sem
        r = np.random.default_rng(seed)
        src = {k: np.array(initial_source[k], float) for k in
               ("x", "y", "z", "u", "v", "w", "energy")}
        n = src["x"].size
        active, cyc = [], []
        for gen in range(generations):
            fp, fx, fy, fz, fu, fv, fw, fE, fc = self._generation(src, maxf)
            k = float(fp.sum()) / n
            mask = np.arange(maxf)[None, :] < fc[:, None]
            bx, by, bz, bu, bv, bw, bE = (a[mask] for a in (fx, fy, fz, fu, fv, fw, fE))
            produced = int(fc.sum())
            cyc.append((gen, k, produced / n, gen >= inactive))
            if gen >= inactive:
                active.append(k)
            if produced == 0:
                break
            pick = r.choice(produced, size=n, replace=produced < n)
            src = {"x": bx[pick], "y": by[pick], "z": bz[pick],
                   "u": bu[pick], "v": bv[pick], "w": bw[pick], "energy": bE[pick]}
        keff, ksem = mean_sem(active)
        return {"k_eff": keff, "k_sem": ksem, "active_k": active, "cycles": cyc}
