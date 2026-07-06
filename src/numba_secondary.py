"""Flat-array continuum secondary energy-angle distributions for the Numba engine.

Ports SecondaryDistribution (ENDF Law 61 correlated, Kalbach-Mann, Law 4
uncorrelated) into offset-indexed flat arrays and njit samplers that reproduce
sample_correlated / _sample_kalbach / sample_energy (unit-base outgoing-energy
interpolation, correlated mu sub-blocks, Kalbach angular law), plus CM->lab and
the Fermi-gas Maxwellian fallback. Used for (n,2n)/(n,3n)/continuum-inelastic
(MT 16/17/91) so the njit criticality engine matches the vector/scalar physics.
"""
import math

import numpy as np

from .cross_section_read import SecondaryDistribution

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*a, **k):
        return (lambda f: f) if not a else a[0]

_EV = 1.60217663e-19
_M_N = 1.674927471e-27

LAW_UNCORR = 0
LAW_CORR = 1
LAW_KALBACH = 2


class SecondaryTable:
    """Flattened continuum secondary distributions, one 'did' per (element, MT).

    Column data (eout/cdf/aux1/aux2) is concatenated across dids; per-did offset
    arrays recover each block. aux1/aux2 hold (r, a) for Kalbach or the two mu
    offset rows for correlated; unused for uncorrelated.
    """

    def __init__(self, base_path, pairs):
        # pairs: list of (element, mt); returns did per pair (or -1 if not loaded)
        self.did_of = {}
        law = []; frame = []
        ein = []; ein_off = [0]
        blk = []; blk_off = [0]
        eout = []; cdf = []; aux1 = []; aux2 = []; dat_off = [0]
        muval = []; mucdf = []; mu_off = [0]
        did = 0
        for (element, mt) in pairs:
            sd = SecondaryDistribution(base_path, element, mt)
            sd.load()
            if not sd.loaded:
                self.did_of[(element, mt)] = -1
                continue
            if sd.dist_type == "correlated":
                lw = LAW_CORR; data = sd.energy_out_data
            elif sd.dist_type == "kalbach-mann":
                lw = LAW_KALBACH; data = sd.distribution_data
            else:
                lw = LAW_UNCORR; data = sd.distribution_data
            law.append(lw)
            frame.append(1 if sd.ref_frame == "cm" else 0)
            ie = np.asarray(sd.incident_energy, float)
            ein.extend(ie.tolist()); ein_off.append(len(ein))
            blk.extend(np.asarray(sd.offsets, np.int64).tolist()); blk_off.append(len(blk))
            eout.extend(data[0].tolist()); cdf.extend(data[2].tolist())
            if lw == LAW_KALBACH:
                aux1.extend(data[3].tolist()); aux2.extend(data[4].tolist())
            elif lw == LAW_CORR:
                aux1.extend(data[3].tolist()); aux2.extend(data[4].tolist())
                m = sd.mu_data
                muval.extend(m[0].tolist()); mucdf.extend(m[2].tolist())
            else:  # uncorrelated: no aux
                aux1.extend([0.0] * data.shape[1]); aux2.extend([0.0] * data.shape[1])
            dat_off.append(len(eout)); mu_off.append(len(muval))
            self.did_of[(element, mt)] = did
            did += 1
        self.n_did = did
        self.law = np.array(law, np.int64); self.frame = np.array(frame, np.int64)
        self.ein = np.array(ein) if ein else np.zeros(1)
        self.ein_off = np.array(ein_off, np.int64)
        self.blk = np.array(blk, np.int64) if blk else np.zeros(1, np.int64)
        self.blk_off = np.array(blk_off, np.int64)
        self.eout = np.array(eout) if eout else np.zeros(1)
        self.cdf = np.array(cdf) if cdf else np.zeros(1)
        self.aux1 = np.array(aux1) if aux1 else np.zeros(1)
        self.aux2 = np.array(aux2) if aux2 else np.zeros(1)
        self.dat_off = np.array(dat_off, np.int64)
        self.muval = np.array(muval) if muval else np.zeros(1)
        self.mucdf = np.array(mucdf) if mucdf else np.zeros(1)
        self.mu_off = np.array(mu_off, np.int64)


@njit(fastmath=True, cache=True)
def _searchsorted(a, lo, hi, x):
    """np.searchsorted(a[lo:hi], x) offset back to absolute index in a."""
    while hi - lo > 0:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(fastmath=True, cache=True)
def _kalbach_mu_nj(a, r):
    if a < 1e-6:
        a = 1e-6
    elif a > 700.0:
        a = 700.0
    if np.random.random() <= r:
        xi = np.random.random()
        return math.log(xi * math.exp(a) + (1.0 - xi) * math.exp(-a)) / a
    T = (2.0 * np.random.random() - 1.0) * math.sinh(a)
    return math.log(T + math.sqrt(T * T + 1.0)) / a


@njit(fastmath=True, cache=True)
def _cm_to_lab_nj(E_cm, mu_cm, E_in, A):
    Ap1 = A + 1.0
    E_lab = E_cm + (E_in + 2.0 * mu_cm * Ap1 * math.sqrt(E_in * E_cm)) / (Ap1 * Ap1)
    if E_lab <= 0.0:
        return max(1e-5, E_cm), mu_cm
    mu_lab = mu_cm * math.sqrt(E_cm / E_lab) + math.sqrt(E_in / E_lab) / Ap1
    if mu_lab > 1.0: mu_lab = 1.0
    elif mu_lab < -1.0: mu_lab = -1.0
    return E_lab, mu_lab


@njit(fastmath=True, cache=True)
def _blk_bounds(blk, b0, sel, ne, datlen):
    bs = blk[b0 + sel]
    be = datlen if sel >= ne - 1 else blk[b0 + sel + 1]
    return bs, be


@njit(fastmath=True, cache=True)
def sample_secondary(did, E_in, A,
                     law, frame, ein, ein_off, blk, blk_off, dat_off,
                     eout, cdf, aux1, aux2, muval, mucdf, mu_off):
    """(E_out_lab, mu_lab, valid). Mirrors SecondaryDistribution.sample_correlated
    / _sample_kalbach / sample_energy including unit-base interpolation; applies
    CM->lab when the distribution frame is CM. valid=0 -> caller falls back."""
    e0 = ein_off[did]; e1 = ein_off[did + 1]; ne = e1 - e0
    if ne < 2:
        return 0.0, 0.0, 0
    idx = _searchsorted(ein, e0, e1, E_in) - 1 - e0
    if idx < 0: idx = 0
    elif idx > ne - 2: idx = ne - 2
    El = ein[e0 + idx]; Eh = ein[e0 + idx + 1]
    f = 0.0 if Eh == El else (E_in - El) / (Eh - El)
    if f < 0.0: f = 0.0
    elif f > 1.0: f = 1.0
    sel = idx if np.random.random() > f else idx + 1

    b0 = blk_off[did]
    dstart = dat_off[did]; datlen = dat_off[did + 1] - dstart
    bs, be = _blk_bounds(blk, b0, sel, ne, datlen)
    if bs >= be:
        return 0.0, 0.0, 0
    gs = dstart + bs; ge = dstart + be
    lw = law[did]

    r_E = np.random.random()
    k = _searchsorted(cdf, gs, ge, r_E) - gs
    if k < 0: k = 0
    elif k > (be - bs) - 1: k = (be - bs) - 1
    if k == 0:
        E_s = eout[gs]
    else:
        c0 = cdf[gs + k - 1]; c1 = cdf[gs + k]
        ea = eout[gs + k - 1]; eb = eout[gs + k]
        E_s = eb if abs(c1 - c0) < 1e-14 else ea + (r_E - c0) / (c1 - c0) * (eb - ea)

    # unit-base interpolation between incident blocks idx and idx+1
    e_min_sel = eout[gs]; e_max_sel = eout[ge - 1]
    span = e_max_sel - e_min_sel
    frac = 0.0 if span <= 0.0 else (E_s - e_min_sel) / span
    ls, le = _blk_bounds(blk, b0, idx, ne, datlen)
    hs, he = _blk_bounds(blk, b0, idx + 1, ne, datlen)
    lo_min = eout[dstart + ls]; lo_max = eout[dstart + le - 1]
    hi_min = eout[dstart + hs]; hi_max = eout[dstart + he - 1]
    e_min = lo_min + f * (hi_min - lo_min)
    e_max = lo_max + f * (hi_max - lo_max)
    E_s = e_min + frac * (e_max - e_min)

    gcol = gs + k
    if lw == LAW_KALBACH:
        mu = _kalbach_mu_nj(aux2[gcol], aux1[gcol])
        have_mu = 1
    elif lw == LAW_CORR:
        v3 = int(aux1[gcol]); v4 = int(aux2[gcol])
        mstart = -1; mlen = 0
        mo0 = mu_off[did]; molen = mu_off[did + 1] - mo0
        if v4 > 10:
            mstart = v4
            mlen = (int(aux2[gcol + 1]) - v4) if (gcol + 1) < dat_off[did + 1] else (molen - v4)
        elif v3 > 10:
            mstart = v3
            mlen = (int(aux1[gcol + 1]) - v3) if (gcol + 1) < dat_off[did + 1] else (molen - v3)
        if mstart < 0 or mlen <= 0:
            return 0.0, 0.0, 0
        ms = mo0 + mstart
        if abs(mucdf[ms + mlen - 1] - 1.0) > 0.01:
            return 0.0, 0.0, 0
        if mlen == 1:
            mu = muval[ms]
        else:
            rr = np.random.random()
            j = _searchsorted(mucdf, ms, ms + mlen, rr)
            if j <= ms:
                mu = muval[ms]
            elif j >= ms + mlen:
                mu = muval[ms + mlen - 1]
            else:
                c0 = mucdf[j - 1]; c1 = mucdf[j]
                mu = muval[j - 1] if c1 == c0 else muval[j - 1] + (rr - c0) / (c1 - c0) * (muval[j] - muval[j - 1])
        have_mu = 1
    else:
        mu = 2.0 * np.random.random() - 1.0
        have_mu = 0

    if have_mu == 0:
        mu = 2.0 * np.random.random() - 1.0
    if frame[did] == 1:
        E_s, mu = _cm_to_lab_nj(E_s, mu, E_in, A)
    return max(1e-5, E_s), mu, 1


@njit(fastmath=True, cache=True)
def maxwellian_nj(E_avail, A):
    """Fermi-gas Maxwellian evaporation energy (eV), capped at E_avail."""
    if E_avail <= 0.0:
        return 1e3
    T = math.sqrt((E_avail / 1e6) / (A / 8.0)) * 1e6
    r1 = np.random.random(); r2 = np.random.random(); r3 = np.random.random()
    E = T * (-math.log(r1) - math.log(r2) * math.cos(math.pi * r3 / 2.0) ** 2)
    return E if E < E_avail else E_avail
