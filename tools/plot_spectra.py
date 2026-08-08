"""Per-case leakage-spectrum panels for the OpenMC comparison suite.

The suite already computes everything this needs -- both spectra, their
batch-means uncertainties and the coarse-group chi-square -- and then throws
the spectra away and keeps one number per case.  This draws them instead, so
a POOR spectrum verdict can be *seen* rather than inferred from a chi-square
in a table column.

Called from ``_common.validate_case`` when ``run_all.py --plots`` is given.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from tools import plotstyle as ps
import matplotlib.pyplot as plt


def _coarse(values, sems, n_groups):
    """Sum a fine spectrum into ``n_groups`` groups, propagating the SEM.

    The same grouping ``compare_spectra`` scores, so the figure and the
    reported chi-square per degree of freedom describe the same histogram.
    """
    L = len(values)
    edges = np.linspace(0, L, n_groups + 1).astype(int)
    v = np.array([values[lo:hi].sum() for lo, hi in zip(edges[:-1], edges[1:])])
    s = np.array([np.sqrt(np.nansum(np.asarray(sems[lo:hi]) ** 2))
                  for lo, hi in zip(edges[:-1], edges[1:])])
    return v, s, edges


def _clip_err(values, sems, log_y):
    """Asymmetric error bars that stay on a log axis.

    Two ways a bar misdraws.  A group with no usable batch variance carries a
    NaN, which matplotlib renders as a bar spanning the whole frame -- the
    uncollided peak of a fast case does this, and the spike reads as a real
    feature.  And a one-sigma bar wider than the value itself reaches zero or
    below, which a log scale cannot draw.  Non-finite becomes no bar, and the
    lower arm is held just short of the value; the upper arm is untouched.
    """
    values = np.asarray(values, dtype=float)
    sems = np.nan_to_num(np.asarray(sems, dtype=float),
                         nan=0.0, posinf=0.0, neginf=0.0)
    if not log_y:
        return sems
    return np.vstack([np.minimum(sems, np.abs(values) * 0.95), sems])


def _energy_scale(e_max):
    """Pick eV / keV / MeV so the axis numbers stay legible."""
    if e_max >= 1e6:
        return 1e6, "MeV"
    if e_max >= 1e3:
        return 1e3, "keV"
    return 1.0, "eV"


def plot_case(case_name, e_bins, omc, pyn, chi2_dof=None, status=None,
              n_groups=20, outdir=None):
    """Draw one case: both spectra above, their ratio below.

    Two stacked panels sharing an energy axis rather than one panel with two
    y-scales -- a second scale would let any pair of curves be made to agree
    by choosing it well.
    """
    O, sO, edges = _coarse(np.asarray(omc["spectrum"]),
                           np.asarray(omc["spectrum_sem"]), n_groups)
    P, sP, _ = _coarse(np.asarray(pyn["spectrum"]),
                       np.asarray(pyn["spectrum_sem"]), n_groups)

    e_bins = np.asarray(e_bins)
    g_edges = e_bins[edges]
    mid = 0.5 * (g_edges[:-1] + g_edges[1:])
    scale, unit = _energy_scale(g_edges[-1])
    x = mid / scale

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.0, 5.6), sharex=True,
        gridspec_kw={"height_ratios": [2.6, 1.0], "hspace": 0.08})

    # Spectra span decades whenever the source is fast; a linear axis would
    # show the first group and nothing else.
    log_y = (np.any(O > 0)
             and np.nanmax(O) / np.nanmin(O[O > 0]) > 50)

    # Drawn as steps, not as a connected line: this is a histogram, and on a
    # fast case the last group holds the uncollided source peak two or three
    # decades above its neighbour -- joining the two with a straight segment
    # puts a near-vertical stroke on the figure that reads as a feature of
    # the spectrum rather than as the gap between two bins.
    ex = g_edges / scale
    ax.stairs(O, ex, color=ps.C_OPENMC, lw=2.0, label="OpenMC")
    ax.stairs(P, ex, color=ps.C_PYNEUT, lw=2.0, label="PyNeut", alpha=0.9)
    ax.errorbar(x, O, yerr=_clip_err(O, sO, log_y), color=ps.C_OPENMC,
                marker="o", ms=4, linestyle="none", capsize=0)
    ax.errorbar(x, P, yerr=_clip_err(P, sP, log_y), color=ps.C_PYNEUT,
                marker="s", ms=4, linestyle="none", capsize=0, alpha=0.9)
    ax.set_ylabel("Leakage per source neutron")
    ax.legend(loc="best", fontsize=9.5)
    if log_y:
        ax.set_yscale("log")

    title = case_name
    if chi2_dof is not None:
        title += "   ·   spectrum $\\chi^2$/dof = {:.2f}".format(chi2_dof)
        if status:
            title += " ({})".format(status)
    ax.set_title(title)

    # Ratio panel. Only groups where OpenMC carries signal are meaningful.
    good = O > 0
    ratio = np.full_like(O, np.nan, dtype=float)
    rerr = np.full_like(O, np.nan, dtype=float)
    ratio[good] = P[good] / O[good]
    # Both codes contribute statistical error to the ratio.
    rerr[good] = ratio[good] * np.sqrt((sP[good] / np.maximum(P[good], 1e-300)) ** 2
                                       + (sO[good] / O[good]) ** 2)

    axr.axhline(1.0, color=ps.AXISLINE, lw=1.4)
    axr.axhspan(0.9, 1.1, color=ps.STATUS_GOOD, alpha=0.09, lw=0)
    axr.errorbar(x, ratio, yerr=rerr, color=ps.INK_SECOND, marker="o", ms=4,
                 lw=1.4, capsize=0, linestyle="none")
    axr.set_ylabel("PyNeut / OpenMC")
    axr.set_xlabel("Escape energy ({})".format(unit))

    # The shared grid runs to at least 1 eV so thermal upscatter is captured,
    # which on a 0.0253 eV case leaves most of the axis empty. Trim to the
    # last group either code put a neutron in.
    signal = np.nonzero((O > 0) | (P > 0))[0]
    if len(signal):
        hi = ex[min(signal[-1] + 2, len(ex) - 1)]
        if hi > ex[0]:
            axr.set_xlim(ex[0], hi)

    finite = np.isfinite(ratio)
    if finite.any():
        span = np.nanmax(np.abs(ratio[finite] - 1.0))
        pad = max(0.15, min(span * 1.35, 1.0))
        axr.set_ylim(1 - pad, 1 + pad)

    name = "spectrum_{}.png".format(case_name)
    return ps.save(fig, name, subdir=outdir or "cases")


def configure():
    """Apply the shared style once per process."""
    ps.use_style()
