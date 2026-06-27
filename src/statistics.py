"""Batch-means uncertainties for Monte Carlo tallies.

Histories are split into n_batches contiguous batches, the tally is formed once
per batch, and the result is the batch mean with the standard error of the mean,
SEM = sqrt(Σ_b (x_b - mean)^2 / (B(B-1))). Ratio estimators (e.g. weighted-mean
escape energy) are formed within each batch before the batches are combined.
"""
import numpy as np


def batch_edges(n_items, n_batches):
    """Contiguous, near-equal batch boundaries over ``range(n_items)``."""
    n_batches = max(1, min(n_batches, n_items))
    return np.linspace(0, n_items, n_batches + 1).astype(int)


def mean_sem(batch_values):
    """Mean and standard error of the mean of a set of batch estimates.

    Batches that produced no estimate (``NaN``) are dropped.  With fewer than
    two valid batches the SEM is undefined and returned as ``nan``.
    """
    vals = np.asarray(batch_values, dtype=float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return 0.0, float('nan')
    mean = float(vals.mean())
    if vals.size < 2:
        return mean, float('nan')
    sem = float(vals.std(ddof=1) / np.sqrt(vals.size))
    return mean, sem


def escape_statistics(leaked_weight, escape_energy, n_source, n_batches=10,
                      e_bins=None):
    """Batch-means leakage fraction and weighted-average escape energy.

    leaked_weight and escape_energy span every source history (0 where no leak).
    If e_bins is given the escape energy is taken from a (energy, weight)
    histogram on bin midpoints — the same estimator OpenMC uses, so the binning
    bias cancels in a comparison. Returns leakage, leakage_sem, avg_energy,
    avg_energy_sem.
    """
    w = np.asarray(leaked_weight, dtype=float)
    e = np.asarray(escape_energy, dtype=float)
    if w.size != n_source or e.size != n_source:
        raise ValueError(
            f"leaked_weight/escape_energy must have length n_source={n_source}; "
            f"got {w.size} and {e.size}")

    if e_bins is not None:
        e_bins = np.asarray(e_bins, dtype=float)
        e_mid = 0.5 * (e_bins[:-1] + e_bins[1:])

    edges = batch_edges(n_source, n_batches)
    leak_batches = []
    energy_batches = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        n_b = hi - lo
        if n_b <= 0:
            continue
        wb = w[lo:hi]
        eb = e[lo:hi]
        leak_batches.append(wb.sum() / n_b)
        if e_bins is not None:
            hist, _ = np.histogram(eb, bins=e_bins, weights=wb)
            hsum = hist.sum()
            energy_batches.append((e_mid * hist).sum() / hsum if hsum > 0 else np.nan)
        else:
            wsum = wb.sum()
            energy_batches.append((eb * wb).sum() / wsum if wsum > 0 else np.nan)

    leakage, leakage_sem = mean_sem(leak_batches)
    avg_energy, avg_energy_sem = mean_sem(energy_batches)
    return {
        'leakage': leakage,
        'leakage_sem': leakage_sem,
        'avg_energy': avg_energy,
        'avg_energy_sem': avg_energy_sem,
    }


def ratio_mean_uncertainty(values, weights, weight_sigmas):
    """Weighted mean E = Σ Eᵢwᵢ / Σwᵢ and its error when the weights carry
    independent uncertainties: σ_E = sqrt(Σ (Eᵢ-E)² σᵢ²) / Σwᵢ. Used for the
    OpenMC leakage spectrum. Returns (mean, sigma); (0, nan) if Σwᵢ <= 0."""
    Ei = np.asarray(values, dtype=float)
    wi = np.asarray(weights, dtype=float)
    si = np.asarray(weight_sigmas, dtype=float)
    B = wi.sum()
    if B <= 0:
        return 0.0, float('nan')
    E = float((Ei * wi).sum() / B)
    var = float((((Ei - E) ** 2) * si ** 2).sum()) / (B ** 2)
    return E, float(np.sqrt(max(0.0, var)))
