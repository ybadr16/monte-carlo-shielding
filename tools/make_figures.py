"""Build the figures the README and the documentation page embed.

Every figure here is derived from data already committed to the repository,
so the script needs neither OpenMC nor the ENDF/B library and runs in a
second or two:

    python tools/make_figures.py

Figures that *do* need a live OpenMC run are produced elsewhere -- the
per-case spectrum panels come from ``run_all.py --plots``, and the geometry
renders from ``tools/plot_geometry.py``.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from tools import plotstyle as ps
import matplotlib.pyplot as plt

CSV_PATH = os.path.join(ps.ROOT, "Validation", "OpenMC_Comparison",
                        "validation_results.csv")

# The band the validation suite judges against: within 2 sigma is agreement,
# 2-4 sigma is a warning, beyond 4 sigma is a significant bias.
Z_OK = 2.0
Z_WARN = 4.0


def _load_rows():
    with open(CSV_PATH, newline="") as fh:
        return list(csv.DictReader(fh))


def _fmt_energy(mev):
    """Label a source energy the way a reactor physicist would say it."""
    ev = float(mev) * 1e6
    if ev < 1.0:
        return "{:.3g} eV".format(ev)
    if ev < 1e3:
        return "{:.3g} eV".format(ev)
    if ev < 1e6:
        return "{:.3g} keV".format(ev / 1e3)
    return "{:.4g} MeV".format(ev / 1e6)


# --------------------------------------------------------------------------
# Validation summary -- every case, both metrics, one image
# --------------------------------------------------------------------------
def fig_validation_summary():
    """One row per case; a marker per metric; the agreement band behind them.

    This is the figure that answers "is it validated?" without the reader
    having to parse a twenty-column table.  The horizontal axis is the
    z-score -- the difference between the two codes measured in units of the
    combined statistical uncertainty -- so cases with wildly different
    leakage fractions are still directly comparable.
    """
    rows = _load_rows()

    def zval(r, key):
        v = r.get(key)
        return float(v) if v not in (None, "", "None") else np.nan

    # Hardest case at the top: the eye lands on the worst result first, which
    # is the honest way round for a validation figure.
    rows.sort(key=lambda r: max(zval(r, "leakage_zscore"),
                                zval(r, "energy_zscore")))

    labels = ["{}  ·  {}  ·  {}".format(r["isotope"], r["geometry"],
                                        _fmt_energy(r["energy_MeV"]))
              for r in rows]
    z_leak = [zval(r, "leakage_zscore") for r in rows]
    z_energy = [zval(r, "energy_zscore") for r in rows]
    y = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(rows) + 1.9))

    # Bands first, so the markers sit on top of them.
    ax.axvspan(0, Z_OK, color=ps.STATUS_GOOD, alpha=0.10, lw=0)
    ax.axvspan(Z_OK, Z_WARN, color=ps.STATUS_WARN, alpha=0.12, lw=0)
    ax.axvline(Z_OK, color=ps.STATUS_GOOD, lw=1.2, alpha=0.55)
    ax.axvline(Z_WARN, color=ps.STATUS_WARN, lw=1.2, alpha=0.55)

    # A hairline leader per row ties the two markers to their label.
    for yi, a, b in zip(y, z_leak, z_energy):
        lo, hi = np.nanmin([a, b]), np.nanmax([a, b])
        ax.plot([lo, hi], [yi, yi], color=ps.GRIDLINE, lw=1.4, zorder=1)

    # A surface-coloured ring keeps the two markers legible where they land
    # on top of one another.
    common = dict(linestyle="none", markeredgecolor=ps.SURFACE,
                  markeredgewidth=1.4, zorder=3)
    ax.plot(z_leak, y, marker="o", color=ps.SERIES[2],
            label="escape fraction", **common)
    ax.plot(z_energy, y, marker="D", color=ps.SERIES[3], markersize=7,
            label="mean escape energy", **common)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    # Scale to the data, not to the 4 sigma threshold: every case sits below
    # 2.6, and stretching the axis out to 4 would squash the whole suite into
    # the left half of the figure to show empty band.
    z_all = [z for z in z_leak + z_energy if np.isfinite(z)]
    ax.set_xlim(0, max(2.5, max(z_all)) * 1.15)
    ax.set_xlabel("Disagreement with OpenMC  (combined standard errors, $z$)")
    ax.grid(axis="y", visible=False)

    # State the actual outcome rather than asserting a clean sweep -- one case
    # (the 14 MeV iron slab) does land in the warning band.
    n_ok = sum(1 for r in rows if max(zval(r, "leakage_zscore"),
                                      zval(r, "energy_zscore")) <= Z_OK)
    ax.set_title("PyNeut vs OpenMC: {} of {} cases agree within "
                 "$2\\sigma$".format(n_ok, len(rows)), pad=34)

    # Label the bands in place rather than adding two more legend entries.
    ax.text(Z_OK / 2, len(rows) - 0.5, "agreement", ha="center",
            va="center", fontsize=9, color=ps.INK_SECOND)
    ax.text(Z_OK + 0.08, len(rows) - 0.5, "warning", ha="left",
            va="center", fontsize=9, color=ps.INK_SECOND)

    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=2,
              fontsize=9.5, handletextpad=0.4, columnspacing=1.8)
    ps.save(fig, "validation_summary.png")


def main():
    ps.use_style()
    print("Building repository figures...")
    fig_validation_summary()
    print("Done.")


if __name__ == "__main__":
    main()
