"""Shared plotting style for every figure the repository ships.

The point of putting this in one place is that the README gallery, the
per-case validation panels and the geometry renders should read as one
system rather than as three scripts that each reached for matplotlib's
defaults.

Nothing here is imported by ``src``.  PyNeut's transport core depends on
numpy and h5py alone, and that stays true -- matplotlib is a tool-side
dependency only.

Colour choices are not free-hand.  The categorical hues below were checked
with a contrast/colour-vision validator on the *all-pairs* list (every pair
on screen at once, which is the situation a material map creates, as opposed
to a line chart where only neighbouring series need separating):

    blue #2a78d6, orange #eb6834, aqua #1baf7a, violet #4a3aa7
    worst pair, colour-deficient view : dE 9.2  (floor 8)
    worst pair, normal vision         : dE 16.3 (floor 15)

Four is the ceiling: every fifth hue tried drops a pair under one of those
floors.  Media past the fourth therefore reuse the four hues with a hatch
overlay, so identity never rests on colour alone.  Aqua sits at 2.74:1
against the surface, just under the 3:1 bar, which obliges visible labels --
the material maps direct-label every region, so that is satisfied.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------
# Figures are rendered on an opaque light surface rather than a transparent
# one: GitHub serves the README on a dark background for a good fraction of
# readers, and a transparent PNG puts black axis text on near-black there.
SURFACE = "#fcfcfb"

INK_PRIMARY = "#0b0b0b"
INK_SECOND = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXISLINE = "#c3c2b7"

# Categorical slots, in fixed order.  Assign by identity and never cycle:
# the same entity keeps its colour across every figure in the gallery.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]

# The two codes under comparison get slot 1 and slot 2 everywhere.
C_OPENMC = SERIES[0]
C_PYNEUT = SERIES[1]

# Reserved for verdicts (agreement bands, pass/warn/fail marks).  Kept apart
# from SERIES so a status colour can never be mistaken for a data series, and
# always paired with a label rather than standing on its own.
STATUS_GOOD = "#0ca30c"
STATUS_WARN = "#fab219"
STATUS_CRIT = "#d03b3b"

# Hatch patterns for media past the fourth (see module docstring).
HATCHES = [None, "///", "\\\\\\", "..."]


def medium_style(i):
    """Return ``(facecolor, hatch)`` for the i-th medium in a material map.

    Cycles colour fastest and hatch slowest, so the first four media are
    plain fills in the validated order and only a fifth medium onwards
    introduces texture.
    """
    return SERIES[i % len(SERIES)], HATCHES[(i // len(SERIES)) % len(HATCHES)]


# --------------------------------------------------------------------------
# Global style
# --------------------------------------------------------------------------
def use_style():
    """Apply the shared rcParams.  Call once at the top of a figure script."""
    plt.rcParams.update({
        "font.size": 11,
        "font.family": "sans-serif",
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        # Chrome recedes so the marks carry the figure.
        "axes.edgecolor": AXISLINE,
        "axes.labelcolor": INK_SECOND,
        "axes.titlecolor": INK_PRIMARY,
        "axes.titlesize": 12,
        "axes.titleweight": "medium",
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRIDLINE,
        "grid.linewidth": 0.8,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_SECOND,
        "ytick.labelcolor": INK_SECOND,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 8,
        # 150 for on-screen review, 200 on disk: enough for a README at 2x
        # without making the repository heavy.
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })
    # Only the two spines that carry an axis; the box adds nothing.
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Tracked, unlike paper_template/, so figures referenced from the README are
# actually present for anyone who clones the repository.
IMAGES = os.path.join(ROOT, "docs", "images")


def save(fig, name, subdir=None):
    """Write ``fig`` into docs/images[/subdir] and report the path."""
    out = IMAGES if subdir is None else os.path.join(IMAGES, subdir)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name)
    fig.savefig(path)
    plt.close(fig)
    print("  [figure] {}".format(os.path.relpath(path, ROOT)))
    return path
