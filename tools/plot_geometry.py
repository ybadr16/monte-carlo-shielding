"""Render a slice through a PyNeut CSG geometry, coloured by medium.

A validation table tells you the answer was right; it does not tell you the
model was the one you meant.  This renders the model itself, which is both a
figure for the documentation and the quickest way to catch a geometry that is
inside out, off centre, or missing a region because a priority was wrong.

Used as a library:

    from tools.plot_geometry import plot_slice
    plot_slice(media, plane="xy", bounds=(-1, 1, -1, 1))

or from the command line to rebuild the documentation figure:

    python tools/plot_geometry.py

Sampling is by point containment, exactly the test the transport uses to
resolve a medium, so what you see is what the tracker sees -- including the
priority ordering, which is why an overlap shows up as the region that
actually wins rather than as a blend.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from tools import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Axis index for each coordinate name, and the pair each slice plane spans.
_AXES = {"x": 0, "y": 1, "z": 2}
_PLANES = {"xy": ("x", "y", "z"), "xz": ("x", "z", "y"), "yz": ("y", "z", "x")}


def resolve_medium(media, x, y, z):
    """Index of the highest-priority medium containing a point, or -1.

    Mirrors the scalar tracker: strictly greater priority wins, so the first
    medium listed keeps a tie.
    """
    best, best_pri = -1, -np.inf
    for i, m in enumerate(media):
        if m.contains(x, y, z) and m.priority > best_pri:
            best, best_pri = i, m.priority
    return best


def medium_map(media, plane="xy", bounds=(-1, 1, -1, 1), offset=0.0, n=400):
    """Sample a plane and return ``(index_grid, extent)``.

    ``bounds`` is ``(lo_1, hi_1, lo_2, hi_2)`` over the two in-plane axes and
    ``offset`` fixes the third.  Grid points sit at pixel centres so a surface
    lying exactly on the sampling line does not bias a whole row.
    """
    if plane not in _PLANES:
        raise ValueError("plane must be one of {}".format(sorted(_PLANES)))
    a1, a2, a3 = _PLANES[plane]
    lo1, hi1, lo2, hi2 = bounds

    e1 = np.linspace(lo1, hi1, n + 1)
    e2 = np.linspace(lo2, hi2, n + 1)
    c1 = 0.5 * (e1[:-1] + e1[1:])
    c2 = 0.5 * (e2[:-1] + e2[1:])

    grid = np.full((n, n), -1, dtype=int)
    point = [0.0, 0.0, 0.0]
    point[_AXES[a3]] = offset
    i1, i2 = _AXES[a1], _AXES[a2]
    # Row index is the second axis so the array displays with axis 1 upward.
    for r, v2 in enumerate(c2):
        point[i2] = v2
        for c, v1 in enumerate(c1):
            point[i1] = v1
            grid[r, c] = resolve_medium(media, *point)
    return grid, (lo1, hi1, lo2, hi2)


def plot_slice(media, plane="xy", bounds=(-1, 1, -1, 1), offset=0.0, n=400,
               title=None, filename=None, subdir=None):
    """Draw a labelled material map of one slice through ``media``."""
    grid, extent = medium_map(media, plane, bounds, offset, n)
    a1, a2, a3 = _PLANES[plane]

    fig, ax = plt.subplots(figsize=(6.0, 5.4))

    present = [i for i in range(len(media)) if np.any(grid == i)]
    handles = []
    for i in present:
        colour, hatch = ps.medium_style(i)
        mask = np.where(grid == i, 1.0, np.nan)
        ax.imshow(mask, extent=extent, origin="lower", vmin=0, vmax=1,
                  cmap=_solid_cmap(colour), interpolation="nearest")
        if hatch is not None:
            # Texture is the secondary channel once the four validated hues
            # are used up, so identity never rests on colour alone.
            ax.contourf(np.where(grid == i, 1.0, 0.0), levels=[0.5, 1.5],
                        extent=extent, colors="none", hatches=[hatch])
        name = media[i].name or "medium {}".format(i)
        handles.append(Patch(facecolor=colour, hatch=hatch or "",
                             edgecolor=ps.SURFACE, label=name))

        # Direct label inside the region.  Three of the four hues sit under
        # the 3:1 contrast bar against the surface, so a visible label is
        # required rather than optional.
        spot = _label_point(grid == i, extent)
        if spot is not None:
            ax.text(spot[0], spot[1], name, ha="center", va="center",
                    fontsize=9, color=ps.INK_PRIMARY,
                    bbox=dict(boxstyle="round,pad=0.22", fc=ps.SURFACE,
                              ec="none", alpha=0.82))

    ax.set_xlabel("{} (cm)".format(a1))
    ax.set_ylabel("{} (cm)".format(a2))
    ax.set_title(title or "Material map, {} = {:g} cm".format(a3, offset))
    ax.set_aspect("equal")
    ax.grid(False)
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=9)

    if filename:
        return ps.save(fig, filename, subdir=subdir)
    return fig


def _solid_cmap(colour):
    """A colormap that paints every non-NaN cell one flat colour."""
    from matplotlib.colors import ListedColormap
    return ListedColormap([colour])


# A label needs roughly this much clear space, as a fraction of the shorter
# image side, before it will fit inside the region without spilling over a
# boundary.  Thin shells -- a pin cell's helium gap is a few pixels across --
# fall below it and are identified by the legend alone.
_LABEL_CLEARANCE = 0.045


def _label_point(mask, extent):
    """Where to put a region's label: its point of greatest clearance.

    The centroid is the obvious choice and the wrong one -- every annulus in
    a pin cell has its centroid at the origin, so all four labels would land
    on top of each other in the middle of the fuel.  The point that maximises
    distance to the region's own boundary sits comfortably inside each shell
    instead.  Returns ``None`` when the region is too thin to label.
    """
    if not mask.any():
        return None
    try:
        from scipy.ndimage import distance_transform_edt
    except ImportError:
        return None            # legend still identifies the region

    # Pad so a region touching the frame is measured against the frame too,
    # otherwise the outer moderator would want its label in a corner.
    dist = distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1]
    r, c = np.unravel_index(np.argmax(dist), dist.shape)
    if dist[r, c] < _LABEL_CLEARANCE * min(mask.shape):
        return None

    lo1, hi1, lo2, hi2 = extent
    x = lo1 + (c + 0.5) / mask.shape[1] * (hi1 - lo1)
    y = lo2 + (r + 0.5) / mask.shape[0] * (hi2 - lo2)
    return x, y


# --------------------------------------------------------------------------
# Documentation figure: the BEAVRS pin cell the validation suite runs
# --------------------------------------------------------------------------
# Dimensions match Validation/OpenMC_Comparison/validate_pincell.py.
PITCH = 1.25984
R_PELLET = 0.39218
R_INNER_CLAD = 0.40005
R_OUTER_CLAD = 0.45720


def pincell_media():
    """The four-region BEAVRS HZP pin cell, geometry only (no materials).

    Containment and priority are all the renderer needs, so this skips the
    nuclide compositions and stays runnable without the ENDF/B library.
    """
    from src.medium import Region, Cylinder, Box

    half = PITCH / 2.0
    pellet = Cylinder("z", R_PELLET, (0, 0, 0))
    inner = Cylinder("z", R_INNER_CLAD, (0, 0, 0))
    outer = Cylinder("z", R_OUTER_CLAD, (0, 0, 0))
    cube = Box(-half, half, -half, half, -half, half,
               boundary_type="reflective")

    inside_inner = Region(surfaces=[inner, cube], operation="intersection")
    inside_outer = Region(surfaces=[outer, cube], operation="intersection")

    return [
        Region(surfaces=[pellet, cube], operation="intersection",
               name="UO$_2$ fuel", priority=4),
        Region(surfaces=[inside_inner, pellet], operation="difference",
               name="He gap", priority=3),
        Region(surfaces=[inside_outer, inner], operation="difference",
               name="Zircaloy-4", priority=2),
        Region(surfaces=[cube, outer], operation="difference",
               name="borated water", priority=1),
    ]


def main():
    ps.use_style()
    half = PITCH / 2.0
    print("Rendering geometry figures...")
    plot_slice(
        pincell_media(), plane="xy",
        bounds=(-half, half, -half, half), offset=0.0, n=500,
        title="BEAVRS HZP pin cell, sampled by PyNeut's own containment test",
        filename="geometry_pincell.png",
    )
    print("Done.")


if __name__ == "__main__":
    main()
