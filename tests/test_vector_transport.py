"""Vectorized transport engine (src/vector_transport.py), milestone 1.

The array-form geometry and cross-section lookups are deterministic, so they
are held to near-bitwise agreement with the scalar implementations they
mirror. Full transport is checked physically: uncollided transmission through
a slab must follow the Beer-Lambert law within statistics.
"""
import os

import numpy as np
import pytest

from src.cross_section_read import CrossSectionReader
from src.geometry import calculate_nearest_crossing
from src.material import Material, Nuclide
from src.medium import Box, Cylinder, Plane, Region, Sphere
from src import vector_transport as vt

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")
HAVE_DATA = os.path.exists(os.path.join(ENDF, "neutron", "Pb208.h5"))

needs_data = pytest.mark.skipif(not HAVE_DATA, reason="ENDF data not present")


def _rays(n, seed, span=4.0):
    """Random points in a cube of half-width `span` with unit directions."""
    rng = np.random.default_rng(seed)
    x, y, z = (rng.uniform(-span, span, n) for _ in range(3))
    mu = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - mu ** 2)
    return x, y, z, s * np.cos(phi), s * np.sin(phi), mu


def _scalar_distances(surface, x, y, z, u, v, w):
    """Loop the scalar nearest_surface_method, mapping None -> inf."""
    out = np.empty(x.size)
    for i in range(x.size):
        d = surface.nearest_surface_method(
            float(x[i]), float(y[i]), float(z[i]),
            float(u[i]), float(v[i]), float(w[i]))
        out[i] = np.inf if d is None else d
    return out


@pytest.mark.parametrize("surface", [
    Plane(1.0, 2.0, -3.0, 4.0),
    Plane(0.0, 0.0, 1.0, 0.5),
    Sphere(center=(0.3, -0.2, 1.0), radius=2.5),
    Cylinder(axis="z", radius=1.7, center=(0.1, 0.2, 0.0)),
    Cylinder(axis="x", radius=2.2, center=(0.0, -0.4, 0.3)),
    Cylinder(axis="y", radius=0.9, center=(0.5, 0.0, -0.1)),
])
def test_surface_distances_match_scalar(surface):
    x, y, z, u, v, w = _rays(4000, seed=7)
    vec = vt.surface_distances(surface, x, y, z, u, v, w)
    ref = _scalar_distances(surface, x, y, z, u, v, w)
    hit_v, hit_r = np.isfinite(vec), np.isfinite(ref)
    assert np.array_equal(hit_v, hit_r)
    np.testing.assert_allclose(vec[hit_v], ref[hit_r], rtol=1e-12, atol=1e-12)


def test_contains_many_matches_scalar():
    inner = Sphere(center=(0.0, 0.0, 0.0), radius=1.5)
    box = Box(-3, 3, -3, 3, -3, 3)
    shell = Region([box, Region([inner])], operation="difference")
    union = Region([Sphere(center=(2, 0, 0), radius=1.0),
                    Sphere(center=(-2, 0, 0), radius=1.0)], operation="union")
    x, y, z, *_ = _rays(4000, seed=11)
    for region in (box, shell, union):
        vec = vt.contains_many(region, x, y, z)
        ref = np.array([region.contains(float(x[i]), float(y[i]), float(z[i]))
                        for i in range(x.size)])
        assert np.array_equal(vec, ref)


def test_nearest_crossings_match_scalar():
    sphere_region = Region([Sphere(center=(0, 0, 0), radius=2.0)], priority=1)
    box_region = Box(-3, 3, -3, 3, -3, 3)
    mediums = [sphere_region, box_region]
    x, y, z, u, v, w = _rays(2000, seed=13, span=2.5)

    dist, sid, surfaces = vt.nearest_crossings(mediums, x, y, z, u, v, w)

    for i in range(x.size):
        state = {"x": float(x[i]), "y": float(y[i]), "z": float(z[i])}
        _, _, ref_d, ref_surf = calculate_nearest_crossing(
            state, mediums, float(u[i]), float(v[i]), float(w[i]))
        if not np.isfinite(ref_d):
            assert not np.isfinite(dist[i])
            continue
        np.testing.assert_allclose(dist[i], ref_d, rtol=1e-12, atol=1e-12)
        assert surfaces[sid[i]] is ref_surf


@needs_data
def test_vector_xs_matches_scalar_reader():
    reader = CrossSectionReader(ENDF)
    mat = Material("Pb208", 11.35, 208, 206.19)
    nuc = Nuclide("Pb208", mat.number_density, 206.19)
    vx = vt.VectorXS(reader, [nuc])

    # above the 10 eV Doppler cutoff the scalar lookup is deterministic
    E = np.logspace(np.log10(11.0), 7.0, 500)
    Sel, Sin, Scap, Sfis, St = vx.macroscopic(E)
    for i, e in enumerate(E):
        s_el, s_in, s_cap, s_fis, s_t = reader.get_cross_sections(
            "Pb208", float(e), nuc.sampler, nuc.number_density,
            nuc.atomic_weight_ratio, rng=None)
        np.testing.assert_allclose(
            [Sel[i], Sin[i], Scap[i], Sfis[i], St[i]],
            [s_el, s_in, s_cap, s_fis, s_t], rtol=1e-12)


@needs_data
def test_beer_lambert_transmission():
    reader = CrossSectionReader(ENDF)
    mat = Material("Pb208", 11.35, 208, 206.19)
    nuc = Nuclide("Pb208", mat.number_density, 206.19)

    energy = 2.0e6
    Sigma_t = reader.get_cross_sections(
        "Pb208", energy, nuc.sampler, nuc.number_density,
        nuc.atomic_weight_ratio, rng=None)[4]
    thickness = 1.0 / Sigma_t  # one mean free path

    slab = Box(-50, 50, -50, 50, 0.0, thickness)
    slab.composition = [nuc]

    n = 100_000
    z0 = 1e-9
    source = {
        "x": np.zeros(n), "y": np.zeros(n), "z": np.full(n, z0),
        "u": np.zeros(n), "v": np.zeros(n), "w": np.ones(n),
        "energy": np.full(n, energy),
    }
    res = vt.run_streaming(reader, [slab], source, np.random.default_rng(42))

    expected = np.exp(-Sigma_t * (thickness - z0))
    sigma = np.sqrt(expected * (1 - expected) / n)
    assert abs(res["leakage"] - expected) < 4 * sigma
    assert res["escaped"].sum() + res["collided"].sum() == n


@needs_data
def test_reflective_plane_returns_beam():
    """A beam fired at a reflective wall through a void must come back out the
    open face: exercises the vectorized reflection + inward-normal step-off."""
    reader = CrossSectionReader(ENDF)
    box = Box(-1, 1, -1, 1, 0.0, 5.0)
    # make the far z-plane reflective, leave the rest transmitting
    for p in box.surfaces:
        if p.C == 1.0:  # the z <= z_max plane
            p.boundary_type = "reflective"
    box.is_void = True

    n = 1000
    source = {
        "x": np.zeros(n), "y": np.zeros(n), "z": np.full(n, 1e-3),
        "u": np.zeros(n), "v": np.zeros(n), "w": np.ones(n),
        "energy": np.full(n, 1.0e6),
    }
    res = vt.run_streaming(reader, [box], source, np.random.default_rng(1))
    assert res["escaped"].all()
    assert res["collided"].sum() == 0
