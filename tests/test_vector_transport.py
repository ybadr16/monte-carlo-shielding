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


# ---------------------------------------------------------------------------
# Milestone 3: collision physics, validated statistically against the scalar
# kernel (same physics, different RNG consumption order).
# ---------------------------------------------------------------------------

def test_rotate_directions_units_and_mean():
    rng = np.random.default_rng(5)
    n = 20000
    mu = rng.uniform(-1, 1, n)
    _, _, _, u, v, w = _rays(n, seed=6)
    un, vn, wn = vt.rotate_directions(u, v, w, mu, rng)
    np.testing.assert_allclose(un**2 + vn**2 + wn**2, 1.0, atol=1e-9)
    # the rotated direction's cosine to the incident one must equal mu
    cosang = u * un + v * vn + w * wn
    np.testing.assert_allclose(cosang, mu, atol=1e-9)


def test_static_elastic_energy_bounds():
    rng = np.random.default_rng(8)
    n = 10000
    A = 55.454
    E = np.full(n, 1.0e6)
    mu = rng.uniform(-1, 1, n)
    E_prime, mu_lab = vt.static_elastic(E, np.full(n, A), mu)
    alpha = ((A - 1) / (A + 1)) ** 2
    assert (E_prime >= alpha * E - 1e-6).all() and (E_prime <= E + 1e-6).all()
    assert (np.abs(mu_lab) <= 1.0).all()


def _scalar_reference(reader, nuc, radius, energy, n, seed):
    """Per-history leaked weight / escape energy from the scalar kernel."""
    from src.random_number_generator import RNGHandler
    from src.settings import Settings
    from src.simulation import simulate_single_particle

    region = Region([Sphere(center=(0, 0, 0), radius=radius)],
                    name="sphere", element=nuc.element)
    settings = Settings("shielding", n)
    src_rng = np.random.default_rng(seed)
    leaked_w = np.zeros(n)
    leaked_E = np.zeros(n)
    for i in range(n):
        mu = 2 * src_rng.random() - 1
        state = {"x": 0.0, "y": 0.0, "z": 0.0,
                 "theta": float(np.arccos(mu)),
                 "phi": float(2 * np.pi * src_rng.random()),
                 "energy": energy, "weight": 1.0, "has_interacted": False}
        res = simulate_single_particle((
            state, reader, [region], nuc.atomic_weight_ratio,
            nuc.number_density, nuc.sampler, None, False,
            RNGHandler(seed + 17 * i), settings))
        leaked_w[i] = res["final_weight"]
        leaked_E[i] = res["final_energy"]
    return leaked_w, leaked_E


def _vector_run(reader, nuc, radius, energy, n, seed):
    region = Region([Sphere(center=(0, 0, 0), radius=radius)],
                    name="sphere", composition=[nuc])
    rng = np.random.default_rng(seed)
    mu = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - mu**2)
    source = {"x": np.zeros(n), "y": np.zeros(n), "z": np.zeros(n),
              "u": s * np.cos(phi), "v": s * np.sin(phi), "w": mu,
              "energy": np.full(n, energy)}
    from src.settings import Settings
    return vt.run_transport(reader, [region], source, rng,
                            Settings("shielding", n))


def _z(a, sa, b, sb):
    return abs(a - b) / np.sqrt(sa**2 + sb**2)


@needs_data
@pytest.mark.parametrize("element,density,amass,awr,radius,energy", [
    ("Pb208", 11.35, 208.0, 206.19, 10.0, 2.0e6),
    ("Fe56", 7.874, 55.845, 55.454, 15.0, 2.5e4),
])
def test_transport_matches_scalar_kernel(element, density, amass, awr,
                                         radius, energy):
    """Elastic-only fast/intermediate spheres: the vector engine must agree
    with the scalar kernel on leakage and mean escape energy within
    statistics, and its approximate lanes must never fire."""
    from src.statistics import escape_statistics

    if not os.path.exists(os.path.join(ENDF, "neutron", f"{element}.h5")):
        pytest.skip(f"{element} data not present")
    reader = CrossSectionReader(ENDF)
    mat = Material(element, density, amass, awr)
    nuc = Nuclide(element, mat.number_density, awr)

    n_s, n_v = 3000, 30000
    lw_s, le_s = _scalar_reference(reader, nuc, radius, energy, n_s, seed=101)
    stats_s = escape_statistics(lw_s, le_s, n_s)

    res = _vector_run(reader, nuc, radius, energy, n_v, seed=202)
    # elastic-only cases: no inelastic, (n,xn) or thermal free-gas lanes fire
    assert res["counters"]["nxn_events"] == 0
    assert res["counters"]["discrete_inelastic_events"] == 0
    assert res["counters"]["freegas_events"] == 0
    stats_v = escape_statistics(res["leaked_weight"], res["escape_energy"], n_v)

    z_leak = _z(stats_v["leakage"], stats_v["leakage_sem"],
                stats_s["leakage"], stats_s["leakage_sem"])
    z_energy = _z(stats_v["avg_energy"], stats_v["avg_energy_sem"],
                  stats_s["avg_energy"], stats_s["avg_energy_sem"])
    assert z_leak < 4.0, (stats_v, stats_s)
    assert z_energy < 4.0, (stats_v, stats_s)


# ---------------------------------------------------------------------------
# Milestone 4: free-gas thermal, discrete inelastic, (n,xn) multiplication.
# ---------------------------------------------------------------------------

@needs_data
@pytest.mark.parametrize("element,density,amass,awr,radius,energy,expect", [
    # thermal graphite sphere: exercises free-gas kinematics + Doppler lookup
    ("C12", 1.7, 12.011, 11.8969, 5.0, 0.0253, "freegas_events"),
    # 1 MeV iron sphere: discrete-level inelastic open (first level 847 keV)
    ("Fe56", 7.874, 55.845, 55.454, 10.0, 1.0e6, "discrete_inelastic_events"),
    # 14 MeV lead sphere: (n,2n) multiplication, leakage > 1, banked children
    ("Pb208", 11.35, 208.0, 206.19, 10.0, 1.4e7, "nxn_events"),
])
def test_transport_full_physics_matches_scalar(element, density, amass, awr,
                                               radius, energy, expect):
    from src.statistics import escape_statistics

    if not os.path.exists(os.path.join(ENDF, "neutron", f"{element}.h5")):
        pytest.skip(f"{element} data not present")
    reader = CrossSectionReader(ENDF)
    mat = Material(element, density, amass, awr)
    nuc = Nuclide(element, mat.number_density, awr)

    n_s, n_v = 2000, 20000
    lw_s, le_s = _scalar_reference(reader, nuc, radius, energy, n_s, seed=303)
    stats_s = escape_statistics(lw_s, le_s, n_s)

    res = _vector_run(reader, nuc, radius, energy, n_v, seed=404)
    assert res["counters"][expect] > 0, res["counters"]
    stats_v = escape_statistics(res["leaked_weight"], res["escape_energy"], n_v)

    z_leak = _z(stats_v["leakage"], stats_v["leakage_sem"],
                stats_s["leakage"], stats_s["leakage_sem"])
    z_energy = _z(stats_v["avg_energy"], stats_v["avg_energy_sem"],
                  stats_s["avg_energy"], stats_s["avg_energy_sem"])
    assert z_leak < 4.0, (element, stats_v, stats_s, res["counters"])
    assert z_energy < 4.0, (element, stats_v, stats_s, res["counters"])


# ---------------------------------------------------------------------------
# Milestone 5: URR band sampling, Watt spectrum, analog k-eigenvalue.
# ---------------------------------------------------------------------------

def test_sample_watt_many_moments():
    from src.physics import WATT_PARAMS
    a, b = WATT_PARAMS["U235"]
    rng = np.random.default_rng(21)
    E = vt.sample_watt_many(rng, 200_000, a, b)
    assert (E >= 0).all()
    mean_expected = 1.5 * a + a * a * b / 4.0   # ~2.03 MeV
    sem = E.std() / np.sqrt(E.size)
    assert abs(E.mean() - mean_expected) < 4 * sem


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ENDF, "neutron", "U238.h5")),
    reason="U238 data not present")
def test_urr_band_sampling_mean_neutral():
    """The probability-weighted URR band average must reproduce the smooth
    (infinite-dilution) cross section — the treatment is unbiased outside the
    self-shielding effect (paper Section on URR)."""
    reader = CrossSectionReader(ENDF)
    nuc = Nuclide("U238", 4.8e22, 236.006)
    vx = vt.VectorXS(reader, [nuc])
    n = 200_000
    E = np.full(n, 5.0e4)  # 50 keV, inside the U238 URR (20-149 keV)
    smooth_el, _, smooth_cap, smooth_fis = (a[0] for a in vx.per_nuclide(E[:1]))

    el, inl, cap, fis = vx.per_nuclide(E)
    el, cap, fis = vt._urr_adjust(reader, vx, el, cap, fis, E,
                                  np.random.default_rng(3))
    # band samples must actually vary (self-shielding is active)...
    assert el[0].std() > 0
    # ...but their mean reproduces the smooth value to a few permille
    assert abs(el[0].mean() / smooth_el[0] - 1.0) < 5e-3
    assert abs(cap[0].mean() / smooth_cap[0] - 1.0) < 5e-3


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ENDF, "neutron", "U235.h5")),
    reason="U235 data not present")
def test_keff_vector_matches_scalar():
    """Bare U235 metal sphere near the critical radius: the vector power
    iteration must agree with the scalar run_keff within statistics."""
    from src.criticality import run_keff, uniform_sphere_source
    from src.random_number_generator import RNGHandler
    from src.settings import Settings

    reader = CrossSectionReader(ENDF)
    rho, amass, awr, radius = 18.74, 235.0, 233.0248, 8.7
    mat = Material("U235", rho, amass, awr)
    nuc = Nuclide("U235", mat.number_density, awr, atomic_mass=amass)

    # scalar reference
    region_s = Region([Sphere((0, 0, 0), radius)], element="U235", priority=1)
    settings = Settings("criticality", 1)
    src_rng = RNGHandler(77)
    init = uniform_sphere_source((0, 0, 0), radius, 300, src_rng)
    ref = run_keff(reader, [region_s], init, awr, nuc.number_density,
                   nuc.sampler, settings, generations=24, inactive=6, seed=5)

    # vector run (larger bank, same physics)
    region_v = Region([Sphere((0, 0, 0), radius)], composition=[nuc],
                      priority=1)
    n = 2000
    rng = np.random.default_rng(88)
    r3 = radius * rng.random(n) ** (1.0 / 3.0)
    mu = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - mu ** 2)
    mu2 = rng.uniform(-1, 1, n)
    phi2 = rng.uniform(0, 2 * np.pi, n)
    s2 = np.sqrt(1 - mu2 ** 2)
    src = {"x": r3 * s * np.cos(phi), "y": r3 * s * np.sin(phi),
           "z": r3 * mu,
           "u": s2 * np.cos(phi2), "v": s2 * np.sin(phi2), "w": mu2,
           "energy": vt.sample_watt_many(rng, n, 0.988e6, 2.249e-6)}
    res = vt.run_keff_vector(reader, [region_v], src, settings,
                             generations=30, inactive=8, seed=9)

    z = _z(res["k_eff"], res["k_sem"], ref["k_eff"], ref["k_sem"])
    assert z < 4.0, (res["k_eff"], res["k_sem"], ref["k_eff"], ref["k_sem"])
    # sanity: near-critical sphere, k in a physically sensible window
    assert 0.95 < res["k_eff"] < 1.10
