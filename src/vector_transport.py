"""Event-based (vectorized) transport engine.

Processes a whole bank of neutrons in lockstep with NumPy array operations —
one "event" (free flight, then boundary crossing or collision) per pass over
the live set — instead of one Python-interpreter iteration per collision per
history. The interpreter overhead that dominates the scalar kernel's runtime
is paid once per event for the whole bank, not once per neutron.

The scalar kernel in simulation.py is unchanged and remains both the readable
reference implementation and the correctness oracle: the array-form geometry
and cross-section lookups here must match it exactly (they are deterministic),
while full transport is compared statistically, because the per-history order
of RNG consumption necessarily differs between the two engines.

Current scope (milestone 1): geometry, cross sections, free flights, and
boundary crossings with transmission/specular reflection. A collision
terminates the history and is scored as "collided", which is exactly the
estimator a Beer-Lambert uncollided-transmission measurement needs. Collision
physics (elastic, free-gas thermal, inelastic, (n,xn)) arrives in the next
milestones.
"""
import numpy as np

from .angular_distribution import AngularDistribution
from .geometry import get_primitive_surfaces
from .medium import Plane, Sphere, Cylinder
from .simulation import THERMAL_KINEMATICS_CUTOFF

EPSILON = 1e-6  # boundary-crossing nudge, matches the scalar kernel


# ---------------------------------------------------------------------------
# Array-form surface distances.
#
# Each function reproduces the corresponding scalar nearest_surface_method
# branch for branch, with "no forward intersection" (scalar None) encoded as
# np.inf. Keeping the arithmetic term-for-term identical to medium.py is what
# lets the unit tests demand agreement to ~1 ulp.
# ---------------------------------------------------------------------------

def _plane_distances(s, x, y, z, u, v, w):
    num = -s.D - s.A * x - s.B * y - s.C * z
    den = s.A * u + s.B * v + s.C * w
    d = np.full(x.shape, np.inf)
    on_surface = np.abs(num) <= 1e-8
    d = np.where(on_surface & (np.abs(den) > 1e-8), 0.0, d)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = num / den
    forward = (~on_surface) & (den != 0) & (t >= 0)
    return np.where(forward, t, d)


def _sphere_distances(s, x, y, z, u, v, w):
    dn = np.sqrt(u * u + v * v + w * w)
    ok = dn != 0
    safe_dn = np.where(ok, dn, 1.0)
    un, vn, wn = u / safe_dn, v / safe_dn, w / safe_dn
    xb, yb, zb = x - s.x0, y - s.y0, z - s.z0
    k = xb * un + yb * vn + zb * wn
    c = xb * xb + yb * yb + zb * zb - s.radius ** 2
    disc = k * k - c
    hit = ok & (disc >= 0)
    sq = np.sqrt(np.where(disc >= 0, disc, 0.0))
    d1 = -k - sq
    d2 = -k + sq
    inside = c < 0
    outside_d = np.where(d1 >= 0, d1, np.where(d2 >= 0, d2, np.inf))
    d = np.where(inside, np.maximum(d1, d2), outside_d)
    return np.where(hit, d, np.inf)


def _cylinder_distances(s, x, y, z, u, v, w):
    if s.axis == "z":
        rb1, rb2 = x - s.x0, y - s.y0
        a = u * u + v * v
        k = rb1 * u + rb2 * v
        du_, dv_ = u, v
    elif s.axis == "x":
        rb1, rb2 = y - s.y0, z - s.z0
        a = v * v + w * w
        k = rb1 * v + rb2 * w
        du_, dv_ = v, w
    else:  # "y"
        rb1, rb2 = x - s.x0, z - s.z0
        a = u * u + w * w
        k = rb1 * u + rb2 * w
        du_, dv_ = u, w
    c = rb1 * rb1 + rb2 * rb2 - s.radius ** 2

    valid = a != 0
    disc = k * k - a * c
    hit = valid & (disc >= 0)
    sq = np.sqrt(np.where(disc >= 0, disc, 0.0))
    safe_a = np.where(valid, a, 1.0)
    d1 = (-k - sq) / safe_a
    d2 = (-k + sq) / safe_a

    inside = c < 0
    on_surface = c == 0
    dot = rb1 * du_ + rb2 * dv_
    on_d = np.where((dot < 0) & (d2 >= 0), d2, np.inf)
    outside_d = np.where(d1 >= 0, d1, np.where(d2 >= 0, d2, np.inf))
    d = np.where(inside, np.maximum(d1, d2), np.where(on_surface, on_d, outside_d))
    return np.where(hit, d, np.inf)


def surface_distances(surface, x, y, z, u, v, w):
    """Forward distances from (x, y, z) along (u, v, w) to `surface`, np.inf
    where the scalar method would return None."""
    if isinstance(surface, Plane):
        return _plane_distances(surface, x, y, z, u, v, w)
    if isinstance(surface, Sphere):
        return _sphere_distances(surface, x, y, z, u, v, w)
    if isinstance(surface, Cylinder):
        return _cylinder_distances(surface, x, y, z, u, v, w)
    raise TypeError(f"Unsupported surface type: {type(surface).__name__}")


def surface_normals(surface, x, y, z):
    """Unit surface normals at the given points, array-form of `normal`."""
    if isinstance(surface, Plane):
        shape = np.shape(x)
        return (np.full(shape, surface.A), np.full(shape, surface.B),
                np.full(shape, surface.C))
    if isinstance(surface, Sphere):
        dx, dy, dz = x - surface.x0, y - surface.y0, z - surface.z0
        mag = np.sqrt(dx * dx + dy * dy + dz * dz)
        mag = np.where(mag == 0, 1.0, mag)
        return dx / mag, dy / mag, dz / mag
    if isinstance(surface, Cylinder):
        zero = np.zeros(np.shape(x))
        if surface.axis == "z":
            dx, dy = x - surface.x0, y - surface.y0
            mag = np.where((m := np.sqrt(dx * dx + dy * dy)) == 0, 1.0, m)
            return dx / mag, dy / mag, zero
        if surface.axis == "x":
            dy, dz = y - surface.y0, z - surface.z0
            mag = np.where((m := np.sqrt(dy * dy + dz * dz)) == 0, 1.0, m)
            return zero, dy / mag, dz / mag
        dx, dz = x - surface.x0, z - surface.z0
        mag = np.where((m := np.sqrt(dx * dx + dz * dz)) == 0, 1.0, m)
        return dx / mag, zero, dz / mag
    raise TypeError(f"Unsupported surface type: {type(surface).__name__}")


# ---------------------------------------------------------------------------
# Array-form CSG containment and boundary search.
# ---------------------------------------------------------------------------

def contains_many(region, x, y, z, tolerance=1e-9):
    """Boolean mask of points inside `region`; array form of Region.contains."""
    from .medium import Region

    op = region.operation
    surfaces = region.surfaces

    def _inside(s):
        if isinstance(s, Region):
            return contains_many(s, x, y, z, tolerance)
        return s.evaluate(x, y, z) <= tolerance

    if op == "intersection":
        res = np.ones(np.shape(x), dtype=bool)
        for s in surfaces:
            res &= _inside(s)
        return res
    if op == "union":
        res = np.zeros(np.shape(x), dtype=bool)
        for s in surfaces:
            res |= _inside(s)
        return res
    if op == "complement":
        return ~_inside(surfaces[0])
    if op == "difference":
        if len(surfaces) != 2:
            raise ValueError("Difference operation requires exactly two surfaces")
        a, b = surfaces
        return _inside(a) & ~_inside(b)
    raise ValueError(f"Unknown operation: {op}")


def _region_primitives(region):
    prims = getattr(region, "_prim_cache", None)
    if prims is None:
        prims = get_primitive_surfaces(region.surfaces)
        try:
            region._prim_cache = prims
        except AttributeError:
            pass
    return prims


def nearest_crossings(regions, x, y, z, u, v, w):
    """Nearest boundary crossing for every particle at once.

    Mirrors geometry.calculate_nearest_crossing: for each region's primitive
    surfaces, a candidate distance counts only if the crossing point lies in
    the region (the CSG check), and the strictly smallest distance wins — the
    same region/surface iteration order as the scalar loop, so tie-breaking is
    identical. Returns (distances, surface_ids, surfaces) where surface_ids
    index into the flattened `surfaces` list (-1 = no crossing / escape).
    """
    n = x.size
    best = np.full(n, np.inf)
    best_sid = np.full(n, -1, dtype=np.int64)
    surfaces = []

    sid = 0
    for region in regions:
        for surface in _region_primitives(region):
            d = surface_distances(surface, x, y, z, u, v, w)
            better = np.isfinite(d) & (d >= 0) & (d < best)
            idx = np.nonzero(better)[0]
            if idx.size:
                px = x[idx] + d[idx] * u[idx]
                py = y[idx] + d[idx] * v[idx]
                pz = z[idx] + d[idx] * w[idx]
                inside = contains_many(region, px, py, pz)
                sel = idx[inside]
                best[sel] = d[sel]
                best_sid[sel] = sid
            surfaces.append(surface)
            sid += 1

    return best, best_sid, surfaces


def resolve_media(mediums, x, y, z):
    """Index of the highest-priority medium containing each point (-1 = none),
    array form of the kernel's current-medium scan (strict >, so the earlier
    medium wins a priority tie, as in the scalar loop)."""
    n = x.size
    midx = np.full(n, -1, dtype=np.int64)
    best_pri = np.full(n, -np.inf)
    for i, m in enumerate(mediums):
        inside = contains_many(m, x, y, z)
        upd = inside & (m.priority > best_pri)
        midx[upd] = i
        best_pri[upd] = m.priority
    return midx


# ---------------------------------------------------------------------------
# Vectorized macroscopic cross sections.
# ---------------------------------------------------------------------------

class VectorXS:
    """Macroscopic cross sections of one medium for an energy array.

    Uses the reader's per-nuclide fast tables with np.interp, which implements
    exactly the clamped lin-lin rule of the scalar _locate/_read pair — so the
    lookup is bit-compatible with reader.get_cross_sections above the thermal
    Doppler cutoff. The sub-10-eV elastic flight-energy shift and the URR
    band sampling are per-particle random processes and are deferred to the
    collision-physics milestones.
    """

    def __init__(self, reader, composition):
        self.nuclides = []
        for iso in composition:
            if iso.element not in reader._fast_tables:
                reader._build_fast_table(iso.element)
            tbl = reader._fast_tables[iso.element]
            if tbl is None:
                raise ValueError(f"No cross-section data for {iso.element}")
            self.nuclides.append((iso, tbl))

    def macroscopic(self, energy):
        E = np.asarray(energy, dtype=float)
        Sel = np.zeros_like(E)
        Sin = np.zeros_like(E)
        Scap = np.zeros_like(E)
        Sfis = np.zeros_like(E)
        for iso, tbl in self.nuclides:
            scale = 1e-24 * iso.number_density
            g = tbl["grid"]
            Sel += np.interp(E, g, tbl["el"]) * scale
            Sin += np.interp(E, g, tbl["in"]) * scale
            Scap += np.interp(E, g, tbl["cap"]) * scale
            Sfis += np.interp(E, g, tbl["fis"]) * scale
        return Sel, Sin, Scap, Sfis, Sel + Sin + Scap + Sfis

    def per_nuclide(self, energy):
        """Per-nuclide (elastic, inelastic, capture, fission) macroscopic
        arrays, each of shape (n_nuclides, n_particles) — the vector analogue
        of the kernel's iso_data list, used for per-collision isotope
        selection."""
        E = np.asarray(energy, dtype=float)
        el, inl, cap, fis = [], [], [], []
        for iso, tbl in self.nuclides:
            scale = 1e-24 * iso.number_density
            g = tbl["grid"]
            el.append(np.interp(E, g, tbl["el"]) * scale)
            inl.append(np.interp(E, g, tbl["in"]) * scale)
            cap.append(np.interp(E, g, tbl["cap"]) * scale)
            fis.append(np.interp(E, g, tbl["fis"]) * scale)
        return np.stack(el), np.stack(inl), np.stack(cap), np.stack(fis)


# ---------------------------------------------------------------------------
# Vectorized collision physics.
# ---------------------------------------------------------------------------

_ELEMENT_MAP = {"C": "C0", "Graphite": "C0", "Be": "Be9",
                "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}


def _get_angular(reader, element):
    """Reader-cached elastic angular distribution, as in get_elastic_mu."""
    actual = _ELEMENT_MAP.get(element, element)
    if actual not in reader.angular_dists:
        dist = AngularDistribution(reader.base_path, actual, mt=2)
        dist.load()
        reader.angular_dists[actual] = dist
    return reader.angular_dists[actual]


def sample_mu_many(dist, energy, rng):
    """CM elastic cosines for an energy array, vector form of
    AngularDistribution.sample_mu.

    Same statistical interpolation: each particle picks one bracketing
    incident-energy block with the linear weight, then inverts that block's
    tabulated CDF. Particles are grouped by selected block so each distinct
    block costs one np.interp over its members rather than one call per
    particle."""
    E = np.asarray(energy, dtype=float)
    m = E.size
    if not dist.loaded:
        return 2.0 * rng.random(m) - 1.0

    grid = dist.energy_grid
    idx = np.clip(np.searchsorted(grid, E) - 1, 0, len(grid) - 2)
    f = np.clip((E - grid[idx]) / (grid[idx + 1] - grid[idx]), 0.0, 1.0)
    block = np.where(rng.random(m) > f, idx, idx + 1)

    r = rng.random(m)
    mu = np.empty(m)
    for b in np.unique(block):
        sel = block == b
        start = dist.offsets[b]
        end = (dist.mu_data.shape[1] if b == len(grid) - 1
               else dist.offsets[b + 1])
        mu[sel] = np.interp(r[sel], dist.mu_data[2, start:end],
                            dist.mu_data[0, start:end])
    return mu


def static_elastic(E, A, mu_cm):
    """Static-target two-body elastic kinematics for arrays, mirroring the
    kernel's E >= 10 eV branch: outgoing energy and lab cosine from the CM
    cosine and the atomic weight ratio."""
    term = A * A + 1.0 + 2.0 * A * mu_cm
    E_prime = np.maximum(1e-5, E * term / ((A + 1.0) ** 2))
    mu_lab = np.clip((1.0 + A * mu_cm) / np.sqrt(term), -1.0, 1.0)
    return E_prime, mu_lab


def rotate_directions(u, v, w, mu_lab, rng):
    """Rotate unit directions by lab cosine mu_lab with a uniform azimuth —
    array form of physics.sample_new_direction_cosines, including the
    near-pole branch and the final renormalization."""
    m = u.size
    phi = 2.0 * np.pi * rng.random(m)
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - mu_lab * mu_lab))

    denom = np.sqrt(np.maximum(1e-12, 1.0 - w * w))
    ratio = sin_theta / denom
    u_new = mu_lab * u + ratio * (u * w * cphi - v * sphi)
    v_new = mu_lab * v + ratio * (v * w * cphi + u * sphi)
    w_new = mu_lab * w - sin_theta * denom * cphi
    norm = np.sqrt(u_new * u_new + v_new * v_new + w_new * w_new)
    u_new, v_new, w_new = u_new / norm, v_new / norm, w_new / norm

    pole = np.abs(w) >= 0.999999
    if pole.any():
        sign = np.where(w > 0, 1.0, -1.0)
        u_new = np.where(pole, sin_theta * cphi, u_new)
        v_new = np.where(pole, sin_theta * sphi, v_new)
        w_new = np.where(pole, sign * mu_lab, w_new)
    return u_new, v_new, w_new


# ---------------------------------------------------------------------------
# Event loop.
# ---------------------------------------------------------------------------

def _apply_crossings(bank, gi, px, py, pz, cu, cv, cw, sids, surfaces):
    """Apply boundary conditions to crossers: transmission nudges past along
    the flight direction; a reflective surface bounces the direction and steps
    off along the inward normal (the grazing-incidence reflection-trap fix,
    mirrored from simulation._advance_across_boundary). Writes positions and
    directions back into the full-size bank arrays at global indices `gi`."""
    x, y, z, u, v, w = bank
    for s_id in np.unique(sids):
        surf = surfaces[s_id]
        sm = sids == s_id
        g = gi[sm]
        if getattr(surf, "boundary_type", "transmission") == "reflective":
            nx, ny, nz = surface_normals(surf, px[sm], py[sm], pz[sm])
            dot = cu[sm] * nx + cv[sm] * ny + cw[sm] * nz
            ru = cu[sm] - 2.0 * dot * nx
            rv = cv[sm] - 2.0 * dot * ny
            rw = cw[sm] - 2.0 * dot * nz
            step = np.where(ru * nx + rv * ny + rw * nz > 0, EPSILON, -EPSILON)
            x[g] = px[sm] + step * nx
            y[g] = py[sm] + step * ny
            z[g] = pz[sm] + step * nz
            u[g], v[g], w[g] = ru, rv, rw
        else:
            x[g] = px[sm] + EPSILON * cu[sm]
            y[g] = py[sm] + EPSILON * cv[sm]
            z[g] = pz[sm] + EPSILON * cw[sm]

class _LegacyIso:
    """Adapter giving a bare `element=` medium the Nuclide attribute interface,
    the vector-engine counterpart of simulation._Iso."""
    __slots__ = ("element", "number_density", "atomic_weight_ratio", "sampler")

    def __init__(self, element, number_density, atomic_weight_ratio, sampler):
        self.element = element
        self.number_density = number_density
        self.atomic_weight_ratio = atomic_weight_ratio
        self.sampler = sampler


def _medium_xs_tables(reader, mediums, legacy=None):
    """One VectorXS per non-void medium (None for voids). `legacy` is an
    optional (A, N, sampler) tuple backing mediums that carry only `element`."""
    tables = []
    for m in mediums:
        if m.is_void:
            tables.append(None)
            continue
        comp = m.composition
        if comp is None:
            if legacy is None:
                raise ValueError(
                    f"Medium {m.name!r} has no composition and no legacy "
                    "(A, N, sampler) fallback was given")
            A, N, sampler = legacy
            comp = [_LegacyIso(m.element, N, A, sampler)]
        tables.append(VectorXS(reader, comp))
    return tables


def run_streaming(reader, mediums, source, rng, legacy=None, max_events=100_000):
    """Vectorized transport in which a collision terminates the history.

    source: dict of equal-length float arrays x, y, z, u, v, w, energy (unit
    direction cosines), optionally weight. rng: a numpy Generator. Returns a
    dict with per-particle fate masks and the leakage fraction — the
    uncollided-transmission estimator (Beer-Lambert observable).
    """
    x = np.array(source["x"], dtype=float)
    y = np.array(source["y"], dtype=float)
    z = np.array(source["z"], dtype=float)
    u = np.array(source["u"], dtype=float)
    v = np.array(source["v"], dtype=float)
    w = np.array(source["w"], dtype=float)
    E = np.array(source["energy"], dtype=float)
    n = x.size
    weight = np.array(source.get("weight", np.ones(n)), dtype=float)

    xs_tables = _medium_xs_tables(reader, mediums, legacy)

    alive = np.ones(n, dtype=bool)
    escaped = np.zeros(n, dtype=bool)
    collided = np.zeros(n, dtype=bool)
    events = 0

    while alive.any():
        events += 1
        if events > max_events:
            raise RuntimeError(f"exceeded max_events={max_events} with "
                               f"{int(alive.sum())} particles still alive")
        idx = np.nonzero(alive)[0]
        xi, yi, zi = x[idx], y[idx], z[idx]
        ui, vi, wi = u[idx], v[idx], w[idx]

        midx = resolve_media(mediums, xi, yi, zi)
        out = midx < 0
        if out.any():
            gone = idx[out]
            escaped[gone] = True
            alive[gone] = False
            keep = ~out
            idx = idx[keep]
            if idx.size == 0:
                continue
            xi, yi, zi = xi[keep], yi[keep], zi[keep]
            ui, vi, wi = ui[keep], vi[keep], wi[keep]
            midx = midx[keep]

        dist, sid, surfaces = nearest_crossings(mediums, xi, yi, zi, ui, vi, wi)
        no_cross = ~np.isfinite(dist)
        if no_cross.any():
            gone = idx[no_cross]
            escaped[gone] = True
            alive[gone] = False

        # flight length: per-medium total cross section, then s = -ln(1-xi)/Sigma_t
        Sigma_t = np.zeros(idx.size)
        for m_i in np.unique(midx):
            tbl = xs_tables[m_i]
            if tbl is None:
                continue  # void: Sigma_t stays 0 -> si = inf -> streams
            in_m = midx == m_i
            Sigma_t[in_m] = tbl.macroscopic(E[idx][in_m])[4]

        si = np.full(idx.size, np.inf)
        absorbing = Sigma_t > 0
        si[absorbing] = -np.log(1.0 - rng.random(int(absorbing.sum()))) \
            / Sigma_t[absorbing]

        live = ~no_cross
        crosses = live & (si > dist)
        collides = live & ~crosses

        # collisions terminate the history in this milestone
        if collides.any():
            hit = idx[collides]
            collided[hit] = True
            alive[hit] = False

        # boundary crossings: transmission nudges across, reflective surfaces
        # bounce and step off along the inward normal (the reflection-trap fix,
        # mirrored from simulation._advance_across_boundary)
        if crosses.any():
            c = np.nonzero(crosses)[0]
            _apply_crossings((x, y, z, u, v, w), idx[c],
                             xi[c] + dist[c] * ui[c],
                             yi[c] + dist[c] * vi[c],
                             zi[c] + dist[c] * wi[c],
                             ui[c], vi[c], wi[c], sid[c], surfaces)

    leaked = float(weight[escaped].sum())
    return {
        "escaped": escaped,
        "collided": collided,
        "leakage": leaked / float(weight.sum()),
        "events": events,
    }


def run_transport(reader, mediums, source, rng, settings, legacy=None,
                  max_events=1_000_000):
    """Vectorized fixed-source transport with survival biasing.

    Physics scope (milestone 3): implicit capture with weight-cutoff Russian
    roulette, per-collision isotope selection weighted by the scattering cross
    section, and static-target elastic scattering with the tabulated CM
    angular distribution — the kernel's E >= 10 eV elastic branch. Two lanes
    are transported approximately and counted so the caller can verify they
    never fired: a lane that selects an inelastic channel scatters
    elastically instead (counters["inelastic_as_elastic"]; milestone 4), and
    a lane below the thermal cutoff stays on the static path instead of
    free-gas (counters["sub10ev_static"]; milestone 4).

    Returns per-source-history leaked weight and escape energy arrays (feed
    them to statistics.escape_statistics for batch-means uncertainties), the
    total absorbed weight, and the counters dict.
    """
    x = np.array(source["x"], dtype=float)
    y = np.array(source["y"], dtype=float)
    z = np.array(source["z"], dtype=float)
    u = np.array(source["u"], dtype=float)
    v = np.array(source["v"], dtype=float)
    w = np.array(source["w"], dtype=float)
    E = np.array(source["energy"], dtype=float)
    n = x.size
    wgt = np.array(source.get("weight", np.ones(n)), dtype=float)

    if not settings.use_implicit_capture:
        raise NotImplementedError(
            "run_transport currently implements survival-biased (shielding) "
            "mode only; analog/criticality arrives with the fission milestone")

    xs_tables = _medium_xs_tables(reader, mediums, legacy)

    leaked_w = np.zeros(n)
    leaked_E = np.zeros(n)
    absorbed = 0.0
    counters = {"collisions": 0, "inelastic_as_elastic": 0,
                "sub10ev_static": 0, "killed": 0}
    alive = np.ones(n, dtype=bool)
    events = 0

    while alive.any():
        events += 1
        if events > max_events:
            raise RuntimeError(f"exceeded max_events={max_events} with "
                               f"{int(alive.sum())} particles still alive")
        idx = np.nonzero(alive)[0]
        xi, yi, zi = x[idx], y[idx], z[idx]
        ui, vi, wi = u[idx], v[idx], w[idx]

        midx = resolve_media(mediums, xi, yi, zi)
        out = midx < 0
        if out.any():
            g = idx[out]
            leaked_w[g] = wgt[g]
            leaked_E[g] = E[g]
            alive[g] = False
            keep = ~out
            idx = idx[keep]
            if idx.size == 0:
                continue
            xi, yi, zi = xi[keep], yi[keep], zi[keep]
            ui, vi, wi = ui[keep], vi[keep], wi[keep]
            midx = midx[keep]

        dist, sid, surfaces = nearest_crossings(mediums, xi, yi, zi, ui, vi, wi)
        no_cross = ~np.isfinite(dist)
        if no_cross.any():
            g = idx[no_cross]
            leaked_w[g] = wgt[g]
            leaked_E[g] = E[g]
            alive[g] = False

        Sigma_t = np.zeros(idx.size)
        for m_i in np.unique(midx):
            tbl = xs_tables[m_i]
            if tbl is None:
                continue  # void streams
            in_m = midx == m_i
            Sigma_t[in_m] = tbl.macroscopic(E[idx][in_m])[4]

        si = np.full(idx.size, np.inf)
        interacting = Sigma_t > 0
        si[interacting] = -np.log(1.0 - rng.random(int(interacting.sum()))) \
            / Sigma_t[interacting]

        live = ~no_cross
        crosses = live & (si > dist)
        collides = live & ~crosses

        if crosses.any():
            c = np.nonzero(crosses)[0]
            _apply_crossings((x, y, z, u, v, w), idx[c],
                             xi[c] + dist[c] * ui[c],
                             yi[c] + dist[c] * vi[c],
                             zi[c] + dist[c] * wi[c],
                             ui[c], vi[c], wi[c], sid[c], surfaces)

        if not collides.any():
            continue

        # advance colliders to the collision site
        c = np.nonzero(collides)[0]
        g_all = idx[c]
        x[g_all] = xi[c] + si[c] * ui[c]
        y[g_all] = yi[c] + si[c] * vi[c]
        z[g_all] = zi[c] + si[c] * wi[c]

        for m_i in np.unique(midx[c]):
            lanes = c[midx[c] == m_i]
            g = idx[lanes]
            # a collision site that left the medium is discarded, as in the
            # scalar kernel ("if not current_medium.contains: continue")
            inside = contains_many(mediums[m_i], x[g], y[g], z[g])
            g = g[inside]
            if g.size == 0:
                continue
            tbl = xs_tables[m_i]
            counters["collisions"] += int(g.size)

            # weight-cutoff Russian roulette, before the implicit absorption
            low = wgt[g] < settings.weight_cutoff
            if low.any():
                p = settings.roulette_survival_prob
                roll = rng.random(int(low.sum()))
                lw = g[low]
                die = lw[roll >= p]
                boost = lw[roll < p]
                wgt[boost] /= p
                alive[die] = False
                counters["killed"] += int(die.size)
                g = g[~low | np.isin(g, boost)]
                if g.size == 0:
                    continue

            Ew = E[g]
            el, inl, cap, fis = tbl.per_nuclide(Ew)
            Ssc = el.sum(axis=0) + inl.sum(axis=0)
            St = Ssc + cap.sum(axis=0) + fis.sum(axis=0)

            # implicit capture: score the absorbed weight, keep the scatter
            p_scatter = np.where(St > 0, Ssc / St, 1.0)
            absorbed += float((wgt[g] * (1.0 - p_scatter)).sum())
            wgt[g] = wgt[g] * p_scatter
            dead = wgt[g] <= 0
            if dead.any():
                alive[g[dead]] = False
                counters["killed"] += int(dead.sum())
                g = g[~dead]
                if g.size == 0:
                    continue
                el, inl = el[:, ~dead], inl[:, ~dead]
                Ew = E[g]

            # isotope selection weighted by the per-nuclide scattering XS;
            # single-isotope media pick directly, as in the scalar kernel
            k = len(tbl.nuclides)
            m = g.size
            if k == 1:
                choice = np.zeros(m, dtype=np.int64)
            else:
                wts = el + inl
                cum = np.cumsum(wts, axis=0)
                roll = rng.random(m) * cum[-1]
                choice = np.minimum((cum < roll).sum(axis=0), k - 1)

            lane_idx = np.arange(m)
            s_el = el[choice, lane_idx]
            s_in = inl[choice, lane_idx]
            denom = s_el + s_in
            p_inel = np.where(denom > 0, s_in / denom, 0.0)
            inelastic = rng.random(m) < p_inel
            counters["inelastic_as_elastic"] += int(inelastic.sum())
            counters["sub10ev_static"] += int(
                (Ew < THERMAL_KINEMATICS_CUTOFF).sum())

            # elastic scatter: tabulated CM cosine per selected nuclide
            mu_cm = np.empty(m)
            A_arr = np.empty(m)
            for ki in np.unique(choice):
                sel = choice == ki
                iso = tbl.nuclides[ki][0]
                dist_obj = _get_angular(reader, iso.element)
                mu_cm[sel] = sample_mu_many(dist_obj, Ew[sel], rng)
                A_arr[sel] = iso.atomic_weight_ratio

            E_prime, mu_lab = static_elastic(Ew, A_arr, mu_cm)
            un, vn, wn = rotate_directions(u[g], v[g], w[g], mu_lab, rng)
            E[g] = E_prime
            u[g], v[g], w[g] = un, vn, wn

    return {
        "leaked_weight": leaked_w,
        "escape_energy": leaked_E,
        "leakage": float(leaked_w.sum()) / float(np.sum(
            np.array(source.get("weight", np.ones(n)), dtype=float))),
        "absorbed_weight": absorbed,
        "counters": counters,
        "events": events,
    }
