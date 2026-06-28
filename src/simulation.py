from .geometry import (
    calculate_nearest_crossing,
    calculate_void_si_max,
    reflect_direction,
)
from .physics import (
    elastic_scattering,
    inelastic_scattering,
    sample_new_direction_cosines,
    sample_maxwellian,
    get_nuclear_temperature,
    sample_watt_spectrum,
    watt_params_for,
)
import numpy as np

m_n = 1.674927471e-27
eV_to_J = 1.60217663e-19

# Below this energy elastic scattering uses free-gas vector kinematics (target
# motion explicit, upscatter allowed); above it the static-target formula with
# the tabulated angular distribution. Scattering is s-wave below the cutoff.
THERMAL_KINEMATICS_CUTOFF = 10.0  # eV

def simulate_single_particle(args):
    """Track a source neutron and its banked descendants to termination."""
    initial_state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings = args

    bank = [initial_state.copy()]

    total_absorbed_weight = 0.0
    absorbed_events = []
    fission_events = 0
    fission_bank = []
    escaped_weight = 0.0
    escaped_energy_sum = 0.0
    region_detections = 0
    full_trajectory = []
    counters = {}
    source_uncollided = False
    is_source = True

    while bank:
        current_state = bank.pop()

        status, child_particles, trajectory_segment, abs_weight, abs_coords = run_particle_kernel(
            current_state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings,
            counters, fission_bank
        )

        if is_source:
            source_uncollided = (status == "escaped"
                                 and not current_state.get("has_interacted", False))
            is_source = False

        total_absorbed_weight += abs_weight
        absorbed_events.extend(abs_coords)

        if track_coordinates and trajectory_segment:
            full_trajectory.extend(trajectory_segment)

        if status == "escaped":
            escaped_weight += current_state["weight"]
            escaped_energy_sum += current_state["energy"] * current_state["weight"]

        elif status == "fission":
            fission_events += 1

        if current_state.get("region_detected", False):
            region_detections += 1

        if child_particles:
            bank.extend(child_particles)

    final_E = escaped_energy_sum / escaped_weight if escaped_weight > 0 else 0.0

    return {
        "result": "simulated",
        "absorbed_weight": total_absorbed_weight,
        "absorbed_coords": absorbed_events,
        "fissioned": fission_events,
        "final_energy": final_E,
        "final_weight": escaped_weight,
        "region_detected": region_detections > 0,
        "trajectory": full_trajectory if track_coordinates else None,
        "physics_counts": counters,
        "source_uncollided": source_uncollided,
        "fission_sites": fission_bank,
    }


def _inc(counters, key, n=1):
    if counters is not None:
        counters[key] = counters.get(key, 0) + n


class _Iso:
    """Single-isotope composition entry for the backward-compatible path.

    Mirrors the attribute interface of `material.Nuclide` so the kernel can
    treat the legacy `element` + global A/N/sampler exactly like one mixture
    component.
    """
    __slots__ = ("element", "number_density", "atomic_weight_ratio", "sampler")

    def __init__(self, element, number_density, atomic_weight_ratio, sampler):
        self.element = element
        self.number_density = number_density
        self.atomic_weight_ratio = atomic_weight_ratio
        self.sampler = sampler


def _select_isotope(iso_data, weight_fn, total, rng):
    """Pick one (iso, s_el, s_in, s_cap, s_fis) entry weighted by `weight_fn`.

    With a single isotope no random number is drawn, so single-material runs
    keep their exact RNG stream (and byte-identical results vs. the legacy path).
    """
    if len(iso_data) == 1:
        return iso_data[0]
    roll = rng.random() * total
    acc = 0.0
    for entry in iso_data:
        acc += weight_fn(entry)
        if roll <= acc:
            return entry
    return iso_data[-1]


def _cm_to_lab(E_cm, mu_cm, E_in, A):
    """CM-to-lab boost of an emitted neutron, mu_cm measured from the incident direction.

        E_lab  = E_cm + [E_in + 2 mu_cm (A+1) sqrt(E_in E_cm)] / (A+1)^2
        mu_lab = mu_cm sqrt(E_cm/E_lab) + sqrt(E_in/E_lab)/(A+1)
    """
    Ap1 = A + 1.0
    E_lab = E_cm + (E_in + 2.0 * mu_cm * Ap1 * np.sqrt(E_in * E_cm)) / (Ap1 * Ap1)
    if E_lab <= 0.0:
        return max(1e-5, E_cm), mu_cm
    mu_lab = mu_cm * np.sqrt(E_cm / E_lab) + np.sqrt(E_in / E_lab) / Ap1
    return E_lab, max(-1.0, min(1.0, mu_lab))


def _emit_secondary_neutron(reader, element, mt, E_in, E_avail, A, rng):
    """Sample (energy_eV, mu, source) for one emitted (n,xn) neutron.

    mu is the lab-frame cosine from the incident direction; source is
    'law61', 'law4' or 'maxwellian'.
    """
    E, mu, ok = reader.get_secondary_correlated_sample(element, mt, E_in, rng)
    src = "law61"
    if not ok:
        E_uncorr = reader.get_secondary_energy(element, mt, E_in, rng)
        if E_uncorr is not None:
            E, mu, ok, src = E_uncorr, None, True, "law4"
    if not ok:
        E = sample_maxwellian(get_nuclear_temperature(E_avail, A), rng)
        return max(1e-5, E), None, "maxwellian"

    # isotropic CM cosine if no angular law, then boost CM data to lab
    if mu is None:
        mu = 2.0 * rng.random() - 1.0
    if reader.get_secondary_frame(element, mt) == "cm":
        E, mu = _cm_to_lab(E, mu, E_in, A)
    return max(1e-5, E), mu, src


def _advance_across_boundary(state, nearest_point, nearest_surface, u, v, w, epsilon):
    """Move the neutron onto the boundary it is about to cross and apply that
    surface's boundary condition.

    For a "reflective" surface the direction is specularly reflected and the
    neutron is nudged back into the current medium; otherwise it transmits across
    (the default vacuum-leakage behaviour). Returns the (possibly reflected)
    direction cosines (u, v, w).
    """
    px, py, pz = nearest_point
    if (nearest_surface is not None
            and getattr(nearest_surface, "boundary_type", "transmission") == "reflective"):
        nx, ny, nz = nearest_surface.normal(px, py, pz)
        u, v, w = reflect_direction(u, v, w, nx, ny, nz)
        # Step off the surface ALONG THE NORMAL by epsilon, in the inward sense
        # (the side the reflected ray goes). Nudging along the reflected
        # direction instead fails at grazing incidence: there the reflected ray
        # is nearly parallel to the surface, the normal clearance is ~0, the
        # point stays within the on-surface tolerance, and the neutron
        # re-triggers a zero-distance crossing indefinitely (a reflection trap).
        s = epsilon if (u * nx + v * ny + w * nz) > 0 else -epsilon
        state["x"] = px + s * nx
        state["y"] = py + s * ny
        state["z"] = pz + s * nz
    else:
        # transmission: nudge past the surface along the direction of travel
        state["x"] = px + epsilon * u
        state["y"] = py + epsilon * v
        state["z"] = pz + epsilon * w
    return u, v, w


def run_particle_kernel(state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings,
                        counters=None, fission_bank=None):
    epsilon = 1e-6
    trajectory = [] if track_coordinates else None
    absorbed_weight_local = 0.0
    absorbed_coords_local = []
    children = []

    u = np.sin(state["theta"]) * np.cos(state["phi"])
    v = np.sin(state["theta"]) * np.sin(state["phi"])
    w = np.cos(state["theta"])

    while True:
        if track_coordinates:
            trajectory.append((state["x"], state["y"], state["z"]))

        # detector region crossing
        if region_bounds:
            x_min, x_max, y_min, y_max, z_min, z_max = region_bounds
            if x_min <= state["x"] <= x_max and y_min <= state["y"] <= y_max and z_min <= state["z"] <= z_max:
                if not state.get("was_in_region", False):
                    state["region_detected"] = True
                    state["was_in_region"] = True

        # highest-priority medium containing the current point
        current_medium = None
        max_priority = -float('inf')
        point_check = (state["x"], state["y"], state["z"])
        for medium in mediums:
            if medium.contains(*point_check) and medium.priority > max_priority:
                current_medium = medium
                max_priority = medium.priority

        if current_medium is None:
            return "escaped", children, trajectory, absorbed_weight_local, absorbed_coords_local

        nearest_point, _, nearest_distance, nearest_surface = calculate_nearest_crossing(
            state, mediums, u, v, w)

        if nearest_point is None:
             return "escaped", children, trajectory, absorbed_weight_local, absorbed_coords_local

        # stream through voids to the next surface
        if current_medium.is_void:
            if nearest_point is not None:
                u, v, w = _advance_across_boundary(
                    state, nearest_point, nearest_surface, u, v, w, epsilon)
            continue

        # per-isotope macroscopic cross sections for this medium; an explicit
        # mixture uses medium.composition, otherwise fall back to the single
        # isotope carried by the global element/A/N/sampler (legacy behaviour)
        composition = current_medium.composition or [
            _Iso(current_medium.element, N, A, sampler)
        ]

        iso_data = []
        Sig_el = Sig_in = Sig_cap = Sig_fis = 0.0
        for iso in composition:
            s_el, s_in, s_cap, s_fis, _ = reader.get_cross_sections(
                iso.element, state["energy"], iso.sampler,
                iso.number_density, iso.atomic_weight_ratio
            )
            iso_data.append((iso, s_el, s_in, s_cap, s_fis))
            Sig_el += s_el
            Sig_in += s_in
            Sig_cap += s_cap
            Sig_fis += s_fis

        Sig_scatter = Sig_el + Sig_in
        Sigma_t = Sig_scatter + Sig_cap + Sig_fis

        # distance to next collision
        if Sigma_t <= 0:
            si = float('inf')
        else:
            si = -np.log(1 - rng.random()) / Sigma_t

        if si > nearest_distance:
            u, v, w = _advance_across_boundary(
                state, nearest_point, nearest_surface, u, v, w, epsilon)
            continue

        state["x"] += si * u
        state["y"] += si * v
        state["z"] += si * w

        if not current_medium.contains(state["x"], state["y"], state["z"]):
            continue

        state["has_interacted"] = True

        # collision estimator of fission production (k-eff): every collision
        # contributes w·Σᵢ νᵢ·Σ_f,ᵢ / Σ_t, lower variance than counting actual
        # fission events. Only scored when a fission bank is active (criticality).
        if fission_bank is not None and Sig_fis > 0.0 and counters is not None:
            nu_sigf = 0.0
            for entry in iso_data:
                if entry[4] > 0.0:
                    nu_sigf += reader.get_nu(entry[0].element, state["energy"]) * entry[4]
            counters["fission_production"] = (
                counters.get("fission_production", 0.0)
                + state["weight"] * nu_sigf / Sigma_t)

        if settings and settings.use_implicit_capture:
            if state["weight"] < settings.weight_cutoff:
                if rng.random() < settings.roulette_survival_prob:
                    state["weight"] /= settings.roulette_survival_prob
                else:
                    return "killed", children, trajectory, absorbed_weight_local, absorbed_coords_local

            if Sigma_t > 0:
                p_absorb = (Sig_cap + Sig_fis) / Sigma_t
                p_scatter = Sig_scatter / Sigma_t
            else:
                p_absorb = 0.0
                p_scatter = 1.0

            w_loss = state["weight"] * p_absorb
            absorbed_weight_local += w_loss
            if w_loss > 0:
                absorbed_coords_local.append((state["x"], state["y"], state["z"], w_loss))

            state["weight"] *= p_scatter
            if state["weight"] <= 0:
                 return "killed", children, trajectory, absorbed_weight_local, absorbed_coords_local

            # survives as a scatter: pick which isotope it scatters off,
            # weighted by each isotope's scattering macroscopic cross section
            iso, s_el_i, s_in_i, _, _ = _select_isotope(
                iso_data, lambda e: e[1] + e[2], Sig_scatter, rng
            )

        else:
            # analog: pick the collision isotope (weighted by total XS), then
            # elastic / inelastic / capture within that isotope
            iso, s_el_i, s_in_i, s_cap_i, s_fis_i = _select_isotope(
                iso_data, lambda e: e[1] + e[2] + e[3] + e[4], Sigma_t, rng
            )
            s_tot_i = s_el_i + s_in_i + s_cap_i + s_fis_i

            interaction_prob = rng.random()
            p_scatter_el = s_el_i / s_tot_i if s_tot_i > 0 else 0
            p_scatter_in = s_in_i / s_tot_i if s_tot_i > 0 else 0

            if interaction_prob < p_scatter_el:
                pass
            elif interaction_prob < (p_scatter_el + p_scatter_in):
                pass
            else:
                # absorption. With a fission bank (criticality / k-eff driver),
                # split fission off capture and emit ν prompt neutrons into the
                # NEXT generation's source; otherwise fission is pure absorption.
                if (fission_bank is not None and s_fis_i > 0.0
                        and rng.random() < s_fis_i / (s_cap_i + s_fis_i)):
                    nu_bar = reader.get_nu(iso.element, state["energy"])
                    n_emit = int(nu_bar + rng.random())  # integer multiplicity
                    a_w, b_w = watt_params_for(iso.element)
                    for _ in range(n_emit):
                        fission_bank.append({
                            "x": state["x"], "y": state["y"], "z": state["z"],
                            "theta": np.arccos(2 * rng.random() - 1),
                            "phi": 2 * np.pi * rng.random(),
                            "energy": sample_watt_spectrum(rng, a_w, b_w),
                            "weight": 1.0, "has_interacted": False,
                        })
                    _inc(counters, "fission_events")
                    return "fission", children, trajectory, absorbed_weight_local, absorbed_coords_local
                absorbed_weight_local += state["weight"]
                absorbed_coords_local.append((state["x"], state["y"], state["z"], state["weight"]))
                return "absorbed", children, trajectory, absorbed_weight_local, absorbed_coords_local

        # the selected isotope drives every downstream lookup / kinematic
        element_i = iso.element
        A_i = iso.atomic_weight_ratio
        sampler_i = iso.sampler
        N_i = iso.number_density

        if (s_el_i + s_in_i) > 0:
            prob_inelastic = s_in_i / (s_el_i + s_in_i)
        else:
            prob_inelastic = 0.0

        if rng.random() < prob_inelastic:
            _, components = reader.get_inelastic_components(element_i, state["energy"], N_i)

            if not components:
                E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A_i, sampler_i, rng)
            else:
                xs_values = [c[1] for c in components]
                total_in_xs = sum(xs_values)
                roll = rng.random() * total_in_xs
                running_sum = 0.0
                selected_q = 0.0
                selected_mt = 0

                for mt, xs, q_val in components:
                    running_sum += xs
                    if roll <= running_sum:
                        selected_q = abs(q_val)
                        selected_mt = mt
                        break

                # continuum / multiplication channels (n,2n), (n,3n), (n,nc)
                if selected_mt in [16, 17, 91]:
                    E_in_snapshot = state["energy"]
                    # emitted cosines are measured from the incident direction
                    u_in, v_in, w_in = u, v, w
                    _inc(counters, "nxn_events")

                    E_avail = max(1.0, E_in_snapshot - selected_q)

                    # parent neutron: Law 61, then Law 4, then evaporation,
                    # rejecting energies above the available budget
                    E_prime = 0.0
                    mu_lab = 0.0
                    valid_sample = False
                    parent_source = None

                    for _ in range(10):
                        E_try, mu_try, success = reader.get_secondary_correlated_sample(
                            element_i, selected_mt, E_in_snapshot, rng
                        )
                        src = "law61"

                        if not success:
                            E_uncorr = reader.get_secondary_energy(element_i, selected_mt, E_in_snapshot, rng)
                            if E_uncorr is not None:
                                E_try = E_uncorr
                                mu_try = None
                                success = True
                                src = "law4"

                        if success:
                            if E_try <= E_avail:
                                E_prime = E_try
                                mu_lab = mu_try if mu_try is not None else (2 * rng.random() - 1)
                                valid_sample = True
                                parent_source = src
                                break

                    if not valid_sample:
                        T_nuc = get_nuclear_temperature(E_avail, A_i)
                        E_prime = sample_maxwellian(T_nuc, rng)
                        E_prime = min(E_prime, E_avail)
                        mu_lab = 2 * rng.random() - 1
                        parent_source = "maxwellian"
                    _inc(counters, f"parent_{parent_source}")

                    # the energy budget is checked in the native frame, so boost
                    # to the lab frame afterwards for CM-frame data
                    if parent_source != "maxwellian" and reader.get_secondary_frame(
                            element_i, selected_mt) == "cm":
                        E_prime, mu_lab = _cm_to_lab(E_prime, mu_lab, E_in_snapshot, A_i)

                    u, v, w, state["phi"] = sample_new_direction_cosines(u_in, v_in, w_in, mu_lab, rng)
                    state["theta"] = np.arccos(w)
                    state["energy"] = max(1e-5, E_prime)

                    # (n,2n) banks one extra neutron, (n,3n) two; each is sampled
                    # independently with no per-neutron energy constraint
                    n_children = {16: 1, 17: 2}.get(selected_mt, 0)
                    if selected_mt == 17:
                        _inc(counters, "n3n_events")
                    for _ in range(n_children):
                        child_E, child_mu, child_source = _emit_secondary_neutron(
                            reader, element_i, selected_mt,
                            E_in_snapshot, E_avail, A_i, rng)
                        _inc(counters, f"child_{child_source}")

                        child = state.copy()
                        child["energy"] = child_E
                        child["has_interacted"] = True
                        child["was_in_region"] = True
                        mu_c = child_mu if child_mu is not None else (2 * rng.random() - 1)
                        uc, vc, wc, child["phi"] = sample_new_direction_cosines(
                            u_in, v_in, w_in, mu_c, rng)
                        child["theta"] = np.arccos(wc)
                        children.append(child)

                    continue

                else:
                    # discrete levels, MT 51-90
                    _inc(counters, "discrete_inelastic_events")
                    E_prime, mu_cm, mu_lab = inelastic_scattering(state["energy"], A_i, selected_q, sampler_i, rng)
                    u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
                    state["theta"] = np.arccos(w)
                    state["energy"] = E_prime

        else:
            if state["energy"] < THERMAL_KINEMATICS_CUTOFF:
                # free-gas kinematics below the cutoff: the moving target lets the
                # neutron upscatter and equilibrate to the medium temperature
                E_prime, _, mu_lab = elastic_scattering(state["energy"], A_i, sampler_i, rng)
            else:
                # static target with the tabulated CM angular distribution
                mu_cm = reader.get_elastic_mu(element_i, state["energy"], rng)

                A_sq = A_i**2
                term = A_sq + 1 + 2*A_i*mu_cm
                E_prime = state["energy"] * term / ((A_i + 1)**2)

                numerator = 1 + A_i * mu_cm
                denominator = np.sqrt(term)
                mu_lab = numerator / denominator

                E_prime = max(1e-5, E_prime)
                mu_lab = max(-1.0, min(1.0, mu_lab))

            u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
            state["theta"] = np.arccos(w)
            state["energy"] = E_prime
