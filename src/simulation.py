from .geometry import calculate_nearest_boundary, calculate_void_si_max
from .physics import (
    elastic_scattering,
    inelastic_scattering,
    sample_new_direction_cosines,
    sample_maxwellian,
    get_nuclear_temperature,
    get_vectors_and_vcm,
    sample_forward_biased_mu
)
import numpy as np

# Constants needed for inline calc
m_n = 1.674927471e-27
eV_to_J = 1.60217663e-19

def simulate_single_particle(args):
    """
    Simulates a single source particle and all its descendants (Bank/Stack method).
    """
    initial_state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings = args

    bank = [initial_state.copy()]

    total_absorbed_weight = 0.0
    absorbed_events = []
    fission_events = 0
    escaped_weight = 0.0
    escaped_energy_sum = 0.0
    region_detections = 0
    full_trajectory = []

    while bank:
        current_state = bank.pop()

        status, child_particles, trajectory_segment, abs_weight, abs_coords = run_particle_kernel(
            current_state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings
        )

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
    }

def run_particle_kernel(state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings):
    epsilon = 1e-6
    trajectory = [] if track_coordinates else None
    absorbed_weight_local = 0.0
    absorbed_coords_local = []
    children = []

    u = np.sin(state["theta"]) * np.cos(state["phi"])
    v = np.sin(state["theta"]) * np.sin(state["phi"])
    w = np.cos(state["theta"])

    while True:
        # --- 1. TRACKING ---
        if track_coordinates:
            trajectory.append((state["x"], state["y"], state["z"]))

        # --- 2. REGION DETECTION ---
        if region_bounds:
            x_min, x_max, y_min, y_max, z_min, z_max = region_bounds
            if x_min <= state["x"] <= x_max and y_min <= state["y"] <= y_max and z_min <= state["z"] <= z_max:
                if not state.get("was_in_region", False):
                    state["region_detected"] = True
                    state["was_in_region"] = True

        # --- 3. LOCATE MEDIUM ---
        current_medium = None
        max_priority = -float('inf')
        point_check = (state["x"], state["y"], state["z"])
        for medium in mediums:
            if medium.contains(*point_check) and medium.priority > max_priority:
                current_medium = medium
                max_priority = medium.priority

        if current_medium is None:
            return "escaped", children, trajectory, absorbed_weight_local, absorbed_coords_local

        # --- 4. DISTANCE TO BOUNDARY ---
        nearest_point, _, nearest_distance = calculate_nearest_boundary(state, mediums, u, v, w)

        if nearest_point is None:
             return "escaped", children, trajectory, absorbed_weight_local, absorbed_coords_local

        # --- 5. VOID HANDLING ---
        if current_medium.is_void:
            if nearest_point is not None:
                state["x"], state["y"], state["z"] = nearest_point
                state["x"] += epsilon * u
                state["y"] += epsilon * v
                state["z"] += epsilon * w
            continue

        # --- 6. CROSS SECTIONS ---
        sigma_s, sigma_in, sigma_a, sigma_f, Sigma_t = reader.get_cross_sections(
            current_medium.element, state["energy"], sampler, N, A
        )
        total_scatter_xs = sigma_s + sigma_in

        # --- 7. SAMPLE DISTANCE ---
        if Sigma_t <= 0:
            si = float('inf')
        else:
            si = -np.log(1 - rng.random()) / Sigma_t

        # --- 8. MOVE ---
        if si > nearest_distance:
            state["x"], state["y"], state["z"] = nearest_point
            state["x"] += epsilon * u
            state["y"] += epsilon * v
            state["z"] += epsilon * w
            continue

        state["x"] += si * u
        state["y"] += si * v
        state["z"] += si * w

        if not current_medium.contains(state["x"], state["y"], state["z"]):
            continue

        state["has_interacted"] = True

        # --- 9. COLLISION KERNEL ---
        if settings and settings.use_implicit_capture:
            if state["weight"] < settings.weight_cutoff:
                if rng.random() < settings.roulette_survival_prob:
                    state["weight"] /= settings.roulette_survival_prob
                else:
                    return "killed", children, trajectory, absorbed_weight_local, absorbed_coords_local

            if Sigma_t > 0:
                p_absorb = (sigma_a + sigma_f) / Sigma_t
                p_scatter = total_scatter_xs / Sigma_t
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

        else:
            # Analog Monte Carlo logic
            interaction_prob = rng.random()
            p_scatter_el = sigma_s / Sigma_t if Sigma_t > 0 else 0
            p_scatter_in = sigma_in / Sigma_t if Sigma_t > 0 else 0

            if interaction_prob < p_scatter_el:
                pass
            elif interaction_prob < (p_scatter_el + p_scatter_in):
                pass
            else:
                absorbed_weight_local += state["weight"]
                absorbed_coords_local.append((state["x"], state["y"], state["z"], state["weight"]))
                return "absorbed", children, trajectory, absorbed_weight_local, absorbed_coords_local

        # === SCATTERING KINEMATICS ===
        if (sigma_s + sigma_in) > 0:
            prob_inelastic = sigma_in / (sigma_s + sigma_in)
        else:
            prob_inelastic = 0.0

        if rng.random() < prob_inelastic:
            # --- INELASTIC SCATTERING ---
            _, components = reader.get_inelastic_components(current_medium.element, state["energy"], N)

            if not components:
                E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler, rng)
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

                # --- TABULATED CONTINUUM / MULTIPLICATION (Law 61 & Law 4) ---
                if selected_mt in [16, 17, 91]:
                    E_in_snapshot = state["energy"]

                    # Calculate available energy budget
                    E_avail = max(1.0, E_in_snapshot - selected_q)

                    # ---------------------------------------------------------
                    # [PHYSICS FIX] PARENT NEUTRON
                    # ---------------------------------------------------------
                    E_prime = 0.0
                    mu_lab = 0.0
                    valid_sample = False

                    # Attempt Rejection Sampling to satisfy Energy Conservation
                    for _ in range(10):
                        # A. Try Correlated (Law 61) - e.g. Be9, Fe56
                        E_try, mu_try, success = reader.get_secondary_correlated_sample(
                            current_medium.element, selected_mt, E_in_snapshot, rng
                        )

                        # B. Try Uncorrelated (Law 4) - e.g. Pb208
                        # If Law 61 failed (success=False), check if Law 4 data exists
                        if not success:
                            E_uncorr = reader.get_secondary_energy(current_medium.element, selected_mt, E_in_snapshot, rng)
                            if E_uncorr is not None:
                                E_try = E_uncorr
                                mu_try = None # Law 4 usually implies isotropic or separate angular file
                                success = True

                        if success:
                            if E_try <= E_avail:
                                E_prime = E_try
                                mu_lab = mu_try if mu_try is not None else (2 * rng.random() - 1)
                                valid_sample = True
                                break

                    # C. Fallback: Evaporation Model
                    if not valid_sample:
                        T_nuc = get_nuclear_temperature(E_avail, A)
                        E_prime = sample_maxwellian(T_nuc, rng)
                        E_prime = min(E_prime, E_avail) # Hard clamp is acceptable in fallback
                        mu_lab = 2 * rng.random() - 1

                    # Update Parent State
                    u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
                    state["theta"] = np.arccos(w)
                    state["energy"] = max(1e-5, E_prime)

                    # ---------------------------------------------------------
                    # [PHYSICS FIX] CHILDREN NEUTRONS
                    # ---------------------------------------------------------

                    if selected_mt == 16:  # (n,2n) -> 1 Child
                        # Sample child independently from the same ENDF distribution as the
                        # parent — no E_c <= E_remaining constraint.  The original constraint
                        # caused systematic rejection when E_remaining was small (parent took
                        # a large share), forcing the Maxwellian fallback and a soft spectrum.
                        # OpenMC samples both neutrons from the distribution independently;
                        # per-event energy conservation is not enforced.
                        child_E, _, succ_c = reader.get_secondary_correlated_sample(
                            current_medium.element, selected_mt, E_in_snapshot, rng
                        )
                        if not succ_c:
                            child_E_uncorr = reader.get_secondary_energy(
                                current_medium.element, selected_mt, E_in_snapshot, rng
                            )
                            if child_E_uncorr is not None:
                                child_E = child_E_uncorr
                                succ_c = True
                        if not succ_c:
                            T_nuc = get_nuclear_temperature(E_avail, A)
                            child_E = sample_maxwellian(T_nuc, rng)

                        child = state.copy()
                        child["energy"] = max(1e-5, child_E)
                        child["has_interacted"] = True
                        child["was_in_region"] = True

                        # New direction
                        mu_c = 2 * rng.random() - 1
                        uc, vc, wc, child["phi"] = sample_new_direction_cosines(u, v, w, mu_c, rng)
                        child["theta"] = np.arccos(wc)
                        children.append(child)

                    elif selected_mt == 17:  # (n,3n) -> 2 Children
                        E_remaining = max(1e-5, E_avail - E_prime)

                        # Simple valid approximation: Split remaining budget 50/50
                        E_c1 = E_remaining * 0.5
                        E_c2 = E_remaining * 0.5

                        for e_c in [E_c1, E_c2]:
                            child = state.copy()
                            child["energy"] = max(1e-5, e_c)
                            child["has_interacted"] = True
                            child["was_in_region"] = True
                            mu_c = 2 * rng.random() - 1
                            uc, vc, wc, child["phi"] = sample_new_direction_cosines(u, v, w, mu_c, rng)
                            child["theta"] = np.arccos(wc)
                            children.append(child)

                    continue

                else:
                    # Discrete Levels (51-90)
                    E_prime, mu_cm, mu_lab = inelastic_scattering(state["energy"], A, selected_q, sampler, rng)
                    u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
                    state["theta"] = np.arccos(w)
                    state["energy"] = E_prime

        else:
            # --- ELASTIC SCATTERING ---
            mu_cm = reader.get_elastic_mu(current_medium.element, state["energy"], rng)

            A_sq = A**2
            term = A_sq + 1 + 2*A*mu_cm
            E_prime = state["energy"] * term / ((A + 1)**2)

            numerator = 1 + A * mu_cm
            denominator = np.sqrt(term)
            mu_lab = numerator / denominator

            E_prime = max(1e-5, E_prime)
            mu_lab = max(-1.0, min(1.0, mu_lab))

            u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
            state["theta"] = np.arccos(w)
            state["energy"] = E_prime
