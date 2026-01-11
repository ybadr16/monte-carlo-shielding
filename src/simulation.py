from .geometry import calculate_nearest_boundary, calculate_void_si_max
# Updated imports to include get_nuclear_temperature
from .physics import (
    elastic_scattering,
    inelastic_scattering,
    sample_new_direction_cosines,
    sample_maxwellian,
    get_nuclear_temperature
)
import numpy as np

def simulate_single_particle(args):
    """
    Simulates a single source particle and all its descendants (Bank/Stack method).
    """
    # Unpack settings
    initial_state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings = args

    # --- THE PARTICLE BANK ---
    bank = [initial_state.copy()]

    # --- TALLY ACCUMULATORS ---
    total_absorbed_weight = 0.0
    absorbed_events = []
    fission_events = 0
    escaped_weight = 0.0
    escaped_energy_sum = 0.0
    region_detections = 0
    full_trajectory = []

    # --- PROCESS BANK UNTIL EMPTY ---
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
        # FIX: Passing 'A' correctly here
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

                # --- TABULATED CONTINUUM / MULTIPLICATION ---
                # MT 91 (Continuum), 16 (n,2n), 17 (n,3n)
                if selected_mt in [16, 17, 91]:
                    # 1. Try Reader
                    E_out = reader.get_secondary_energy(current_medium.element, selected_mt, state["energy"], rng)

                    # 2. Fallback (Fermi Gas) if missing
                    if E_out is None:
                        E_avail = max(1.0, state["energy"] - selected_q)
                        T_nuc = get_nuclear_temperature(E_avail, A)
                        E_out = sample_maxwellian(T_nuc, rng)

                    if selected_mt == 91:
                        E_prime = E_out
                        mu_lab = 2 * rng.random() - 1

                    elif selected_mt == 16:
                        E1 = E_out
                        # For child, read again or use fallback
                        E2 = reader.get_secondary_energy(current_medium.element, 16, state["energy"], rng)
                        if E2 is None:
                            E_avail = max(1.0, state["energy"] - selected_q)
                            T_nuc = get_nuclear_temperature(E_avail, A)
                            E2 = sample_maxwellian(T_nuc, rng)

                        # Conserve
                        E_avail = max(1.0, state["energy"] - selected_q)
                        if (E1 + E2) > E_avail:
                            scale = E_avail / (E1 + E2)
                            E1 *= scale; E2 *= scale

                        E_prime = E1
                        mu_lab = 2 * rng.random() - 1

                        # Child
                        child = state.copy()
                        child["energy"] = E2
                        child["has_interacted"] = True
                        child["was_in_region"] = True
                        # Randomize child direction
                        u2,v2,w2 = 2*rng.random()-1, 2*rng.random()-1, 2*rng.random()-1
                        norm = np.sqrt(u2**2+v2**2+w2**2)
                        child["theta"] = np.arccos(w2/norm)
                        child["phi"] = np.arctan2(v2, u2)
                        children.append(child)

                    elif selected_mt == 17:
                        E1 = E_out
                        # Simplified for n,3n robustness
                        E2 = E1; E3 = E1
                        E_avail = max(1.0, state["energy"] - selected_q)
                        if (E1+E2+E3) > E_avail:
                            scale = E_avail/(E1+E2+E3)
                            E1*=scale; E2*=scale; E3*=scale

                        E_prime = E1
                        mu_lab = 2 * rng.random() - 1

                        for e_c in [E2, E3]:
                            child = state.copy()
                            child["energy"] = e_c
                            child["has_interacted"] = True
                            child["was_in_region"] = True
                            u2,v2,w2 = 2*rng.random()-1, 2*rng.random()-1, 2*rng.random()-1
                            norm = np.sqrt(u2**2+v2**2+w2**2)
                            child["theta"] = np.arccos(w2/norm)
                            child["phi"] = np.arctan2(v2, u2)
                            children.append(child)

                else:
                    # Discrete Levels (51-90)
                    E_prime, mu_cm, mu_lab = inelastic_scattering(state["energy"], A, selected_q, sampler, rng)

        else:
            # --- ELASTIC SCATTERING ---
            from .physics import get_vectors_and_vcm
            vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(state["energy"], A, sampler, rng)
            Vn_speed = np.linalg.norm(vec_Vn)

            mu_cm = reader.get_elastic_mu(current_medium.element, state["energy"], rng)

            phi_cm = 2 * np.pi * rng.random()
            sin_cm = np.sqrt(max(0.0, 1.0 - mu_cm**2))

            u_prime = sin_cm * np.cos(phi_cm)
            v_prime = sin_cm * np.sin(phi_cm)
            w_prime = mu_cm

            vec_Vn_prime = Vn_speed * np.array([u_prime, v_prime, w_prime])
            vec_vn_prime = vec_Vn_prime + vec_vcm
            vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
            E_prime = 0.5 * 1.674927471e-27 * vn_prime_speed_sq / 1.60217663e-19

            vn_prime_speed = np.sqrt(vn_prime_speed_sq)
            if vn_prime_speed > 0:
                mu_lab = vec_vn_prime[2] / vn_prime_speed
            else:
                mu_lab = 1.0

            E_prime = max(1e-5, E_prime)

        u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
        state["theta"] = np.arccos(w)
        state["energy"] = E_prime
