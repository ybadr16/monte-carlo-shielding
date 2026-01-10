# src/simulation.py
from .geometry import calculate_nearest_boundary, calculate_void_si_max
# Import the new inelastic function
from .physics import elastic_scattering, inelastic_scattering, sample_new_direction_cosines
import numpy as np

def simulate_single_particle(args):
    """
    Simulates a single particle and returns partial tally results.
    """
    # Unpack settings
    state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings = args

    # Call the simulation kernel
    result, absorbed_coords, fissioned, _, final_energy, region_count, trajectory, total_absorbed_weight = simulate_particle(
        state, reader, mediums, A, N, sampler, region_bounds, track_coordinates=track_coordinates, rng=rng, settings=settings
    )

    # Return results
    return {
        "result": result,
        "absorbed_weight": total_absorbed_weight,
        "absorbed_coords": absorbed_coords,
        "fissioned": fissioned,
        "final_energy": final_energy,
        "final_weight": state["weight"],
        "region_detected": region_count > 0,
        "trajectory": trajectory if track_coordinates else None,
    }

def simulate_particle(state, reader, mediums, A, N, sampler, region_bounds=None, track_coordinates=False, rng=None, settings=None):
    """
    Simulate the trajectory of a single particle.
    """
    epsilon = 1e-6
    region_count = 0
    absorbed_coordinates = []
    fission_coordinates = []
    trajectory = [] if track_coordinates else None
    total_absorbed_weight = 0.0

    # Geometry State Initialization
    x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]
    u = np.sin(state["theta"]) * np.cos(state["phi"])
    v = np.sin(state["theta"]) * np.sin(state["phi"])
    w = np.cos(state["theta"])

    while True:
        # 1. TRACKING
        if track_coordinates:
            trajectory.append((state["x"], state["y"], state["z"]))

        # 2. REGION DETECTION
        if region_bounds:
            x_min, x_max, y_min, y_max, z_min, z_max = region_bounds
            if x_min <= state["x"] <= x_max and y_min <= state["y"] <= y_max and z_min <= state["z"] <= z_max:
                if not state.get("was_in_region", False):
                    region_count += 1
                    state["was_in_region"] = True

        # 3. LOCATE CURRENT MEDIUM
        current_medium = None
        max_priority = -float('inf')
        point_check = (state["x"], state["y"], state["z"])
        for medium in mediums:
            if medium.contains(*point_check) and medium.priority > max_priority:
                current_medium = medium
                max_priority = medium.priority

        if current_medium is None:
            return "escaped", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

        # 4. NEAREST BOUNDARY
        nearest_point, nearest_medium, nearest_distance = calculate_nearest_boundary(state, mediums, u, v, w)

        if nearest_point is None:
            return "escaped", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

        # 5. VOID HANDLING
        if current_medium.is_void:
            if nearest_point is not None:
                state["x"], state["y"], state["z"] = nearest_point
                state["x"] += epsilon * u
                state["y"] += epsilon * v
                state["z"] += epsilon * w
            continue

        # 6. GET CROSS SECTIONS (Updated for 5 returns)
        # Sigma_in is the SUM of all inelastic levels (MT 51-91)
        sigma_s, sigma_in, sigma_a, sigma_f, Sigma_t = reader.get_cross_sections(
            current_medium.element, state["energy"], sampler, N
        )

        # Combine Elastic (sigma_s) and Inelastic (sigma_in) for total scattering
        total_scatter_xs = sigma_s + sigma_in

        # 7. SAMPLE DISTANCE
        if Sigma_t <= 0:
            si = float('inf')
        else:
            si = -np.log(1 - rng.random()) / Sigma_t

        # 8. MOVE PARTICLE
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

        # --- PHYSICS KERNEL ---
        state["has_interacted"] = True

        # === IMPLICIT CAPTURE (Variance Reduction) ===
        if settings and settings.use_implicit_capture:
            # Roulette
            if state["weight"] < settings.weight_cutoff:
                if rng.random() < settings.roulette_survival_prob:
                    state["weight"] /= settings.roulette_survival_prob
                else:
                    return "killed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

            # Weight Reduction
            if Sigma_t > 0:
                p_absorb = (sigma_a + sigma_f) / Sigma_t
                p_scatter = total_scatter_xs / Sigma_t # Elastic + Inelastic
            else:
                p_absorb = 0.0
                p_scatter = 1.0

            weight_loss = state["weight"] * p_absorb
            total_absorbed_weight += weight_loss

            if weight_loss > 0:
                absorbed_coordinates.append((state["x"], state["y"], state["z"], weight_loss))

            state["weight"] *= p_scatter

            if state["weight"] <= 0:
                 return "killed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

        else:
            # === ANALOG MONTE CARLO ===
            interaction_prob = rng.random()
            if Sigma_t > 0:
                p_scatter_el = sigma_s / Sigma_t
                p_scatter_in = sigma_in / Sigma_t
                p_absorb = (sigma_a + sigma_f) / Sigma_t
            else:
                return "escaped", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

            # Sample Reaction Type
            current_prob = 0.0

            # Check Elastic
            current_prob += p_scatter_el
            if interaction_prob < current_prob:
                # ELASTIC SCATTERING LOGIC
                pass # Fall through to scattering section below

            # Check Inelastic
            elif interaction_prob < (current_prob + p_scatter_in):
                # INELASTIC SCATTERING LOGIC
                # We need to flag this to handle it differently below
                pass

            # Check Absorption
            else:
                weight_deposited = state["weight"]
                total_absorbed_weight += weight_deposited
                absorbed_coordinates.append((state["x"], state["y"], state["z"], weight_deposited))
                return "absorbed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory, total_absorbed_weight

        # === SCATTERING KINEMATICS ===

        # Decide between Elastic and Inelastic based on their relative probabilities
        # relative_prob_inelastic = sigma_in / (sigma_s + sigma_in)

        if (sigma_s + sigma_in) > 0:
            prob_inelastic = sigma_in / (sigma_s + sigma_in)
        else:
            prob_inelastic = 0.0

        if rng.random() < prob_inelastic:
            # --- INELASTIC SCATTERING ---
            # 1. Ask the reader: "Which elevator floor did we hit?"
            # This returns a list of (MT, XS, Q_value) tuples
            _, components = reader.get_inelastic_components(current_medium.element, state["energy"], N)

            if not components:
                # Fallback to elastic if something goes wrong
                E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler, rng)
            else:
                # 2. Sample which level (MT) occurred
                # We build a localized CDF for these levels
                xs_values = [c[1] for c in components]
                total_in_xs = sum(xs_values)

                roll = rng.random() * total_in_xs
                running_sum = 0.0
                selected_q = 0.0

                for mt, xs, q_val in components:
                    running_sum += xs
                    if roll <= running_sum:
                        selected_q = abs(q_val) # Q is usually negative in files, we need positive magnitude
                        break

                # 3. Perform Inelastic Physics
                # (Uses the new function we added to physics.py)
                E_prime, mu_cm, mu_lab = inelastic_scattering(state["energy"], A, selected_q, sampler, rng)

        else:
            # --- ELASTIC SCATTERING ---
            E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler, rng)

        # Update Direction and Energy
        u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)
        state["theta"] = np.arccos(w)
        state["energy"] = E_prime
