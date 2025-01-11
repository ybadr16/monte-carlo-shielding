# simulation.py
from geometry import calculate_direction_cosines, calculate_nearest_boundary_distance
from physics import elastic_scattering, sample_new_direction_cosines
import math
import numpy as np

def simulate_single_particle(args):
    """
    Simulates a single particle and returns partial tally results.
    """
    state, reader, mediums, A, N, sampler, region_bounds, track_coordinates = args

    # Call the updated simulate_particle function
    result, absorbed, fissioned, _, final_energy, region_count, trajectory = simulate_particle(
        state, reader, mediums, A, N, sampler, region_bounds, track_coordinates=track_coordinates
    )

    # Return results including trajectory if tracking was enabled
    return {
        "result": result,
        "absorbed": absorbed,
        "fissioned": fissioned,
        "final_energy": final_energy,
        "region_detected": region_count > 0,
        "trajectory": trajectory if track_coordinates else None,
    }

def simulate_particle(state, reader, mediums, A, N, sampler, region_bounds=None, track_coordinates=False):
    """
    Simulate the trajectory of a single particle through multiple mediums.

    Args:
        state: Dictionary representing the particle's initial state.
        reader: CrossSectionReader instance for cross-section data.
        mediums: List of Medium instances defining the geometry.
        A: Mass number of the target.
        sampler: VelocitySampler instance for thermal effects.
        region_bounds: Optional boundaries for a specific region.
        track_coordinates: If True, records all particle coordinates during the simulation.

    Returns:
        result: Final state of the particle ("escaped", "absorbed", "fission").
        absorbed_coordinates: List of absorption points.
        fission_coordinates: List of fission points.
        new_particles: List of newly generated particles (for fission).
        final_energy: Particle's final energy.
        region_count: Number of times the particle was detected in the region.
        trajectory: List of all coordinates if track_coordinates is True, otherwise None.
    """
    region_count = 0
    absorbed_coordinates = []
    fission_coordinates = []
    trajectory = [] if track_coordinates else None
    x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]
    #sigma_s, sigma_a, sigma_f, Sigma_t = 0, 0, 0, 0
    while True:
        # Track coordinates if enabled

        if track_coordinates:
            trajectory.append((state["x"], state["y"], state["z"])) #, sigma_s, sigma_a, sigma_f, Sigma_t


        # Determine the current medium based on highest priority
        if region_bounds:
            x_min, x_max, y_min, y_max, z_min, z_max = region_bounds
            if x_min <= state["x"] <= x_max and y_min <= state["y"] <= y_max and z_min <= state["z"] <= z_max:
                region_count += 1
                state["within_region"] = True

        current_medium = None
        max_priority = -float('inf')
        for medium in mediums:
            if medium.contains(state["x"], state["y"], state["z"]) and medium.priority > max_priority:
                current_medium = medium
                max_priority = medium.priority

        if current_medium is None:
            # Particle is outside all defined mediums, escapes
            return "escaped", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory

        if current_medium.is_void:
            # Calculate the free propagation distance without interaction
            #nearest_distance = calculate_nearest_boundary_distance(state, mediums[0]) #be careful of this approach, its not dynamic enough to accomodate multiple mediums. Only works for this case.
            si = np.random.uniform(0.1, 10.0)  # Arbitrary step size for free streaming
            #si = nearest_distance
            state["x"] += si * math.sin(state["theta"]) * math.cos(state["phi"])
            state["y"] += si * math.sin(state["theta"]) * math.sin(state["phi"])
            state["z"] += si * math.cos(state["theta"])
            u, v, w = calculate_direction_cosines(state["x"], state["y"], state["z"], x_prev, y_prev, z_prev)
            state["has_interacted"] = True
            #print(state["x"], state["y"], state["z"])
            continue

        # Get cross-sections for the current medium
        sigma_s, sigma_a, sigma_f, Sigma_t = reader.get_cross_sections(
            current_medium.element, state["energy"], sampler, N
        )

        # Calculate distance to next interaction
        si = -math.log(1 - np.random.rand()) / Sigma_t

        # Update position
        if not state["has_interacted"] and not current_medium.is_void:
            state["x"] += si * math.sin(state["theta"]) * math.cos(state["phi"])
            state["y"] += si * math.sin(state["theta"]) * math.sin(state["phi"])
            state["z"] += si * math.cos(state["theta"])
            u, v, w = calculate_direction_cosines(state["x"], state["y"], state["z"], x_prev, y_prev, z_prev)
            state["has_interacted"] = True
        else:
            state["x"] += si * u
            state["y"] += si * v
            state["z"] += si * w
            x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]

        # Check if particle has left the current medium
        if not current_medium.contains(state["x"], state["y"], state["z"]):
            continue

        # Determine type of interaction
        interaction_prob = np.random.rand()
        if interaction_prob < sigma_s / Sigma_t:  # Scattering
            #print(f"Theta before scattering: {state['theta']}")
            state["has_interacted"] = True
            E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler)
            state["energy"] = E_prime
            u, v, w, state["phi"], state["theta"] = sample_new_direction_cosines(u, v, w, mu_lab)
            #print(f"Theta after scattering: {state['theta']}")


        elif interaction_prob < (sigma_s + sigma_a) / Sigma_t:  # Absorption
            state["has_interacted"] = True
            absorbed_coordinates.append((state["x"], state["y"], state["z"]))
            return "absorbed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory
        elif interaction_prob < (sigma_s + sigma_a + sigma_f) / Sigma_t:  # Fission
            state["has_interacted"] = True
            fission_coordinates.append((state["x"], state["y"], state["z"]))
            new_particles = [
                {"x": state["x"], "y": state["y"], "z": state["z"],
                 "theta": np.random.uniform(0, math.pi),
                 "phi": np.random.uniform(0, 2 * math.pi),
                 "energy": state["energy"]}
                for _ in range(2)
            ]
            return "fission", absorbed_coordinates, fission_coordinates, new_particles, state["energy"], region_count, trajectory
