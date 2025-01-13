# simulation.py
from geometry import calculate_nearest_boundary, calculate_void_si_max
from physics import elastic_scattering, sample_new_direction_cosines
import numpy as np

def simulate_single_particle(args):
    """
    Simulates a single particle and returns partial tally results.
    """
    state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng = args

    # Call the updated simulate_particle function
    result, absorbed, fissioned, _, final_energy, region_count, trajectory = simulate_particle(
        state, reader, mediums, A, N, sampler, region_bounds, track_coordinates=track_coordinates, rng=rng
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

def simulate_particle(state, reader, mediums, A, N, sampler, region_bounds=None, track_coordinates=False, rng= None):
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
    epsilon=1e-6
    region_count = 0
    absorbed_coordinates = []
    fission_coordinates = []
    trajectory = [] if track_coordinates else None
    x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]
    u = np.sin(state["theta"]) * np.cos(state["phi"])
    v = np.sin(state["theta"]) * np.sin(state["phi"])
    w = np.cos(state["theta"])
    while True:
        # Track coordinates if enabled

        if track_coordinates:
            trajectory.append((state["x"], state["y"], state["z"]))#You can add other stuff for tracking, ex: sigma_s, sigma_a, Energy, etc..

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

        nearest_point, nearest_medium, nearest_distance = calculate_nearest_boundary(state, mediums, u, v, w)

        if current_medium.is_void:
            state["x"], state["y"], state["z"] = nearest_point
            # Apply a small offset in the direction of travel
            state["x"] += epsilon * u
            state["y"] += epsilon * v
            state["z"] += epsilon * w
            continue

        # Get cross-sections for the current medium
        sigma_s, sigma_a, sigma_f, Sigma_t = reader.get_cross_sections(
            current_medium.element, state["energy"], sampler, N
        )

        # Calculate distance to next interaction
        si = -np.log(1 - rng.random()) / Sigma_t

        if si > nearest_distance:
            state["x"], state["y"], state["z"] = nearest_point
            state["x"] += epsilon * u
            state["y"] += epsilon * v
            state["z"] += epsilon * w
            continue

        state["x"] += si * u
        state["y"] += si * v
        state["z"] += si * w
        x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]

        # Check if particle has left the current medium
        if not current_medium.contains(state["x"], state["y"], state["z"]):
            continue

        # Determine type of interaction
        interaction_prob = rng.random()
        if interaction_prob < sigma_s / Sigma_t:  # Scattering
            state["has_interacted"] = True
            E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler, rng)
            state["theta"] = np.arccos(mu_lab)
            state["energy"] = E_prime
            u, v, w, state["phi"] = sample_new_direction_cosines(u, v, w, mu_lab, rng)

        elif interaction_prob < (sigma_s + sigma_a) / Sigma_t:  # Absorption
            state["has_interacted"] = True
            absorbed_coordinates.append((state["x"], state["y"], state["z"]))
            return "absorbed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count, trajectory
