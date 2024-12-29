# simulation.py
from geometry import calculate_direction_cosines
from cross_sections import get_cross_sections
from physics import elastic_scattering, sample_new_direction_cosines
import random
import math


def simulate_particle(state, reader, mediums, A, sampler, region_bounds=None):
    """
    Simulate the trajectory of a single particle through multiple mediums.

    Args:
        reader: CrossSectionReader instance
        mediums: List of Medium instances defining the geometry.
        initial_energy: Initial neutron energy (eV).
        A: Mass number of the target.
        sampler: VelocitySampler instance for thermal effects.

    Returns:
        result: Final state of the particle ("escaped", "absorbed", "fission").
        absorbed_coordinates: List of absorption points.
        fission_coordinates: List of fission points.
        new_particles: List of newly generated particles (for fission).
        final_energy: Particle's final energy.
    """

    region_count = 0

    absorbed_coordinates = []
    fission_coordinates = []
    x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]

    while True:
        # Determine the current medium based on highest priority
        if region_bounds:
            x_min, x_max, y_min, y_max, z_min, z_max = region_bounds
            if x_min <= state["x"] <= x_max and y_min <= state["y"] <= y_max and z_min <= state["z"] <= z_max:
                region_count += 1


        current_medium = None
        max_priority = -float('inf')  # Start with lowest possible priority
        for medium in mediums:
            if medium.contains(state["x"], state["y"], state["z"]) and medium.priority > max_priority:
                current_medium = medium
                max_priority = medium.priority

        if current_medium is None:
            # Particle is outside all defined mediums, escapes
            return "escaped", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count

        # If the particle is in a void medium, propagate freely
        if current_medium.is_void:
            # Calculate the free propagation distance without interaction
            si = random.uniform(0.1, 10.0)  # Arbitrary step size for free streaming
            state["x"] += si * math.sin(state["theta"]) * math.cos(state["phi"])
            state["y"] += si * math.sin(state["theta"]) * math.sin(state["phi"])
            state["z"] += si * math.cos(state["theta"])

            print("particle in void")
            # Check if the particle has exited the void medium
            if not current_medium.contains(state["x"], state["y"], state["z"]):
                continue

            continue  # Continue propagating in the void
        print(f"particle in lead, {state['x']}, {state['y']}, {state['z']}")
        # Get cross-sections for the current medium
        sigma_s, sigma_a, sigma_f, Sigma_t = get_cross_sections(
            reader, current_medium.element, state["energy"], sampler
        )

        # Calculate distance to next interaction
        si = -math.log(1 - random.random()) / Sigma_t

        # Update position
        if not state["has_interacted"]:
            state["x"] += si * math.sin(state["theta"]) * math.cos(state["phi"])
            state["y"] += si * math.sin(state["theta"]) * math.sin(state["phi"])
            state["z"] += si * math.cos(state["theta"])
            u, v, w = calculate_direction_cosines(state["x"], state["y"], state["z"], x_prev, y_prev, z_prev)
        else:
            state["x"] += si * u
            state["y"] += si * v
            state["z"] += si * w
            x_prev, y_prev, z_prev = state["x"], state["y"], state["z"]

        # Check if particle has left the current medium
        if not current_medium.contains(state["x"], state["y"], state["z"]):
            continue

        # Determine type of interaction
        interaction_prob = random.random()
        if interaction_prob < sigma_s / Sigma_t:  # Scattering
            state["has_interacted"] = True
            E_prime, mu_cm, mu_lab = elastic_scattering(state["energy"], A, sampler)
            state["energy"] = E_prime
            u, v, w = sample_new_direction_cosines(u, v, w, mu_lab)

        elif interaction_prob < (sigma_s + sigma_a) / Sigma_t:  # Absorption
            absorbed_coordinates.append((state["x"], state["y"], state["z"]))
            return "absorbed", absorbed_coordinates, fission_coordinates, None, state["energy"], region_count
        elif interaction_prob < (sigma_s + sigma_a + sigma_f) / Sigma_t:  # Fission
            fission_coordinates.append((state["x"], state["y"], state["z"]))
            new_particles = [
                {"x": state["x"], "y": state["y"], "z": state["z"],
                 "theta": random.uniform(0, math.pi),
                 "phi": random.uniform(0, 2 * math.pi),
                 "energy": state["energy"]}
                for _ in range(2)
            ]
            return "fission", absorbed_coordinates, fission_coordinates, new_particles, state["energy"], region_count
