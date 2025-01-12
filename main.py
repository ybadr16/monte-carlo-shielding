# main.py
from collections import deque
from cross_section_read import CrossSectionReader
from vt_calc import VelocitySampler
from simulation import simulate_single_particle
from material import Material
from medium import Medium
from multiprocessing import Pool
from tally import Tally
from random_number_generator import RNGHandler
import math
import json
import time

def main():
    # Initialize cross-section reader and sampler
    base_path = "./endfb"
    reader = CrossSectionReader(base_path)

    lead = Material(name="Lead", density=11.35, atomic_mass=208, atomic_weight_ratio=2.5)
    cadmium = Material(name="Cadmium", density=8.65, atomic_mass=112, atomic_weight_ratio=2.3) #example for defining another element

    N = lead.number_density
    A = lead.atomic_weight_ratio
    mass_in_kg = lead.kg_mass
    sampler = VelocitySampler(mass=mass_in_kg)  # mass for Pb

    # Define mediums
    mediums = [
        Medium(name="Lead Shield", x_bounds=(-10, 10), y_bounds=(-20, 20), z_bounds=(-10, 10), element="Pb208", priority=1),
        Medium(name="Void", x_bounds=(-20, 20), y_bounds=(-40, 40), z_bounds=(-20, 20), is_void=True, priority=0)
    ]


    # Simulate particles
    num_particles = 100_000

    rngs = [RNGHandler(seed=12345 + i) for i in range(num_particles)]

    # Initialize Tally object
    tally = Tally()

    # Initialize particle queue
    particle_states = [
        {
            "x": -10.5, "y": 0.0, "z": 0.0,
            "theta": rng.uniform(0, math.pi),
            "phi": rng.uniform(0, 2 * math.pi),
            "has_interacted": False,
            "energy": 1e6,  # eV
        }
        for rng in rngs
    ]

    # Dictionary to store trajectories for all particles
    all_trajectories = {}
    # For tracking if a particle existed within a certain region (x_min, x_max, y_min, y_max, z_min, z_max)
    region_bounds = (14.9, 15.1, -15, 15, -15, 15)

    track = False
    # Prepare arguments for multiprocessing
    args = [
        (state, reader, mediums, A, N, sampler, region_bounds, track, rng)
        for state, rng in zip(particle_states, rngs)
    ]

    # Use multiprocessing to simulate particles
    sim_start_time = time.perf_counter()
    with Pool() as pool:
        partial_results = pool.map(simulate_single_particle, args)
    sim_end_time = time.perf_counter()

    # Merge partial results into the main tally
    for idx, result in enumerate(partial_results):
        tally.merge_partial_results(result)

        # Append trajectory to the dictionary
        if result["trajectory"]:
            all_trajectories[idx + 1] = result["trajectory"]

    # Output all trajectories
    with open("all_trajectories.json", "w") as f:
        json.dump(all_trajectories, f, indent=4)

    # Output results
    tally.print_summary(num_particles)

    print(f"Particle simulation time: {sim_end_time - sim_start_time:.2f} seconds")
    print(f"Average time per particle: {(sim_end_time - sim_start_time) / num_particles:.4f} seconds")

if __name__ == "__main__":
    main()
