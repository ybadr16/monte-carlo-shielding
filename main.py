# main.py
from collections import deque
from cross_section_read import CrossSectionReader
from vt_calc import VelocitySampler
from simulation import simulate_single_particle
from material import Material
from medium import Region, Plane, Cylinder
from multiprocessing import Pool
from tally import Tally
from random_number_generator import RNGHandler
import json
import time
import numpy as np

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

    mediums = [
        Region(
            surfaces=[
                Cylinder(radius=10, x0=0, y0=0, axis="z"),
                Plane(0, 0, -1, 10),  # z >= -10 (D=-10, no absolute)
                Plane(0, 0, 1, 10)     # z <= 10
            ],
            name="Cylinder",
            priority=1,
            element="Pb208"
        ),
        Region(
            surfaces=[
                Plane(-1, 0, 0, 20),  # x >= -20
                Plane(1, 0, 0, 20),    # x <= 20
                Plane(0, -1, 0, 40),  # y >= -40
                Plane(0, 1, 0, 40),    # y <= 40
                Plane(0, 0, -1, 20),  # z >= -20
                Plane(0, 0, 1, 20)     # z <= 20
            ],
            name="Void",
            priority=0,
            is_void=True
        )
    ]



    # Simulate particles
    num_particles = 100

    rngs = [RNGHandler(seed=12345 + i) for i in range(num_particles)]

    # Initialize Tally object
    tally = Tally()

    # Initialize particle queue
    particle_states = [
        {
            "x": -10.5, "y": 0.0, "z": 0.0,
            "theta": rng.uniform(0, np.pi),
            "phi": rng.uniform(0, 2 * np.pi),
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
