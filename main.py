# main.py
from collections import deque
from cross_section_read import CrossSectionReader
from vt_calc import VelocitySampler
from simulation import simulate_single_particle
from material import calculate_number_density
from medium import Medium
from multiprocessing import Pool
from tally import Tally  # Import the Tally class
import numpy as np
import math
import json




def main():
    # Initialize cross-section reader and sampler
    base_path = "./endfb"
    reader = CrossSectionReader(base_path)
    sampler = VelocitySampler(mass=3.44e-25, temperature=294)  # mass for Pb
    #sampler = VelocitySampler(mass=1.864e-25, temperature=294) #mass for cd

    # Define mediums
    mediums = [
        Medium(name="Lead Shield", x_bounds=(-10, 10), y_bounds=(-20, 20), z_bounds=(-10, 10), element="Pb208", priority=1),
        Medium(name="Void", x_bounds=(-20, 20), y_bounds=(-40, 40), z_bounds=(-20, 20), is_void=True, priority=0)
    ]


    # Atomic weight ratio for Pb
    A = 2.5

    #for cadmium
    #A = 2.3

    N = calculate_number_density(11.35, 208) #for Pb as well

    #N = calculate_number_density(8.65, 112) # for CD
    # Simulate particles
    num_particles = 100

    # Initialize Tally object
    tally = Tally()

    # Initialize particle queue
    particle_states = [
        {
            "x": -10.0, "y": 0.0, "z": 0.0,
            "theta": np.random.uniform(0, math.pi),
            "phi": np.random.uniform(0, 2 * math.pi),
            "has_interacted": False,
            "energy": 1e6,  # eV
        }
        for _ in range(num_particles)
    ]

    # Dictionary to store trajectories for all particles
    all_trajectories = {}
    # For tracking if a particle existed within a certain region (x_min, x_max, y_min, y_max, z_min, z_max)
    region_bounds = (14.9, 15.1, -15, 15, -15, 15)

    track = True
    # Prepare arguments for multiprocessing
    args = [(state, reader, mediums, A, N, sampler, region_bounds, track) for state in particle_states]

    # Use multiprocessing to simulate particles
    with Pool() as pool:
        partial_results = pool.map(simulate_single_particle, args)

    # Merge partial results into the main tally
    for idx, result in enumerate(partial_results):
        tally.merge_partial_results(result)

        # Append trajectory to the dictionary
        if result["trajectory"]:
            all_trajectories[idx + 1] = result["trajectory"]

    # Output all trajectories
    with open("all_trajectories.json", "w") as f:
        json.dump(all_trajectories, f, indent=4)
    #print("All particle trajectories:")
    #for particle_id, trajectory in all_trajectories.items():
    #    print(f"Particle {particle_id} trajectory: {trajectory}")

    # Output results
    tally.print_summary(num_particles)

if __name__ == "__main__":
    main()
