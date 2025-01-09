import json
import random
import matplotlib.pyplot as plt

# Load JSON data
with open("all_trajectories.json", "r") as f:
    data = json.load(f)

# Define shield bounds
shield_x_bounds = (-10, 10)
shield_y_bounds = (-20, 20)
shield_z_bounds = (-10, 10)

# Function to check if a position is within the shield bounds
def within_bounds(position):
    x, y, z = position
    return (
        shield_x_bounds[0] < x <= shield_x_bounds[1]
        and shield_y_bounds[0] < y <= shield_y_bounds[1]
        and shield_z_bounds[0] < z <= shield_z_bounds[1]
    )

# Function to classify a particle's trajectory
def is_absorbed(history):
    # Absorbed if the particle stays within bounds and ends inside
    return all(within_bounds(pos) for pos in history) and within_bounds(history[-1])

def is_scattered(history):
    # Scattered if the particle leaves the bounds at any point
    return any(not within_bounds(pos) for pos in history)

# Filter the particles
histories = list(data.values())
absorbed_particles = []
scattered_particles = []

for history in histories:
    if is_absorbed(history):
        absorbed_particles.append(history)
    elif is_scattered(history):
        scattered_particles.append(history)

# Randomly select a few particles for plotting
random_absorbed = random.sample(absorbed_particles, min(5, len(absorbed_particles)))
random_scattered = random.sample(scattered_particles, min(20, len(scattered_particles)))

# Plot the selected trajectories
plt.figure(figsize=(10, 8))

# Plot absorbed particles
for i, history in enumerate(random_absorbed):
    x_coords = [pos[0] for pos in history]
    y_coords = [pos[1] for pos in history]
    plt.plot(x_coords, y_coords, label=f"Absorbed {i+1}", linestyle="dotted", color="blue")

# Plot scattered particles
for i, history in enumerate(random_scattered):
    x_coords = [pos[0] for pos in history]
    y_coords = [pos[1] for pos in history]
    plt.plot(x_coords, y_coords, label=f"Scattered {i+1}", linestyle="solid", color="red")

# Add shield bounds for reference
plt.axvline(x=shield_x_bounds[0], color="gray", linestyle="--", label="Shield X Bound")
plt.axvline(x=shield_x_bounds[1], color="gray", linestyle="--")
plt.axhline(y=shield_y_bounds[0], color="gray", linestyle="--", label="Shield Y Bound")
plt.axhline(y=shield_y_bounds[1], color="gray", linestyle="--")

# Formatting the plot
plt.title("2D Trajectories of Particles (Scattered and Absorbed)")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.legend()
plt.grid()
plt.show()
