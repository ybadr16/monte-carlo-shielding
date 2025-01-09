# geometry.py
import math

def calculate_direction_cosines(x, y, z, x_prev, y_prev, z_prev):
    delta_x = x - x_prev
    delta_y = y - y_prev
    delta_z = z - z_prev
    delta_s = math.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
    return delta_x / delta_s, delta_y / delta_s, delta_z / delta_s

def count_coordinates_in_boundary(coordinates, x_bounds, y_bounds, z_bounds):
    """
    Counts the number of coordinates within the specified boundaries.

    Args:
        coordinates (list of tuples): List of (x, y, z) coordinate tuples.
        x_bounds (tuple): Min and max values for the x-coordinate.
        y_bounds (tuple): Min and max values for the y-coordinate.
        z_bounds (tuple): Min and max values for the z-coordinate.

    Returns:
        int: Number of coordinates within the specified boundaries.
    """
    return sum(
        x_bounds[0] <= x <= x_bounds[1] and
        y_bounds[0] <= y <= y_bounds[1] and
        z_bounds[0] <= z <= z_bounds[1]
        for x, y, z in coordinates
    )


def calculate_nearest_boundary_distance(state, medium):
    """
    Calculates the nearest distance from the current particle position to the boundary of a medium.

    Args:
        state: Dictionary containing the particle's current state (keys: "x", "y", "z").
        medium: Medium instance representing the medium of interest.

    Returns:
        The nearest distance to the boundary of the medium.
    """
    # Extract the current position of the particle from the state
    x, y, z = state["x"], state["y"], state["z"]

    # Extract the bounds of the medium
    x_min, x_max = medium.x_bounds
    y_min, y_max = medium.y_bounds
    z_min, z_max = medium.z_bounds

    # Calculate the distance to each boundary in all three dimensions
    distances = [
        max(0, x - x_min),  # Distance to the x_min boundary
        max(0, x_max - x),  # Distance to the x_max boundary
        max(0, y - y_min),  # Distance to the y_min boundary
        max(0, y_max - y),  # Distance to the y_max boundary
        max(0, z - z_min),  # Distance to the z_min boundary
        max(0, z_max - z)   # Distance to the z_max boundary
    ]

    # Return the smallest positive distance to any boundary
    return min(distances)
