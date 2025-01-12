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


def move_to_nearest_boundary(state, mediums, epsilon=1e-6):
    """
    Moves the particle to the nearest boundary and applies a small offset to ensure it's within the next medium.

    Args:
        state: Dictionary representing the particle's state (x, y, z, theta, phi).
        mediums: List of Medium instances.
        epsilon: Small offset to move the particle slightly inside the new medium.

    Returns:
        state: Updated particle state with new (x, y, z) position.
        nearest_medium: The medium the particle is now inside.
    """
    x, y, z = state["x"], state["y"], state["z"]
    u = math.sin(state["theta"]) * math.cos(state["phi"])
    v = math.sin(state["theta"]) * math.sin(state["phi"])
    w = math.cos(state["theta"])

    nearest_distance = float('inf')
    nearest_point = None
    nearest_medium = None

    def compute_intersection(coord, boundary, direction):
        """Helper to compute intersection distance and point."""
        if direction == 0:  # No movement in this direction
            return None
        t = (boundary - coord) / direction  # Distance to boundary
        if t <= 0:  # Only consider forward direction
            return None
        return t

    for medium in mediums:
        x_min, x_max = medium.x_bounds
        y_min, y_max = medium.y_bounds
        z_min, z_max = medium.z_bounds

        # Check intersections with all 6 planes of the medium
        intersections = []

        # X boundaries
        if u != 0:
            t_x_min = compute_intersection(x, x_min, u)
            t_x_max = compute_intersection(x, x_max, u)
            if t_x_min:
                intersections.append((t_x_min, (x_min, y + t_x_min * v, z + t_x_min * w), medium))
            if t_x_max:
                intersections.append((t_x_max, (x_max, y + t_x_max * v, z + t_x_max * w), medium))

        # Y boundaries
        if v != 0:
            t_y_min = compute_intersection(y, y_min, v)
            t_y_max = compute_intersection(y, y_max, v)
            if t_y_min:
                intersections.append((t_y_min, (x + t_y_min * u, y_min, z + t_y_min * w), medium))
            if t_y_max:
                intersections.append((t_y_max, (x + t_y_max * u, y_max, z + t_y_max * w), medium))

        # Z boundaries
        if w != 0:
            t_z_min = compute_intersection(z, z_min, w)
            t_z_max = compute_intersection(z, z_max, w)
            if t_z_min:
                intersections.append((t_z_min, (x + t_z_min * u, y + t_z_min * v, z_min), medium))
            if t_z_max:
                intersections.append((t_z_max, (x + t_z_max * u, y + t_z_max * v, z_max), medium))

        # Filter valid intersections
        for t, point, medium in intersections:
            px, py, pz = point
            # Check if the intersection point is within the medium's bounds
            if medium.x_bounds[0] <= px <= medium.x_bounds[1] and \
               medium.y_bounds[0] <= py <= medium.y_bounds[1] and \
               medium.z_bounds[0] <= pz <= medium.z_bounds[1]:
                # Check if this is the nearest valid intersection
                if t < nearest_distance or (t == nearest_distance and medium.priority > (nearest_medium.priority if nearest_medium else -1)):
                    nearest_distance = t
                    nearest_point = point
                    nearest_medium = medium

    # Move particle to the nearest boundary
    if nearest_point:
        state["x"], state["y"], state["z"] = nearest_point
        # Apply a small offset in the direction of travel
        state["x"] += epsilon * u
        state["y"] += epsilon * v
        state["z"] += epsilon * w

    return state, nearest_medium
