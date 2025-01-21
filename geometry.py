# geometry.py
import numpy as np
from medium import Region


def calculate_direction_cosines(x, y, z, x_prev, y_prev, z_prev):
    delta_x = x - x_prev
    delta_y = y - y_prev
    delta_z = z - z_prev
    delta_s = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
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

def calculate_nearest_boundary(state, regions, u, v, w, epsilon=1e-14):
    """
    Moves the particle to the nearest boundary in a CSG region hierarchy, considering priority.

    Args:
        state: Dictionary representing the particle's state (x, y, z, theta, phi).
        regions: List of regions ordered from most global to most local.
        u, v, w: Direction vector components.
        epsilon: Small number for floating-point comparisons and offset.

    Returns:
        nearest_point: Nearest boundary point (x, y, z).
        nearest_region: Region the particle is now inside.
        nearest_distance: Distance to the nearest boundary.
    """
    x, y, z = state["x"], state["y"], state["z"]
    nearest_distance = float('inf')
    nearest_point = None
    nearest_region = None
    max_priority = -1  # Initialize to the lowest possible priority value

    def solve_boundary_equation(surface, x, y, z, u, v, w):
        """
        Solve for the distance to the surface boundary: f(x + du, y + dv, z + dw) = 0.
        """
        return surface.nearest_surface_method(x, y, z, u, v, w)

    # Iterate through regions from global to local
    for region in regions:
        for surface in region.surfaces:
            # Check if it's a nested region
            if isinstance(surface, Region):
                # Recursively check nested regions
                calculate_nearest_boundary(state, [surface], u, v, w, epsilon)
            else:
                # Compute the distance to the surface
                distance = solve_boundary_equation(surface, x, y, z, u, v, w)
                if distance is not None and distance >= 0:
                    # Handle floating-point equality
                    if abs(distance - nearest_distance) / max(nearest_distance, epsilon) < epsilon:
                        # Prefer surfaces from higher priority regions if coincident
                        continue
                    if distance < nearest_distance:
                        # Compute the intersection point
                        point = (x + distance * u, y + distance * v, z + distance * w)
                        # Check if the point is within the region
                        if region.contains(*point):
                            # If the region has higher priority, select it
                            if region.priority > max_priority:
                                nearest_distance = distance
                                nearest_point = point
                                nearest_region = region
                                max_priority = region.priority

    if nearest_point:
        # Apply a small offset to move across the boundary
        offset_point = (
            nearest_point[0] + epsilon * u,
            nearest_point[1] + epsilon * v,
            nearest_point[2] + epsilon * w
        )
        state.update({"x": offset_point[0], "y": offset_point[1], "z": offset_point[2]})

    return nearest_point, nearest_region, nearest_distance


def calculate_void_si_max(mediums):
    # Find the void medium
    void_medium = next((m for m in mediums if m.is_void), None)
    if void_medium:
        x_range = void_medium.x_bounds[1] - void_medium.x_bounds[0]
        y_range = void_medium.y_bounds[1] - void_medium.y_bounds[0]
        z_range = void_medium.z_bounds[1] - void_medium.z_bounds[0]
        # Diagonal of the bounding box
        si_max = np.sqrt(x_range**2 + y_range**2 + z_range**2)
        return si_max
    return float('inf')  # No void medium, no limit
