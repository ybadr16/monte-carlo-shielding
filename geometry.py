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

def calculate_nearest_boundary(state, regions, u, v, w, epsilon=1e-10):
    """
    Calculate the nearest boundary and determine the new region/medium.

    Args:
        state: Dictionary representing the particle's state (x, y, z, theta, phi).
        regions: List of regions ordered from global to local.
        u, v, w: Direction vector components (particle's velocity).
        epsilon: Small number for floating-point comparisons.

    Returns:
        nearest_point: Coordinates of the nearest boundary point (x, y, z).
        nearest_region: Region the particle will move into.
        nearest_distance: Distance to the nearest boundary.
    """
    x, y, z = state["x"], state["y"], state["z"]
    nearest_distance = float("inf")
    nearest_point = None
    nearest_region = None
    nearest_surface = None  # Track the nearest surface
    max_priority = -1

    def solve_boundary_equation(surface, x, y, z, u, v, w):
        """
        Solve for the distance to the surface boundary f(x + du, y + dv, z + dw) = 0.
        Returns the distance to the surface from the current position.
        """
        return surface.nearest_surface_method(x, y, z, u, v, w)

    for region in regions:
        for surface in region.surfaces:
            # Compute distance to the surface
            distance = solve_boundary_equation(surface, x, y, z, u, v, w)

            # Validate the distance
            if distance is not None and distance >= 0:  # Forward intersections only
                # Compute the intersection point
                point = (x + distance * u, y + distance * v, z + distance * w)

                # Use the region's `contains` method to validate the point
                if region.contains(*point):
                    # Handle priority for coincident boundaries
                    if abs(distance - nearest_distance) < epsilon:
                        if region.priority > max_priority:
                            nearest_distance = distance
                            nearest_point = point
                            nearest_region = region
                            max_priority = region.priority
                    elif distance < nearest_distance:
                        nearest_distance = distance
                        nearest_point = point
                        nearest_region = region
                        max_priority = region.priority

    if nearest_point is not None and nearest_region is not None and nearest_surface is not None:
        new_x = nearest_point[0]
        new_y = nearest_point[1]
        new_z = nearest_point[2]

        # Find the new region
        for region in regions:
            if region.contains(new_x, new_y, new_z):
                nearest_region = region
                break
        nearest_point = (new_x, new_y, new_z)

    # Handle case where no boundary is found
    if nearest_point is None:
        # Particle escapes geometry; handle accordingly
        return None, None, float('inf')

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
