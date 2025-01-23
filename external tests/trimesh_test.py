import trimesh
import numpy as np
import math
from trimesh import transformations

def benchmark_box(bounds, origin, direction):
    """
    Compare intersections between Trimesh's Box and our custom Box region.
    :param bounds: List [[x_min, y_min, z_min], [x_max, y_max, z_max]]
    :param origin: Ray origin (x, y, z)
    :param direction: Ray direction (u, v, w)
    :return: Trimesh distance, User distance
    """
    # Create Trimesh Box from bounds using extents and transform
    min_bound = np.array(bounds[0])
    max_bound = np.array(bounds[1])
    extents = max_bound - min_bound
    center = (max_bound + min_bound) / 2
    transform = np.eye(4)
    transform[:3, 3] = center  # Translate box to the correct center
    trimesh_box = trimesh.primitives.Box(extents=extents, transform=transform)

    # Create custom Box region
    x_min, y_min, z_min = bounds[0]
    x_max, y_max, z_max = bounds[1]
    custom_box = Box(x_min, x_max, y_min, y_max, z_min, z_max)

    # Normalize direction for both methods
    direction = np.array(direction, dtype=np.float64)
    dir_norm = np.linalg.norm(direction)
    if dir_norm < 1e-12:
        return None, None
    direction = direction / dir_norm

    # Compute Trimesh intersection
    locations, _, _ = trimesh_box.ray.intersects_location(
        ray_origins=[origin],
        ray_directions=[direction]
    )
    trimesh_dist = min(np.linalg.norm(loc - origin) for loc in locations) if len(locations) > 0 else None

    # Compute custom Box intersection
    state = {"x": origin[0], "y": origin[1], "z": origin[2]}
    _, _, user_dist = calculate_nearest_boundary(state, [custom_box], *direction)

    return trimesh_dist, user_dist

def ray_plane_intersection(A, B, C, D, origin, direction):
    # Normalize direction vector
    direction = np.array(direction, dtype=np.float64)
    dir_norm = np.linalg.norm(direction)
    if dir_norm < 1e-12:
        return None
    direction = direction / dir_norm

    # Create line segment with correct shape (2, 1, 3)
    line_end = origin + direction * 1e6
    endpoints = np.array([origin, line_end])[None].transpose(1, 0, 2)  # Shape (2, 1, 3)

    # Create plane normal and origin
    plane_normal = np.array([A, B, C], dtype=np.float64)
    plane_normal /= np.linalg.norm(plane_normal)

    if A != 0:
        plane_origin = np.array([-D/A, 0, 0], dtype=np.float64)
    elif B != 0:
        plane_origin = np.array([0, -D/B, 0], dtype=np.float64)
    elif C != 0:
        plane_origin = np.array([0, 0, -D/C], dtype=np.float64)
    else:
        raise ValueError("Invalid plane equation")

    # Get intersections
    points, valid = trimesh.intersections.plane_lines(
        plane_origin=plane_origin,
        plane_normal=plane_normal,
        endpoints=endpoints,
        line_segments=True
    )

    if not np.any(valid):
        return None

    # Calculate distance for valid intersections
    distances = [np.dot(p - origin, direction) for p in points[valid]]
    return min([d for d in distances if d >= 0], default=None)


def ray_cylinder_intersection(axis, radius, center, origin, direction, height=1e6):
    direction = np.array(direction)
    dir_norm = np.linalg.norm(direction)
    if dir_norm == 0:
        return None
    direction = direction / dir_norm

    if axis == "z":
        transform = np.eye(4)
    elif axis == "x":
        transform = transformations.rotation_matrix(np.pi/2, [0, 1, 0])
    elif axis == "y":
        transform = transformations.rotation_matrix(np.pi/2, [1, 0, 0])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")

    transform[:3, 3] = center
    cylinder = trimesh.primitives.Cylinder(
        radius=radius,
        height=height,
        transform=transform
    )

    locations, _, _ = cylinder.ray.intersects_location([origin], [direction])
    if len(locations) == 0:
        return None

    t_values = [np.dot(loc - origin, direction) for loc in locations]
    valid_t = [t for t in t_values if t >= 0]
    return min(valid_t) if valid_t else None

class Region:
    def __init__(self, surfaces=None, operation="intersection", name=None, priority=0, is_void=False, element=None):
        """
        Create a region combining surfaces using boolean operations.
        :param surfaces: List of surfaces defining the region.
        :param operation: Boolean operation ('intersection', 'union', 'complement').
        :param name: Name of the region.
        :param priority: Priority of the region (higher value means higher priority).
        :param is_void: Whether the region is a void (no material interaction).
        :param element: The material element associated with the region.
        """
        self.surfaces = surfaces if surfaces else []
        self.operation = operation
        self.name = name
        self.priority = priority
        self.is_void = is_void
        self.element = element

    def contains(self, x, y, z):
        evaluations = []
        for surface in self.surfaces:
            if isinstance(surface, Region):
                evaluations.append(surface.contains(x, y, z))
            else:
                eval_result = surface.evaluate(x, y, z) <= 0
                evaluations.append(eval_result)

        if self.operation == "intersection":
            return all(evaluations)
        elif self.operation == "union":
            return any(evaluations)
        elif self.operation == "complement":
            return not evaluations[0]
        elif self.operation == "difference":
            # A - B = A ∩ ¬B
            if len(self.surfaces) != 2:
                raise ValueError("Difference operation requires exactly two surfaces")
            a = evaluations[0]
            b = evaluations[1]
            return a and not b
        else:
            raise ValueError(f"Unknown operation: {self.operation}")

    def add_surface(self, surface):
        """Add a surface to the region."""
        self.surfaces.append(surface)



class Plane:
    def __init__(self, A, B, C, D):
        # Normalize the normal vector (A, B, C)
        norm = np.sqrt(A**2 + B**2 + C**2)
        if np.isclose(norm, 0):
            raise ValueError("Plane normal cannot be a zero vector")

        self.A = A / norm
        self.B = B / norm
        self.C = C / norm
        self.D = -D / norm  # Adjusted for normalized normal

    def evaluate(self, x, y, z):
        #print(self.A * x + self.B * y + self.C * z - self.D)
        return self.A * x + self.B * y + self.C * z + self.D

    def nearest_surface_method(self, x, y, z, u, v, w):
        numerator = self.D - self.A*x - self.B*y - self.C*z
        denominator = self.A*u + self.B*v + self.C*w

        if np.isclose(numerator, 0, atol=1e-8):
            return 0.0 if not np.isclose(denominator, 0) else None

        if denominator == 0:
            return None

        t = numerator / denominator
        return t if t >= 0 else None

    def normal(self, x, y, z):
        mag = (self.A**2 + self.B**2 + self.C**2) ** 0.5
        if mag == 0:
            return (0.0, 0.0, 0.0)  # Handle invalid plane
        return (self.A / mag, self.B / mag, self.C / mag)

class Cylinder:
    def __init__(self, axis, radius, center):
        self.axis = axis
        self.radius = radius
        self.x0, self.y0, self.z0 = center

    def evaluate(self, x, y, z):
        if self.axis == "z":
            return (x - self.x0) ** 2 + (y - self.y0) ** 2 - self.radius ** 2
        elif self.axis == "x":
            return (y - self.y0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2
        elif self.axis == "y":
            return (x - self.x0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2

    def nearest_surface_method(self, x, y, z, u, v, w):
        if self.axis == "z":
            x_bar = x - self.x0
            y_bar = y - self.y0
            a = u**2 + v**2
            k = x_bar * u + y_bar * v
            c = x_bar**2 + y_bar**2 - self.radius**2
            radial_u, radial_v = x_bar, y_bar  # Radial vector components
            dir_u, dir_v = u, v                # Direction components
        elif self.axis == "x":
            y_bar = y - self.y0
            z_bar = z - self.z0
            a = v**2 + w**2
            k = y_bar * v + z_bar * w
            c = y_bar**2 + z_bar**2 - self.radius**2
            radial_u, radial_v = y_bar, z_bar
            dir_u, dir_v = v, w
        elif self.axis == "y":
            x_bar = x - self.x0
            z_bar = z - self.z0
            a = u**2 + w**2
            k = x_bar * u + z_bar * w
            c = x_bar**2 + z_bar**2 - self.radius**2
            radial_u, radial_v = x_bar, z_bar
            dir_u, dir_v = u, w

        if a == 0:
            return None  # Direction parallel to cylinder axis

        discriminant = k**2 - a * c
        if discriminant < 0:
            return None  # No real intersection

        sqrt_d = math.sqrt(discriminant)
        d1 = (-k - sqrt_d) / a
        d2 = (-k + sqrt_d) / a

        if c < 0:
            # Inside cylinder: return exit point (max root)
            return max(d1, d2)
        elif c == 0:
            # On the surface: check direction relative to radial vector
            dot = radial_u * dir_u + radial_v * dir_v
            if dot < 0:
                # Moving inward: return exit point (d2 if valid)
                return d2 if d2 >= 0 else None
            else:
                # Moving outward/tangentially: no future intersection
                return None
        else:
            # Outside: return first valid intersection (min root)
            valid = [d for d in (d1, d2) if d >= 0]
            return min(valid) if valid else None

        def normal(self, x, y, z):
            if self.axis == 'z':
                dx = x - self.x0
                dy = y - self.y0
                mag = (dx**2 + dy**2)**0.5
                if mag == 0:
                    return (0.0, 0.0, 0.0)  # Handle edge case
                return (dx/mag, dy/mag, 0.0)
            elif self.axis == 'x':
                dy = y - self.y0
                dz = z - self.z0
                mag = (dy**2 + dz**2)**0.5
                if mag == 0:
                    return (0.0, 0.0, 0.0)
                return (0.0, dy/mag, dz/mag)
            elif self.axis == 'y':
                dx = x - self.x0
                dz = z - self.z0
                mag = (dx**2 + dz**2)**0.5
                if mag == 0:
                    return (0.0, 0.0, 0.0)
                return (dx/mag, 0.0, dz/mag)

class Sphere:
    def __init__(self, center, radius):
        self.x0, self.y0, self.z0 = center
        self.radius = radius

    def evaluate(self, x, y, z):
        return (x - self.x0) ** 2 + (y - self.y0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2

    def nearest_surface_method(self, x, y, z, u, v, w):
        dir_norm = (u**2 + v**2 + w**2)**0.5
        if dir_norm == 0:
            return None

        u_norm = u / dir_norm
        v_norm = v / dir_norm
        w_norm = w / dir_norm

        x_bar = x - self.x0
        y_bar = y - self.y0
        z_bar = z - self.z0

        k = x_bar * u_norm + y_bar * v_norm + z_bar * w_norm
        c = x_bar**2 + y_bar**2 + z_bar**2 - self.radius**2

        discriminant = k**2 - c
        if discriminant < 0:
            return None

        sqrt_d = discriminant**0.5
        d1 = -k - sqrt_d
        d2 = -k + sqrt_d

        if c < 0:
            return max(d1, d2)
        else:
            return d1 if d1 >= 0 else (d2 if d2 >= 0 else None)

    def normal(self, x, y, z):
        dx = x - self.x0
        dy = y - self.y0
        dz = z - self.z0
        mag = (dx**2 + dy**2 + dz**2) ** 0.5
        if mag == 0:
            return (0.0, 0.0, 0.0)  # Handle center point
        return (dx / mag, dy / mag, dz / mag)

def benchmark_plane(plane_params, origin, direction):
    A, B, C, D = plane_params
    direction = np.array(direction)
    dir_norm = np.linalg.norm(direction)
    if dir_norm == 0:
        return None, None
    direction = direction / dir_norm

    plane = Plane(A, B, C, D)
    trimesh_t = ray_plane_intersection(A, B, C, D, origin, direction)
    user_t = plane.nearest_surface_method(*origin, *direction)
    return trimesh_t, user_t

def benchmark_cylinder(cylinder_params, origin, direction):
    axis, radius, center = cylinder_params
    direction = np.array(direction)
    dir_norm = np.linalg.norm(direction)
    if dir_norm == 0:
        return None, None
    direction = direction / dir_norm

    cylinder = Cylinder(axis, radius, center)
    trimesh_t = ray_cylinder_intersection(axis, radius, center, origin, direction)
    user_t = cylinder.nearest_surface_method(*origin, *direction)
    return trimesh_t, user_t

def ray_sphere_intersection(center, radius, origin, direction):
    direction = np.array(direction) / np.linalg.norm(direction)
    # Create a high-resolution sphere with subdivisions=5
    sphere = trimesh.primitives.Sphere(radius=radius, center=center, subdivisions=5)
    locations, index_rays, _ = sphere.ray.intersects_location([origin], [direction])
    return np.min(np.linalg.norm(locations - origin, axis=1)) if len(locations) > 0 else None


def benchmark_sphere_methods(center, radius, origin, direction):
    sphere = Sphere(center, radius)
    trimesh_distance = ray_sphere_intersection(center, radius, origin, direction)
    user_distance = sphere.nearest_surface_method(*origin, *direction)

    # Handle None comparisons
    if trimesh_distance is None or user_distance is None:
        assert trimesh_distance == user_distance, f"Trimesh: {trimesh_distance}, User: {user_distance}"
    else:
        # Increase tolerance to 0.1% to account for mesh approximation
        assert np.isclose(trimesh_distance, user_distance, rtol=0.001), f"Trimesh: {trimesh_distance}, User: {user_distance}"

    # Safe formatting for None values
    trimesh_str = f"{trimesh_distance:.2f}" if trimesh_distance is not None else "None"
    user_str = f"{user_distance:.2f}" if user_distance is not None else "None"
    print(f"Trimesh: {trimesh_str}, User: {user_str} ✓")

def benchmark_region(region, origin, direction):
    # Convert direction to components
    u, v, w = direction
    dir_norm = np.linalg.norm(direction)
    if dir_norm == 0:
        return None, None
    u, v, w = direction / dir_norm

    # Calculate using our method
    state = {"x": origin[0], "y": origin[1], "z": origin[2]}
    _, _, user_dist = calculate_nearest_boundary(state, [region], u, v, w)

    # Calculate using Trimesh (needs to handle multiple surfaces)
    trimesh_dists = []
    for surface in region.surfaces:
        if isinstance(surface, Cylinder):
            dist = ray_cylinder_intersection(
                surface.axis, surface.radius,
                [surface.x0, surface.y0, surface.z0],
                origin, direction
            )
        elif isinstance(surface, Plane):
            dist = ray_plane_intersection(
                surface.A, surface.B, surface.C, surface.D,
                origin, direction
            )
        if dist is not None:
            trimesh_dists.append(dist)

    trimesh_dist = min(trimesh_dists) if trimesh_dists else None

    return trimesh_dist, user_dist


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



class Box(Region):
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        """
        Define a box region bounded by planes at x_min, x_max, y_min, y_max, z_min, and z_max.
        """
        planes = [
            Plane(-1, 0, 0, -x_min),  # x >= x_min: -x + x_min ≤ 0 → x ≥ x_min
            Plane(1, 0, 0, x_max),    # x <= x_max: x - x_max ≤ 0 → x ≤ x_max
            Plane(0, -1, 0, -y_min),  # y >= y_min: -y + y_min ≤ 0 → y ≥ y_min
            Plane(0, 1, 0, y_max),    # y <= y_max: y - y_max ≤ 0 → y ≤ y_max
            Plane(0, 0, -1, -z_min),  # z >= z_min: -z + z_min ≤ 0 → z ≥ z_min
            Plane(0, 0, 1, z_max),    # z <= z_max: z - z_max ≤ 0 → z ≤ z_max
        ]
        super().__init__(surfaces=planes, operation="intersection", name="Box")


Bounded_cylinder = Region(
    surfaces=[
        Cylinder("z", 10, (0, 0, 0)),
        Plane(0, 0, 1, 500000),  # z = -500000
        Plane(0, 0, -1, 500000)  # z = 500000
    ],
    operation="intersection",
    name="Finite Cylinder"
)


def run_tests():
    # Sphere Test Cases
    print("\n=== Testing Sphere Intersections ===")
    sphere_center = [0, 0, 0]
    sphere_radius = 5
    sphere_tests = [
        ([10, 0, 0], [-1, 0, 0], "Starts outside, intersects"),
        ([0, 0, 0], [1, 0, 0], "Starts at the center, moves outward"),
        ([0, 0, 4], [0, 0, -1], "Starts on the sphere, moves inward"),
        ([0, 0, 5], [1, 0, 0], "Starts on the surface, tangent"),
        ([10, 10, 10], [-1, -1, -1], "Starts outside, misses sphere"),
        ([0, 0, 3], [0, 0, 1], "Inside sphere, moves outward"),
        ([5, 0, 0], [0, 1, 0], "Starts on surface, perpendicular direction"),
        ([-5, 0, 0], [1, 0, 0], "Starts on surface, moves outward"),
        ([10, 10, 10], [1, 1, 1], "Completely misses the sphere"),
        ([5, 5, 5], [-1, -1, -1], "Starts outside, intersects at an angle"),
        ([0, 0, 0], [-1, 1, 0], "Starts at the center, non-axial direction"),
        ([6, 6, 0], [-1, -1, 0], "Starts outside, grazes tangent"),
        ([4.99, 0, 0], [1, 0, 0], "Starts just inside the surface, moves outward"),
        ([5, 0, 0], [-1, 0, 0], "Starts exactly on the surface, moves inward"),
        ([0, 0, -6], [0, 0, 1], "Starts just outside along z-axis"),
    ]

    for i, (origin, direction, desc) in enumerate(sphere_tests):
        print(f"Test {i+1}: {desc}")
        try:
            benchmark_sphere_methods(sphere_center, sphere_radius, origin, direction)
        except AssertionError as e:
            print(f"  → Failed: {e}")

    # Plane Test Cases
    print("\n=== Testing Plane Intersections ===")
    plane_tests = [
        ((0, 0, 1, -5), [0, 0, 0], [0, 0, 1], "Z=5 plane from origin"),
        ((1, 0, 0, -3), [5, 0, 0], [-1, 0, 0], "X=3 plane from right"),
        ((0, 1, 0, 2), [0, -5, 0], [0, 1, 0], "Y=-2 plane from below"),
        ((1, 1, 0, 0), [2, 2, 0], [-1, -1, 0], "Diagonal plane from right-top"),
        ((0, 0, 1, 0), [0, 0, 5], [0, 0, -1], "Z=0 plane from above"),
        ((1, 1, 1, 3), [1, 1, 1], [1, 1, 1], "3D plane from inside"),
        ((0, 0, 1, 5), [0, 0, 5], [0, 0, 1], "On-plane movement"),
        ((1, 0, 0, -5), [10, 0, 0], [-1, 0, 0], "X=5 plane from far right"),
        ((0, 1, 0, 0), [0, -5, 0], [0, 1, 1], "Y=0 plane with angled approach"),
        ((2, -1, 3, 4), [1, 2, 3], [2, -1, 3], "Arbitrary plane from nearby")
    ]

    # In the plane test loop:
    for i, (params, origin, direction, desc) in enumerate(plane_tests):
        print(f"\nPlane Test {i+1}: {desc}")
        try:
            trimesh_dist, user_dist = benchmark_plane(params, origin, direction)

            # Safe formatting for None values
            trimesh_str = f"{trimesh_dist:.2f}" if trimesh_dist is not None else "None"
            user_str = f"{user_dist:.2f}" if user_dist is not None else "None"
            print(f"  ✓ Trimesh: {trimesh_str}, User: {user_str}")

            # Numerical comparison only if both values exist
            if trimesh_dist is not None and user_dist is not None:
                assert np.isclose(trimesh_dist, user_dist, rtol=0.001), \
                       f"Discrepancy: {trimesh_dist:.2f} vs {user_dist:.2f}"
            else:
                assert trimesh_dist == user_dist, "Mismatch in None results"

        except AssertionError as e:
            print(f"  → Failed: {e}")

    # Cylinder Test Cases
    print("\n=== Testing Cylinder Intersections ===")
    cylinder_tests = [
        (("z", 3, [0, 0, 0]), [5, 0, 0], [-1, 0, 0], "Z-axis cylinder from right"),
        (("x", 2, [0, 0, 0]), [0, 5, 0], [0, -1, 0], "X-axis cylinder from top"),
        (("y", 4, [0, 0, 0]), [0, 0, 5], [0, 0, -1], "Y-axis cylinder from front"),
        (("z", 2, [0, 0, 0]), [1, 1, 0], [1, 1, 0], "Z-cylinder diagonal approach"),
        (("x", 3, [0, 0, 0]), [0, 0, 5], [0, 0, -1], "X-cylinder along axis"),
        (("y", 2.5, [0, 0, 0]), [3, 0, 0], [-1, 0, 0], "Y-cylinder radial approach"),
        ((Bounded_cylinder), [0, 0, 0], [0, 0, 1], "Bounded cylinder from center"),
        (("x", 4, [0, 0, 0]), [5, 0, 0], [-1, 1, 0], "X-cylinder angled approach"),
        (("y", 3, [0, 0, 0]), [2, 0, 2], [0, 0, -1], "Y-cylinder edge case"),
        (("z", 2, [0, 0, 0]), [3, 0, 0], [-1, 0, 0], "Z-cylinder tangent test")
    ]

    for i, (region, origin, direction, desc) in enumerate(cylinder_tests):
        print(f"\nCylinder Test {i+1}: {desc}")
        try:
            # For bounded cylinder test
            if isinstance(region, Region):
                trimesh_dist, user_dist = benchmark_region(region, origin, direction)
            else:
                 # Regular cylinder test
                trimesh_dist, user_dist = benchmark_cylinder(region, origin, direction)

            # Safe formatting for None values
            trimesh_str = f"{trimesh_dist:.2f}" if trimesh_dist is not None else "None"
            user_str = f"{user_dist:.2f}" if user_dist is not None else "None"
            print(f"  ✓ Trimesh: {trimesh_str}, User: {user_str}")

            # Numerical comparison only if both values exist
            if trimesh_dist is not None and user_dist is not None:
                assert np.isclose(trimesh_dist, user_dist, rtol=0.01), \
                       f"Discrepancy: {trimesh_dist:.2f} vs {user_dist:.2f}"
            else:
                assert trimesh_dist == user_dist, "Mismatch in None results"

        except AssertionError as e:
            print(f"  → Failed: {e}")


    print("\n=== Testing Box Regions ===")
    box_bounds = [[-2, -3, -4], [2, 3, 4]]  # Bounds for a 4x6x8 box centered at (0,0,0)
    box_tests = [
        ([5, 0, 0], [-1, 0, 0], "X-axis entry from right"),
        ([-5, 0, 0], [1, 0, 0], "X-axis entry from left"),
        ([0, 5, 0], [0, -1, 0], "Y-axis entry from top"),
        ([0, -5, 0], [0, 1, 0], "Y-axis entry from bottom"),
        ([0, 0, 6], [0, 0, -1], "Z-axis entry from front"),
        ([0, 0, -6], [0, 0, 1], "Z-axis entry from back"),
        ([3, 4, 5], [-1, -1, -1], "Diagonal entry from outside"),
        ([0, 0, 0], [1, 0, 0], "Starts inside, exits right"),
        ([0, 0, 0], [0, 1, 0], "Starts inside, exits top"),
        ([2.1, 0, 0], [-1, 0, 0], "Starts just outside right face"),
        ([-2.1, 0, 0], [1, 0, 0], "Starts just outside left face"),
        ([1, 2, 3], [0.5, 0.5, 0.5], "Angled exit from inside"),
        ([5, 5, 5], [-1, -1, -1], "Diagonal miss (no intersection)")
    ]

    for i, (origin, direction, desc) in enumerate(box_tests):
        print(f"\nBox Test {i+1}: {desc}")
        try:
            trimesh_dist, user_dist = benchmark_box(box_bounds, origin, direction)

            trimesh_str = f"{trimesh_dist:.2f}" if trimesh_dist is not None else "None"
            user_str = f"{user_dist:.2f}" if user_dist is not None else "None"
            print(f"  ✓ Trimesh: {trimesh_str}, User: {user_str}")

            if trimesh_dist is not None and user_dist is not None:
                assert np.isclose(trimesh_dist, user_dist, rtol=0.001), \
                    f"Discrepancy: {trimesh_dist:.2f} vs {user_dist:.2f}"
            else:
                assert trimesh_dist == user_dist, "Mismatch in None results"
        except AssertionError as e:
            print(f"  → Failed: {e}")
# Run all tests
run_tests()
