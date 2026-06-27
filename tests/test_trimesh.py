import pytest
import numpy as np
import math

# Try importing trimesh; skip tests if not installed
try:
    import trimesh
    from trimesh import transformations
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False

# Import your actual source code
from src.medium import Region, Plane, Sphere, Cylinder, Box
from src.geometry import calculate_nearest_boundary

# ==========================================
# Trimesh Helper Functions (The "Ground Truth")
# ==========================================

def get_trimesh_box_dist(bounds, origin, direction):
    min_bound = np.array(bounds[0])
    max_bound = np.array(bounds[1])
    extents = max_bound - min_bound
    center = (max_bound + min_bound) / 2
    transform = np.eye(4)
    transform[:3, 3] = center

    box = trimesh.primitives.Box(extents=extents, transform=transform)

    # Normalize direction
    direction = np.array(direction, dtype=np.float64)
    norm = np.linalg.norm(direction)
    if norm < 1e-12: return None
    direction = direction / norm

    locations, _, _ = box.ray.intersects_location(
        ray_origins=[origin], ray_directions=[direction]
    )

    if len(locations) > 0:
        return min(np.linalg.norm(loc - origin) for loc in locations)
    return None

def get_trimesh_plane_dist(A, B, C, D, origin, direction):
    direction = np.array(direction, dtype=np.float64)
    norm = np.linalg.norm(direction)
    if norm < 1e-12: return None
    direction = direction / norm

    # Create a line segment long enough to hit the plane
    line_end = origin + direction * 1e6
    endpoints = np.array([origin, line_end])[None].transpose(1, 0, 2)

    plane_normal = np.array([A, B, C], dtype=np.float64)
    plane_normal /= np.linalg.norm(plane_normal)

    # Find a point on the plane
    if A != 0: plane_origin = [-D/A, 0, 0]
    elif B != 0: plane_origin = [0, -D/B, 0]
    elif C != 0: plane_origin = [0, 0, -D/C]
    else: return None

    points, valid = trimesh.intersections.plane_lines(
        plane_origin=plane_origin,
        plane_normal=plane_normal,
        endpoints=endpoints,
        line_segments=True
    )

    if not np.any(valid): return None

    distances = [np.dot(p - origin, direction) for p in points[valid]]
    valid_dists = [d for d in distances if d >= 0]
    return min(valid_dists) if valid_dists else None

def get_trimesh_cylinder_dist(axis, radius, center, origin, direction):
    direction = np.array(direction)
    if np.linalg.norm(direction) == 0: return None
    direction = direction / np.linalg.norm(direction)

    if axis == "z": transform = np.eye(4)
    elif axis == "x": transform = transformations.rotation_matrix(np.pi/2, [0, 1, 0])
    elif axis == "y": transform = transformations.rotation_matrix(np.pi/2, [1, 0, 0])

    transform[:3, 3] = center
    cyl = trimesh.primitives.Cylinder(radius=radius, height=1e6, transform=transform)

    locations, _, _ = cyl.ray.intersects_location([origin], [direction])
    if len(locations) == 0: return None

    t_values = [np.dot(loc - origin, direction) for loc in locations]
    valid_t = [t for t in t_values if t >= 0]
    return min(valid_t) if valid_t else None

def get_trimesh_sphere_dist(center, radius, origin, direction):
    direction = np.array(direction) / np.linalg.norm(direction)
    # High resolution for accuracy
    sphere = trimesh.primitives.Sphere(radius=radius, center=center, subdivisions=4)
    locations, _, _ = sphere.ray.intersects_location([origin], [direction])

    if len(locations) > 0:
        return np.min(np.linalg.norm(locations - origin, axis=1))
    return None

# ==========================================
# Pytest Test Cases
# ==========================================

@pytest.mark.skipif(not TRIMESH_AVAILABLE, reason="Trimesh not installed")
class TestTrimeshComparisons:

    @pytest.mark.parametrize("origin, direction, desc", [
        ([10, 0, 0], [-1, 0, 0], "Starts outside, intersects"),
        ([0, 0, 0], [1, 0, 0], "Starts at center, moves outward"),
        ([0, 0, 4], [0, 0, -1], "Starts on sphere, moves inward"),
        ([0, 0, 5], [1, 0, 0], "Starts on surface, tangent"),
        ([10, 10, 10], [-1, -1, -1], "Starts outside, misses sphere"),
        ([0, 0, 3], [0, 0, 1], "Inside sphere, moves outward"),
        ([5, 0, 0], [0, 1, 0], "Starts on surface, perpendicular"),
        ([-5, 0, 0], [1, 0, 0], "Starts on surface, moves outward"),
        ([10, 10, 10], [1, 1, 1], "Completely misses sphere"),
        ([4.99, 0, 0], [1, 0, 0], "Starts just inside, moves outward"),
    ])
    def test_sphere_vs_trimesh(self, origin, direction, desc):
        center = [0, 0, 0]
        radius = 5

        # 1. Trimesh Calculation
        trimesh_dist = get_trimesh_sphere_dist(center, radius, origin, direction)

        # 2. Your Calculation
        sphere = Sphere(center, radius)
        # Normalize direction for your method
        d = np.array(direction) / np.linalg.norm(direction)
        user_dist = sphere.nearest_surface_method(*origin, *d)

        # 3. Compare
        if trimesh_dist is None or user_dist is None:
            assert trimesh_dist == user_dist
        else:
            # Relaxed tolerance for spheres due to mesh discretization
            assert np.isclose(trimesh_dist, user_dist, rtol=0.01)

    @pytest.mark.parametrize("params, origin, direction, desc", [
        ((0, 0, 1, -5), [0, 0, 0], [0, 0, 1], "Z=5 plane from origin"),
        ((1, 0, 0, -3), [5, 0, 0], [-1, 0, 0], "X=3 plane from right"),
        ((0, 1, 0, 2), [0, -5, 0], [0, 1, 0], "Y=-2 plane from below"),
        ((1, 1, 0, 0), [2, 2, 0], [-1, -1, 0], "Diagonal plane"),
        ((0, 0, 1, 0), [0, 0, 5], [0, 0, -1], "Z=0 plane from above"),
        ((1, 1, 1, 3), [1, 1, 1], [1, 1, 1], "3D plane from inside"),
        ((0, 0, 1, 5), [0, 0, 5], [0, 0, 1], "On-plane movement"),
    ])
    def test_plane_vs_trimesh(self, params, origin, direction, desc):
        A, B, C, D = params

        trimesh_dist = get_trimesh_plane_dist(A, B, C, D, origin, direction)

        # The trimesh helper uses the standard half-space convention
        # Ax+By+Cz+D=0, whereas PyNeut's Plane constructor stores -D
        # (so its inside test is Ax+By+Cz <= D).  Negate D so both describe
        # the same physical plane before comparing.
        plane = Plane(A, B, C, -D)
        d = np.array(direction) / np.linalg.norm(direction)
        user_dist = plane.nearest_surface_method(*origin, *d)

        if trimesh_dist is None or user_dist is None:
            assert trimesh_dist == user_dist
        else:
            assert np.isclose(trimesh_dist, user_dist, rtol=0.001)

    @pytest.mark.parametrize("params, origin, direction, desc", [
        (("z", 3, [0, 0, 0]), [5, 0, 0], [-1, 0, 0], "Z-axis cylinder from right"),
        (("x", 2, [0, 0, 0]), [0, 5, 0], [0, -1, 0], "X-axis cylinder from top"),
        (("y", 4, [0, 0, 0]), [0, 0, 5], [0, 0, -1], "Y-axis cylinder from front"),
        (("z", 2, [0, 0, 0]), [1, 1, 0], [1, 1, 0], "Z-cylinder diagonal"),
        (("x", 3, [0, 0, 0]), [0, 0, 5], [0, 0, -1], "X-cylinder along axis"),
    ])
    def test_cylinder_vs_trimesh(self, params, origin, direction, desc):
        axis, radius, center = params

        trimesh_dist = get_trimesh_cylinder_dist(axis, radius, center, origin, direction)

        cylinder = Cylinder(axis, radius, center)
        d = np.array(direction) / np.linalg.norm(direction)
        user_dist = cylinder.nearest_surface_method(*origin, *d)

        if trimesh_dist is None or user_dist is None:
            assert trimesh_dist == user_dist
        else:
            assert np.isclose(trimesh_dist, user_dist, rtol=0.001)

    @pytest.mark.parametrize("origin, direction, desc", [
        ([5, 0, 0], [-1, 0, 0], "X-axis entry from right"),
        ([-5, 0, 0], [1, 0, 0], "X-axis entry from left"),
        ([0, 5, 0], [0, -1, 0], "Y-axis entry from top"),
        ([0, 0, 6], [0, 0, -1], "Z-axis entry from front"),
        ([3, 4, 5], [-1, -1, -1], "Diagonal entry"),
        ([0, 0, 0], [1, 0, 0], "Starts inside, exits right"),
        ([2.1, 0, 0], [-1, 0, 0], "Starts just outside right face"),
        ([1, 2, 3], [0.5, 0.5, 0.5], "Angled exit from inside"),
        ([5, 5, 5], [-1, -1, -1], "Diagonal miss"),
    ])
    def test_box_vs_trimesh(self, origin, direction, desc):
        # Box bounds: [-2, -3, -4] to [2, 3, 4]
        x_min, x_max = -2, 2
        y_min, y_max = -3, 3
        z_min, z_max = -4, 4
        bounds = [[x_min, y_min, z_min], [x_max, y_max, z_max]]

        # 1. Trimesh
        trimesh_dist = get_trimesh_box_dist(bounds, origin, direction)

        # 2. User
        custom_box = Box(x_min, x_max, y_min, y_max, z_min, z_max)
        state = {"x": origin[0], "y": origin[1], "z": origin[2]}
        d = np.array(direction) / np.linalg.norm(direction)

        # Using calculate_nearest_boundary as per your original script
        _, _, user_dist = calculate_nearest_boundary(
            state, [custom_box], d[0], d[1], d[2]
        )

        # Handle inf vs None
        if user_dist == float('inf'):
            user_dist = None

        if trimesh_dist is None or user_dist is None:
            assert trimesh_dist == user_dist
        else:
            assert np.isclose(trimesh_dist, user_dist, rtol=0.001)
