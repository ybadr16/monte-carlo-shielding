import pytest
from src.medium import Cylinder
import math

def direction_cosines(x, y, z, x_o, y_o, z_o):
    # Calculate the displacement vector components
    dx = x - x_o
    dy = y - y_o
    dz = z - z_o

    # Compute the magnitude of the displacement vector
    magnitude = math.sqrt(dx**2 + dy**2 + dz**2)

    # Calculate the direction cosines (normalized components of the displacement vector)
    u = dx / magnitude
    v = dy / magnitude
    w = dz / magnitude

    return u, v, w


def test_evaluate():
    # Test evaluate method for a z-axis cylinder
    cylinder = Cylinder("z", 1, (0,0,0))
    assert cylinder.evaluate(0, 0, 0) == -1  # Inside cylinder
    assert cylinder.evaluate(1, 0, 0) == 0   # On surface
    assert cylinder.evaluate(2, 0, 0) == 3   # Outside cylinder

    # Test evaluate method for an x-axis cylinder
    cylinder = Cylinder("x", 2, (0,1,0))
    assert cylinder.evaluate(0, 3, 0) == 0   # On surface
    assert cylinder.evaluate(0, 4, 0) == 5   # Outside cylinder
    assert cylinder.evaluate(0, 1, 0) == -4  # Inside cylinder

def test_nearest_surface():
    # Case 1: Particle outside, moving toward z-axis cylinder
    cylinder = Cylinder("z", 1, (0,0,0))
    assert abs(cylinder.nearest_surface_method(2, 0, 0, -1, 0, 0) - 1) < 1e-6  # Should hit at x=1

    # Case 2: Particle inside, moving outward
    assert abs(cylinder.nearest_surface_method(0.5, 0, 0, 1, 0, 0) - 0.5) < 1e-6  # Hits at x=1

    # Case 3: Particle parallel to the cylinder, no intersection
    assert cylinder.nearest_surface_method(2, 0, 0, 0, 0, 1) is None

    # Case 4: Particle starting on the surface, moving inward
    assert abs(cylinder.nearest_surface_method(1, 0, 0, -1, 0, 0) - 2) < 1e-6  # Crosses center and hits other side

    # Case 5: Particle on the axis of the cylinder
    assert cylinder.nearest_surface_method(0, 0, 0, 0, 0, -1) is None  # Parallel, no intersection


    u, v, w = direction_cosines(0, 0, 0, 0.5, 0.5, 0)

    A = u**2 + v**2
    B = 2 * (0.5 * u + 0.5 * v)
    C = 0.5**2 + 0.5**2 - 1

    # Solve quadratic equation A*t^2 + B*t + C = 0
    discriminant = B**2 - 4 * A * C
    assert discriminant >= 0, "No real intersection!"

    t1 = (-B + math.sqrt(discriminant)) / (2 * A)
    t2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # Use the positive root (moving outward)
    expected_t = max(t1, t2)
    # Case 6: Particle starting inside the cylinder
    assert abs(cylinder.nearest_surface_method(0.5, 0.5, 0, u, v, 0) - expected_t) < 1e-4  # Diagonal movement out

    # Case 7: x-axis cylinder
    cylinder = Cylinder("x", 1, (0, 0, 0))
    assert abs(cylinder.nearest_surface_method(0, 2, 0, 0, -1, 0) - 1) < 1e-6  # Should hit at y=1

    # Case 8: y-axis cylinder
    cylinder = Cylinder("y", 2, (0, 0, 0))
    assert abs(cylinder.nearest_surface_method(2, 0, 2, -1, 0, 0) - 2) < 1e-6  # Hits at x=0


def test_edge_cases():
    # Case 1: Particle starting on the surface, moving tangentially
    cylinder = Cylinder("z", 1, (0,0,0))
    assert cylinder.nearest_surface_method(1, 0, 0, 0, 1, 0) is None  # Tangential, no intersection

    # Case 2: Particle exactly at the cylinder axis
    assert cylinder.nearest_surface_method(0, 0, 0, 1, 0, 0) == 1  # Hits at radius

    # Case 3: Large radius cylinder with diagonal movement
    large_cylinder = Cylinder("z", 1000, (0, 0, 0))
    u, v, w = direction_cosines(0, 0, 0, 500, 500, 0)

    # Solve the intersection distance
    # Equation: (500 + t * u)^2 + (500 + t * v)^2 = 1000^2
    A = u**2 + v**2
    B = 2 * (500 * u + 500 * v)
    C = 500**2 + 500**2 - 1000**2

    discriminant = B**2 - 4 * A * C
    assert discriminant >= 0, "No real intersection!"

    t1 = (-B + math.sqrt(discriminant)) / (2 * A)
    t2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # Use the positive root
    expected_t = max(t1, t2)

    # Run the test
    assert abs(large_cylinder.nearest_surface_method(500, 500, 0, u, v, w) - expected_t) < 1e-3

    # Case 4: Zero radius cylinder (degenerate case)
    degenerate_cylinder = Cylinder("z", 0, (0,0,0))
    assert degenerate_cylinder.nearest_surface_method(0, 0, 0, 1, 0, 0) is None  # No surface to hit

