import pytest
from src.medium import Region, Plane, Sphere, Cylinder, Box
from src.geometry import calculate_nearest_boundary  # Replace with the actual file name where your function resides

def test_box_region_boundary():
    # Define a box region from (0, 0, 0) to (10, 10, 10)
    planes = [
        Plane(-1, 0, 0, 0),   # x >= 0
        Plane(1, 0, 0, 10),   # x <= 10
        Plane(0, -1, 0, 0),   # y >= 0
        Plane(0, 1, 0, 10),   # y <= 10
        Plane(0, 0, -1, 0),   # z >= 0
        Plane(0, 0, 1, 10),   # z <= 10
    ]
    box_region = Region(surfaces=planes, operation="intersection")
    print(box_region)
    state = {"x": -5, "y": 5, "z": 5}
    u, v, w = 1, 0, 0  # Moving along +x direction
    point, medium, distance = calculate_nearest_boundary(state, [box_region], u, v, w)
    assert pytest.approx(distance, rel=1e-6) == 5  # Distance to the nearest boundary #Doesnt work in the negative direction, need to return absolute distance
    assert point == (0, 5, 5)
    assert medium == box_region


def test_union_of_box_and_sphere():
    # Define a box region
    planes = [
        Plane(-1, 0, 0, 0),   # x >= 0
        Plane(1, 0, 0, 10),   # x <= 10
        Plane(0, -1, 0, 0),   # y >= 0
        Plane(0, 1, 0, 10),   # y <= 10
        Plane(0, 0, -1, 0),   # z >= 0
        Plane(0, 0, 1, 10),   # z <= 10
    ]
    box_region = Region(surfaces=planes, operation="intersection")

    # Define a sphere
    sphere = Sphere(x0=15, y0=0, z0=0, radius=5)

    # Create a union of the box and the sphere
    union_region = Region(surfaces=[box_region, sphere], operation="union")

    state = {"x": 15, "y": 0, "z": 0}
    u, v, w = -1, 0, 0  # Moving towards the union region
    point, medium, distance = calculate_nearest_boundary(state, [union_region], u, v, w)
    assert pytest.approx(distance, rel=1e-6) == 5  # Distance to the sphere boundary
    assert pytest.approx(point[0], rel=1e-6) == 10  # Nearest point is on the sphere boundary
    assert medium == union_region


def test_intersection_of_box_and_sphere():
    # Define a box region
    planes = [
        Plane(-1, 0, 0, 0),   # x >= 0
        Plane(1, 0, 0, 10),   # x <= 10
        Plane(0, -1, 0, 0),   # y >= 0
        Plane(0, 1, 0, 10),   # y <= 10
        Plane(0, 0, -1, 0),   # z >= 0
        Plane(0, 0, 1, 10),   # z <= 10
    ]
    box_region = Region(surfaces=planes, operation="intersection")

    # Define a sphere intersecting with the box
    sphere = Sphere(x0=0, y0=5, z0=0, radius=5)

    # Create an intersection of the box and the sphere
    intersection_region = Region(surfaces=[box_region, sphere], operation="intersection")

    state = {"x": 0, "y": 15, "z": 0}
    u, v, w = 0, -1, 0  # Moving towards the region
    point, medium, distance = calculate_nearest_boundary(state, [intersection_region], u, v, w)
    assert pytest.approx(distance, rel=1e-6) == 5  # Distance to the nearest boundary
    assert pytest.approx(point[1], rel=1e-6) == 10  # Nearest point is within the intersection
    assert medium == intersection_region
