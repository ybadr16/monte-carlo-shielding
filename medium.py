# medium.py
import numpy as np
import math

class Surface:
    def __init__(self):
        pass

    def evaluate(self, x, y, z):
        """
        Evaluate the surface equation at a given point.
        :return: A value indicating whether the point is inside (negative), outside (positive), or on (zero) the surface.
        """
        raise NotImplementedError("Subclasses must implement this method.")


class Region:
    def __init__(self, surfaces=None, operation="intersection"):
        """
        Create a region combining surfaces using boolean operations.
        :param surfaces: List of surfaces defining the region.
        :param operation: Boolean operation ('intersection', 'union', 'complement').
        """
        self.surfaces = surfaces if surfaces else []
        self.operation = operation

    def contains(self, x, y, z):
        """
        Check if a point is within the region based on the boolean operation.
        :return: True if the point is inside the region, False otherwise.
        """
        evaluations = []
        for surface in self.surfaces:
            if isinstance(surface, Region):
                # Recursively check nested regions
                evaluations.append(surface.contains(x, y, z))
            else:
                evaluations.append(surface.evaluate(x, y, z) <= 0)

        if self.operation == "intersection":
            return all(evaluations)
        elif self.operation == "union":
            return any(evaluations)
        elif self.operation == "complement":
            # Complement of the first surface
            return not evaluations[0]
        else:
            raise ValueError(f"Unknown operation: {self.operation}")

    def add_surface(self, surface):
        """Add a surface to the region."""
        self.surfaces.append(surface)


class Plane(Surface):
    def __init__(self, A, B, C, D):
        """
        Initialize a generic plane Ax + By + Cz = D.
        :param A: Coefficient for x
        :param B: Coefficient for y
        :param C: Coefficient for z
        :param D: Offset from the origin
        """
        super().__init__()
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def evaluate(self, x, y, z):
        return self.A * x + self.B * y + self.C * z - self.D

    def nearest_surface_method(self, x, y, z, u, v, w):
        """
        Compute the distance to the plane from a point moving in a given direction.
        :param x, y, z: Coordinates of the starting point
        :param u, v, w: Direction vector components
        :return: Distance to the plane (or None if parallel)
        """
        numerator = self.D - (self.A * x + self.B * y + self.C * z)
        denominator = self.A * u + self.B * v + self.C * w

        if denominator == 0:
            return None  # Parallel to the plane

        return numerator / denominator

class Sphere(Surface):
    def __init__(self, x0, y0, z0, radius):
        super().__init__()
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.radius = radius

    def evaluate(self, x, y, z):
        return (x - self.x0) ** 2 + (y - self.y0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2

    def nearest_surface_method(self, x, y, z, u, v, w):
        """
        Compute the nearest distance to the surface of the sphere along the given direction.
        :param x, y, z: Starting coordinates of the particle.
        :param u, v, w: Direction vector components.
        :return: Distance to the nearest surface (or None if no intersection).
        """
        x_bar = x - self.x0
        y_bar = y - self.y0
        z_bar = z - self.z0

        k = x_bar * u + y_bar * v + z_bar * w
        c = x_bar ** 2 + y_bar ** 2 + z_bar ** 2 - self.radius ** 2

        discriminant = k ** 2 - c

        if discriminant < 0:
            return None  # No intersection

        sqrt_discriminant = discriminant ** 0.5
        d1 = -k - sqrt_discriminant
        d2 = -k + sqrt_discriminant

        if c < 0:
            # Particle is inside the sphere
            return max(d1, d2)  # Return the positive distance
        else:
            # Particle is outside the sphere
            if d1 > 0:
                return d1  # Return the smaller positive distance
            elif d2 > 0:
                return d2
            else:
                return None  # Both distances are negative

class Cylinder(Surface):
    def __init__(self, radius, x0 = None, y0 = None, z0 = None, axis="z"):
        super().__init__()
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.radius = radius
        self.axis = axis.lower()
        if self.axis not in {"x", "y", "z"}:
            raise ValueError("Axis must be 'x', 'y', or 'z'.")

    def evaluate(self, x, y, z):
        if self.axis == "z":
            return (x - self.x0) ** 2 + (y - self.y0) ** 2 - self.radius ** 2
        elif self.axis == "x":
            return (y - self.y0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2
        elif self.axis == "y":
            return (x - self.x0) ** 2 + (z - self.z0) ** 2 - self.radius ** 2

    def nearest_surface_method(self, x, y, z, u, v, w):
        """
        Compute the nearest distance to the surface of the cylinder along the given direction.
        :param x, y, z: Starting coordinates of the particle.
        :param u, v, w: Direction vector components.
        :return: Distance to the nearest surface (or None if no intersection).
        """
        if self.axis == "z":
            y_bar = y - self.y0
            x_bar = x - self.x0
            a = v ** 2 + u ** 2
            k = y_bar * v + x_bar * u
            c = y_bar ** 2 + x_bar ** 2 - self.radius ** 2
        elif self.axis == "x":
            z_bar = z - self.z0
            y_bar = y - self.y0
            a = w ** 2 + v ** 2
            k = z_bar * w + y_bar * v
            c = z_bar ** 2 + y_bar ** 2 - self.radius ** 2
        elif self.axis == "y":
            x_bar = x - self.x0
            z_bar = z - self.z0
            a = u ** 2 + w ** 2
            k = x_bar * u + z_bar * w
            c = x_bar ** 2 + z_bar ** 2 - self.radius ** 2

        if a == 0:
            return None  # Parallel to the cylinder

        discriminant = k ** 2 - a * c

        if discriminant < 0:
            return None  # No intersection

        sqrt_discriminant = discriminant ** 0.5
        d1 = (-k - sqrt_discriminant) / a
        d2 = (-k + sqrt_discriminant) / a

        if c < 0:
            # Particle is inside the cylinder
            return max(d1, d2)  # Return the positive distance
        else:
            # Particle is outside the cylinder
            if d1 > 0:
                return d1  # Return the smaller positive distance
            elif d2 > 0:
                return d2
            else:
                return None  # Both distances are negative



class Box(Region):
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        """
        Define a box region bounded by planes at x_min, x_max, y_min, y_max, z_min, and z_max.
        :param x_min: Minimum x-coordinate
        :param x_max: Maximum x-coordinate
        :param y_min: Minimum y-coordinate
        :param y_max: Maximum y-coordinate
        :param z_min: Minimum z-coordinate
        :param z_max: Maximum z-coordinate
        """
        # Create the six planes defining the box
        planes = [
            Plane(-1, 0, 0, x_min),  # x >= x_min
            Plane(1, 0, 0, x_max),  # x <= x_max
            Plane(0, -1, 0, y_min),  # y >= y_min
            Plane(0, 1, 0, y_max),  # y <= y_max
            Plane(0, 0, -1, z_min),  # z >= z_min
            Plane(0, 0, 1, z_max),  # z <= z_max
        ]
        # Initialize the Region with these planes and the intersection operation
        super().__init__(surfaces=planes, operation="intersection")



# Define six planes for a box from (0, 0, 0) to (10, 10, 10)
planes = [
    Plane(-1, 0, 0, 0),   # x >= 0
    Plane(1, 0, 0, 10),   # x <= 10
    Plane(0, -1, 0, 0),   # y >= 0
    Plane(0, 1, 0, 10),   # y <= 10
    Plane(0, 0, -1, 0),   # z >= 0
    Plane(0, 0, 1, 10),   # z <= 10
]


# Combine planes into a box region
box_region = Region(surfaces=planes, operation="intersection")

# Check if points are inside the box
point_inside = (5, 5, 5)
point_outside = (12, 0, 0)

print("Point", point_inside, "is inside the box:", box_region.contains(*point_inside))
print("Point", point_outside, "is inside the box:", box_region.contains(*point_outside))


# Define a sphere
sphere = Sphere(10, 0, 0, 3)  # Sphere centered at (5, 5, 5) with radius 3

# Create a union of the box and the sphere
union_region = Region(surfaces=[box_region, sphere], operation="union")

print("Point", point_inside, "is inside the union region:", union_region.contains(*point_inside))
print("Point", point_outside, "is inside the union region:", union_region.contains(*point_outside))
