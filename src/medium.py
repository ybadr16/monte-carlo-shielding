# medium.py
import numpy as np
import math

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
        elif self.axis == "x":
            y_bar = y - self.y0
            z_bar = z - self.z0
            a = v**2 + w**2
            k = y_bar * v + z_bar * w
            c = y_bar**2 + z_bar**2 - self.radius**2
        elif self.axis == "y":
            x_bar = x - self.x0
            z_bar = z - self.z0
            a = u**2 + w**2
            k = x_bar * u + z_bar * w
            c = x_bar**2 + z_bar**2 - self.radius**2

        if a == 0:
            return None

        discriminant = k**2 - a * c
        if discriminant < 0:
            return None

        sqrt_d = np.sqrt(discriminant)
        d1 = (-k - sqrt_d) / a
        d2 = (-k + sqrt_d) / a

        if c < 0:
            return max(d1, d2)
        else:
            valid = [d for d in [d1, d2] if d >= 0]
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
