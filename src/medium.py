# medium.py
import numpy as np
import math

class Region:
    def __init__(self, surfaces=None, operation="intersection", name=None, priority=0, is_void=False, element=None, composition=None):
        """
        Create a region combining surfaces using boolean operations.
        :param surfaces: List of surfaces defining the region.
        :param operation: Boolean operation ('intersection', 'union', 'complement').
        :param name: Name of the region.
        :param priority: Priority of the region (higher value means higher priority).
        :param is_void: Whether the region is a void (no material interaction).
        :param element: Single-isotope data key (uses the global A/N/sampler).
        :param composition: Optional list of Nuclide objects for a material
            mixture. When given it overrides `element` and the global A/N/sampler;
            when None the region is treated as the single isotope in `element`.
        """
        self.surfaces = surfaces if surfaces else []
        self.operation = operation
        self.name = name
        self.priority = priority
        self.is_void = is_void
        self.element = element
        self.composition = composition
        self._prim_cache = None  # flattened primitive surfaces (lazy, see geometry)

    def contains(self, x, y, z, tolerance=1e-9):
        # Short-circuit each boolean op: bail on the first decisive surface
        # rather than evaluating all of them into a list (this is the geometry
        # hot path, called once per candidate boundary crossing). The per-surface
        # test is inlined to avoid an extra call layer; a nested Region recurses,
        # a primitive is "inside" when evaluate(x,y,z) <= tolerance (the
        # tolerance treats a point sitting on the surface as inside).
        op = self.operation
        surfaces = self.surfaces

        if op == "intersection":
            for s in surfaces:
                if isinstance(s, Region):
                    if not s.contains(x, y, z, tolerance):
                        return False
                elif s.evaluate(x, y, z) > tolerance:
                    return False
            return True

        if op == "union":
            for s in surfaces:
                inside = (s.contains(x, y, z, tolerance) if isinstance(s, Region)
                          else s.evaluate(x, y, z) <= tolerance)
                if inside:
                    return True
            return False

        if op == "complement":
            s = surfaces[0]
            inside = (s.contains(x, y, z, tolerance) if isinstance(s, Region)
                      else s.evaluate(x, y, z) <= tolerance)
            return not inside

        if op == "difference":
            # A - B = A ∩ ¬B
            if len(surfaces) != 2:
                raise ValueError("Difference operation requires exactly two surfaces")
            a, b = surfaces
            a_in = (a.contains(x, y, z, tolerance) if isinstance(a, Region)
                    else a.evaluate(x, y, z) <= tolerance)
            if not a_in:
                return False
            b_in = (b.contains(x, y, z, tolerance) if isinstance(b, Region)
                    else b.evaluate(x, y, z) <= tolerance)
            return not b_in

        raise ValueError(f"Unknown operation: {op}")

    def add_surface(self, surface):
        """Add a surface to the region."""
        self.surfaces.append(surface)
        self._prim_cache = None  # invalidate the flattened-primitive cache



class Plane:
    def __init__(self, A, B, C, D, boundary_type="transmission"):
        # Normalize the normal vector (A, B, C)
        norm = np.sqrt(A**2 + B**2 + C**2)
        if np.isclose(norm, 0):
            raise ValueError("Plane normal cannot be a zero vector")

        # store as plain Python floats: the per-collision hot path evaluates
        # these millions of times and np.float64 scalar arithmetic is several
        # times slower than float for no benefit here.
        self.A = float(A / norm)
        self.B = float(B / norm)
        self.C = float(C / norm)
        self.D = float(-D / norm)  # Adjusted for normalized normal
        # "transmission" (default) lets a neutron cross and leave the geometry
        # (vacuum leakage); "reflective" specularly reflects it back in.
        self.boundary_type = boundary_type

    def evaluate(self, x, y, z):
        return self.A * x + self.B * y + self.C * z + self.D

    def nearest_surface_method(self, x, y, z, u, v, w):
        numerator = -self.D - self.A*x - self.B*y - self.C*z
        denominator = self.A*u + self.B*v + self.C*w

        # scalar abs comparisons reproduce np.isclose(.,0,atol=1e-8) at a tiny
        # fraction of the cost (np.isclose allocates arrays + runs a ufunc
        # reduce inside an errstate context — ~100x slower for a scalar).
        if abs(numerator) <= 1e-8:
            return 0.0 if abs(denominator) > 1e-8 else None

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
    def __init__(self, axis, radius, center, boundary_type="transmission"):
        self.axis = axis
        self.radius = radius
        self.x0, self.y0, self.z0 = center
        self.boundary_type = boundary_type

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
    def __init__(self, center, radius, boundary_type="transmission"):
        self.x0, self.y0, self.z0 = center
        self.radius = radius
        self.boundary_type = boundary_type

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
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max,
                 boundary_type="transmission"):
        """
        Define a box region bounded by planes at x_min, x_max, y_min, y_max, z_min, and z_max.

        `boundary_type` is applied to all six bounding planes; pass "reflective"
        for an infinite-lattice (e.g. pin-cell) calculation.
        """
        bt = boundary_type
        planes = [
            Plane(-1, 0, 0, -x_min, bt),  # x >= x_min: -x + x_min ≤ 0 → x ≥ x_min
            Plane(1, 0, 0, x_max, bt),    # x <= x_max: x - x_max ≤ 0 → x ≤ x_max
            Plane(0, -1, 0, -y_min, bt),  # y >= y_min: -y + y_min ≤ 0 → y ≥ y_min
            Plane(0, 1, 0, y_max, bt),    # y <= y_max: y - y_max ≤ 0 → y ≤ y_max
            Plane(0, 0, -1, -z_min, bt),  # z >= z_min: -z + z_min ≤ 0 → z ≥ z_min
            Plane(0, 0, 1, z_max, bt),    # z <= z_max: z - z_max ≤ 0 → z ≤ z_max
        ]
        super().__init__(surfaces=planes, operation="intersection", name="Box")
