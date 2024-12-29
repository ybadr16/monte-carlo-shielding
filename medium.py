# medium.py

class Medium:
    def __init__(self, name, x_bounds, y_bounds, z_bounds, element=None, is_void=False, priority=0):
        """
        Represents a medium with specific boundaries and cross-section properties.

        Args:
            name: Name of the medium (e.g., "U235", "Void").
            x_bounds: Tuple (x_min, x_max) defining x boundaries.
            y_bounds: Tuple (y_min, y_max) defining y boundaries.
            z_bounds: Tuple (z_min, z_max) defining z boundaries.
            element: Name of the element associated with the medium (if applicable).
            is_void: Boolean indicating if the medium is a void (no interactions).
            priority: Integer priority of the medium. Higher priority = more specific medium.
        """
        self.name = name
        self.x_bounds = x_bounds
        self.y_bounds = y_bounds
        self.z_bounds = z_bounds
        self.element = element
        self.is_void = is_void
        self.priority = priority

    def contains(self, x, y, z):
        """Check if a point (x, y, z) is inside the medium."""
        x_min, x_max = self.x_bounds
        y_min, y_max = self.y_bounds
        z_min, z_max = self.z_bounds
        return x_min <= x <= x_max and y_min <= y <= y_max and z_min <= z <= z_max
