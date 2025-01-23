class Material:
    def __init__(self, name, density, atomic_mass, atomic_weight_ratio=None):
        """
        Initialize a Material object.

        Parameters:
        name (str): Name of the material (e.g., "Lead", "Cadmium").
        density (float): Material density in g/cm³.
        atomic_mass (float): Atomic mass in g/mol.
        atomic_weight_ratio (float, optional): Atomic weight ratio. Default is None.
        """
        self.name = name
        self.density = density
        self.atomic_mass = atomic_mass
        self.atomic_weight_ratio = atomic_weight_ratio
        self.number_density = self.calculate_number_density()
        self.kg_mass = self.calculate_atomic_mass_kg()

    def calculate_number_density(self):
        """
        Calculate number density (atoms/cm³) from material density and atomic mass.

        Returns:
        float: Number density in atoms/cm³.
        """
        avogadro = 6.022e23  # Avogadro's number (atoms/mol)
        return (self.density * avogadro) / self.atomic_mass

    def calculate_atomic_mass_kg(self):
        """
        Calculate the mass of a single atom in kilograms.

        Returns:
        float: Mass of a single atom in kilograms.
        """
        avogadro = 6.022e23  # Avogadro's number (atoms/mol)
        return (self.atomic_mass / avogadro) * 1e-3  # Convert g to kg

    def __repr__(self):
        """
        String representation of the Material object for debugging.
        """
        return (
            f"Material(name={self.name}, density={self.density} g/cm³, "
            f"atomic_mass={self.atomic_mass} g/mol, "
            f"atomic_weight_ratio={self.atomic_weight_ratio}, "
            f"number_density={self.number_density:.2e} atoms/cm³)"
        )
