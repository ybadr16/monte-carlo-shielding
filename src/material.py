from .vt_calc import VelocitySampler

NEUTRON_MASS_AMU = 1.00866491588
AVOGADRO = 6.022e23  # atoms/mol


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

    @classmethod
    def mixture(cls, name, density, components, temperature=294):
        """Build a multi-isotope composition from a mass density and atom fractions.

        Parameters:
        name (str): label for the mixture (e.g. "HEU", "Stainless Steel").
        density (float): mixture mass density in g/cm³.
        components (list): one tuple per isotope, either
            (element, atomic_weight_ratio, atom_fraction) or
            (element, atomic_weight_ratio, atom_fraction, atomic_mass).
            Atom fractions need not be normalised; atomic_mass (g/mol) defaults
            to atomic_weight_ratio × neutron mass when omitted.

        Returns:
        list[Nuclide]: pass directly to Region(composition=...). Each nuclide's
        number density is density·N_A·(frac_i/Σfrac) / M̄, with M̄ the
        atom-fraction-weighted mean molar mass.
        """
        parsed = []
        frac_sum = 0.0
        for comp in components:
            if len(comp) == 4:
                element, awr, frac, amass = comp
            else:
                element, awr, frac = comp
                amass = awr * NEUTRON_MASS_AMU
            parsed.append((element, awr, frac, amass))
            frac_sum += frac

        if frac_sum <= 0:
            raise ValueError("atom fractions must sum to a positive value")

        mbar = sum((frac / frac_sum) * amass for _, _, frac, amass in parsed)
        total_N = density * AVOGADRO / mbar  # atoms/cm³

        return [
            Nuclide(element, total_N * (frac / frac_sum), awr,
                    atomic_mass=amass, temperature=temperature)
            for element, awr, frac, amass in parsed
        ]


class Nuclide:
    """One isotope within a region's material composition.

    Carries everything the transport kernel needs to treat the isotope on its
    own: the ENDF data key, its number density, the atomic-weight ratio used for
    collision kinematics, and a free-gas velocity sampler built from its mass.
    """

    def __init__(self, element, number_density, atomic_weight_ratio,
                 atomic_mass=None, temperature=294):
        """
        Parameters:
        element (str): ENDF data key (e.g. "U235", "Pb208").
        number_density (float): atoms/cm³ of this isotope in the region.
        atomic_weight_ratio (float): A = m_nuclide / m_neutron (kinematics).
        atomic_mass (float, optional): molar mass in g/mol for the free-gas
            sampler; defaults to atomic_weight_ratio × neutron mass.
        temperature (float): sampler temperature in K (default 294).
        """
        self.element = element
        self.number_density = number_density
        self.atomic_weight_ratio = atomic_weight_ratio
        self.atomic_mass = (atomic_mass if atomic_mass is not None
                            else atomic_weight_ratio * NEUTRON_MASS_AMU)
        self.kg_mass = (self.atomic_mass / AVOGADRO) * 1e-3
        self.sampler = VelocitySampler(mass=self.kg_mass, temperature=temperature)

    def __repr__(self):
        return (f"Nuclide(element={self.element}, "
                f"number_density={self.number_density:.3e} atoms/cm³, "
                f"A={self.atomic_weight_ratio})")
