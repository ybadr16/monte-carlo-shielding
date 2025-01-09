def calculate_number_density(density, atomic_mass):
    """
    Calculate number density (atoms/cm³) from material density and atomic mass

    Parameters:
    density (float): Material density in g/cm³
    atomic_mass (float): Atomic mass in g/mol

    Returns:
    float: Number density in atoms/cm³
    """
    avogadro = 6.022e23  # Avogadro's number (atoms/mol)
    N = (density * avogadro) / atomic_mass
    return N
