# cross_sections.py
from physics import calculate_E_cm_prime

def get_cross_sections(reader, element, energy, sampler):
    """
    Get energy-dependent cross sections for a given element and energy.

    Args:
        reader: CrossSectionReader instance
        element: Name of the element (e.g., "U235")
        energy: Neutron energy (eV)

    Returns:
        sigma_s: Scattering cross-section
        sigma_a: Absorption cross-section
        sigma_f: Fission cross-section (0 if not fissionable)
        Sigma_t: Total cross-section
    """
    # Get cross-sections from the reader. If fission isn't available, set it to 0.
    energy_cm = calculate_E_cm_prime(energy, 2.5, sampler)#2.5 is A for now
    sigma_s = reader.get_cross_section(element, 2, energy_cm)  # Scattering
    sigma_a = reader.get_cross_section(element, 102, energy)  # Radiative capture
    try:
        sigma_f = reader.get_cross_section(element, 18, energy)  # Fission
    except RuntimeError:
        sigma_f = 0  # If fission isn't available, set it to zero.

    Sigma_t = sigma_s + sigma_a + sigma_f
    return sigma_s, sigma_a, sigma_f, Sigma_t
