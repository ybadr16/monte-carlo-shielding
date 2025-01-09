# cross_sections.py
from physics import calculate_E_cm_prime

def get_cross_sections(reader, element, energy, sampler, number_density):
    """
    Get energy-dependent macroscopic cross sections for a given element and energy.

    Args:
        reader: CrossSectionReader instance
        element: Name of the element (e.g., "U235")
        energy: Neutron energy (eV)
        sampler: Sampler instance for CM energy calculation
        number_density: Number density in atoms/cm³

    Returns:
        Sigma_s: Macroscopic scattering cross-section (cm⁻¹)
        Sigma_a: Macroscopic absorption cross-section (cm⁻¹)
        Sigma_f: Macroscopic fission cross-section (cm⁻¹, 0 if not fissionable)
        Sigma_t: Total macroscopic cross-section (cm⁻¹)
    """
    # Get microscopic cross-sections and convert to macroscopic
    energy_cm = calculate_E_cm_prime(energy, 2.5, sampler)  # 2.5 is A for now

    # Calculate each macroscopic cross section directly
    Sigma_s = reader.get_macroscopic_xs(element, 2, energy_cm, number_density)    # Scattering
    Sigma_a = reader.get_macroscopic_xs(element, 102, energy, number_density)     # Radiative capture

    try:
        Sigma_f = reader.get_macroscopic_xs(element, 18, energy, number_density)  # Fission
    except RuntimeError:
        Sigma_f = 0  # If fission isn't available, set it to zero

    Sigma_t = Sigma_s + Sigma_a + Sigma_f

    return Sigma_s, Sigma_a, Sigma_f, Sigma_t
