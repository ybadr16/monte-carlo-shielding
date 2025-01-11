from physics import calculate_E_cm_prime
import os
import h5py
import numpy as np

class CrossSectionReader:
    def __init__(self, base_path: str):
        """
        Initialize the CrossSectionReader with the base path to the data files.
        :param base_path: Base directory where HDF5 files are located.
        """
        self.base_path = base_path

    def get_cross_section(self, element: str, mt: int, energy: float) -> float:
        """
        Get the cross-section for a specific nuclide, reaction, and energy.
        :param element: The nuclide in the format "U235" or "Pb208".
        :param mt: The reaction MT number (e.g., 18 for fission).
        :param energy: The energy in eV for which the cross-section is required.
        :return: The interpolated cross-section at the specified energy.
        """
        # Validate inputs
        if not element.isalnum():
            raise ValueError("Invalid element format. Use alphanumeric characters (e.g., U235, Pb208).")
        if not (1 <= mt <= 999):
            raise ValueError("MT number must be between 1 and 999.")

        # Format MT number as a 3-digit string
        mt_str = f"{mt:03}"

        # Construct file path
        file_path = os.path.join(self.base_path, f"neutron/{element}.h5")
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"HDF5 file for {element} not found at {file_path}.")

        # Define HDF5 dataset paths
        reaction_group_path = f"{element}/reactions/reaction_{mt_str}/294K"
        energy_path = f"{element}/energy/294K"

        # Read the HDF5 file
        try:
            with h5py.File(file_path, 'r') as f:
                # Load energy data
                if energy_path not in f:
                    raise KeyError(f"Energy data path '{energy_path}' not found in HDF5 file.")
                energy_data = f[energy_path][:]

                # Load cross-section data
                if reaction_group_path not in f:
                    raise KeyError(f"Reaction group path '{reaction_group_path}' not found in HDF5 file.")
                if f"{reaction_group_path}/xs" not in f:
                    raise KeyError(f"Cross-section data not found at '{reaction_group_path}/xs'.")
                xs_data = f[f"{reaction_group_path}/xs"][:]
                threshold_idx = f[f"{reaction_group_path}/xs"].attrs.get('threshold_idx', 0)

                # Validate threshold index
                if not (0 <= threshold_idx < len(energy_data)):
                    raise ValueError("Invalid threshold index in the HDF5 file.")

                # Construct the full cross-section array
                xs_full = np.zeros_like(energy_data)
                xs_full[threshold_idx:threshold_idx + len(xs_data)] = xs_data

                # Interpolate cross-section at the requested energy
                if energy < energy_data[threshold_idx]:
                    return 0.0  # Below threshold energy
                cross_section = np.interp(energy, energy_data, xs_full)
                return cross_section

        except (OSError, KeyError, ValueError) as e:
            raise RuntimeError(f"Error while reading HDF5 file: {e}") from e

    def calculate_macroscopic_xs(self, microscopic_xs: float, number_density: float) -> float:
        """
        Calculate macroscopic cross section from microscopic cross section and number density.

        :param microscopic_xs: Microscopic cross section in barns (10⁻²⁴ cm²)
        :param number_density: Number density in atoms/cm³
        :return: Macroscopic cross section in cm⁻¹
        :raises ValueError: If inputs are negative
        """
        if microscopic_xs < 0:
            raise ValueError("Microscopic cross section cannot be negative")
        if number_density < 0:
            raise ValueError("Number density cannot be negative")

        microscopic_xs_cm2 = microscopic_xs * 1e-24  # Convert barns to cm²
        return microscopic_xs_cm2 * number_density

    def get_macroscopic_xs(self, element: str, mt: int, energy: float, number_density: float) -> float:
        """
        Convenience method to get macroscopic cross section directly from element, MT number, and energy.

        :param element: The nuclide in the format "U235" or "Pb208"
        :param mt: The reaction MT number
        :param energy: The energy in eV
        :param number_density: Number density in atoms/cm³
        :return: Macroscopic cross section in cm⁻¹
        """
        microscopic_xs = self.get_cross_section(element, mt, energy)
        return self.calculate_macroscopic_xs(microscopic_xs, number_density)

    def get_cross_sections(self, element, energy, sampler, number_density):
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
        Sigma_s = self.get_macroscopic_xs(element, 2, energy_cm, number_density)    # Scattering
        Sigma_a = self.get_macroscopic_xs(element, 102, energy, number_density)     # Radiative capture

        try:
            Sigma_f = self.get_macroscopic_xs(element, 18, energy, number_density)  # Fission
        except RuntimeError:
            Sigma_f = 0  # If fission isn't available, set it to zero

        Sigma_t = Sigma_s + Sigma_a + Sigma_f

        return Sigma_s, Sigma_a, Sigma_f, Sigma_t

