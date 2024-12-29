import h5py
import numpy as np
import os

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
        reaction_group_path = f"{element}/reactions/reaction_{mt_str}/1200K"
        energy_path = f"{element}/energy/1200K"

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

'''
# Example usage
if __name__ == "__main__":
    base_path = "./"
    reader = CrossSectionReader(base_path)

    try:
        element = "U235"
        mt = 18  # Example reaction MT number
        energy = 0.01  # Example energy in eV
        cross_section = reader.get_cross_section(element, mt, energy)
        print(f"Cross-section for {element} (MT={mt}) at {energy} eV: {cross_section:.4e} barns")
    except Exception as e:
        print(f"An error occurred: {e}")
'''
