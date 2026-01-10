from .physics import calculate_E_cm_prime
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
        # Cache structure: {(element, mt): {'energy': np.array, 'xs': np.array, 'q_value': float}}
        self._cache = {}

        # Store a list of available inelastic MTs per element so we don't check empty ones repeatedly
        self._available_inelastic_mts = {}

    def _load_data_to_cache(self, element: str, mt: int):
        """
        Internal method to load data from HDF5 into memory cache.
        """
        # Validate inputs
        if not element.isalnum():
            raise ValueError("Invalid element format. Use alphanumeric characters (e.g., U235, Pb208).")

        mt_str = f"{mt:03}"
        file_path = os.path.join(self.base_path, f"neutron/{element}.h5")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"HDF5 file for {element} not found at {file_path}.")

        reaction_group_path = f"{element}/reactions/reaction_{mt_str}"
        # Some HDF5 structures have the temperature directly under reaction, others don't.
        # OpenMC format usually: reaction_XXX/294K/xs
        # We assume 294K (room temp) for now.
        temp_group = "294K"
        reaction_temp_path = f"{reaction_group_path}/{temp_group}"
        energy_path = f"{element}/energy/{temp_group}"

        try:
            with h5py.File(file_path, 'r') as f:
                # Load energy data
                if energy_path not in f:
                    raise KeyError(f"Energy data path '{energy_path}' not found in HDF5 file.")
                energy_data = f[energy_path][:]

                # Check if reaction exists
                if reaction_group_path not in f:
                    # Reaction doesn't exist (e.g., MT=18 for Lead).
                    # Cache a dummy zero entry to prevent repeated file opens.
                    self._cache[(element, mt)] = None
                    return

                # Load Q-Value (Energy released/consumed)
                # In OpenMC HDF5, Q-value is often an attribute of the reaction group
                q_value = f[reaction_group_path].attrs.get('Q_value', 0.0)

                # Load cross-section data
                if f"{reaction_temp_path}/xs" not in f:
                    # Exists but no XS data? Cache None.
                    self._cache[(element, mt)] = None
                    return

                xs_dataset = f[f"{reaction_temp_path}/xs"]
                xs_data = xs_dataset[:]
                threshold_idx = xs_dataset.attrs.get('threshold_idx', 0)

                # Construct the full cross-section array in memory
                xs_full = np.zeros_like(energy_data)

                # Safety bounds check
                if threshold_idx + len(xs_data) <= len(xs_full):
                    xs_full[threshold_idx:threshold_idx + len(xs_data)] = xs_data
                else:
                    # Handle rare edge cases where index exceeds bounds
                    available_len = len(xs_full) - threshold_idx
                    xs_full[threshold_idx:] = xs_data[:available_len]

                # Store in cache
                self._cache[(element, mt)] = {
                    'energy': energy_data,
                    'xs': xs_full,
                    'threshold_energy': energy_data[threshold_idx] if threshold_idx < len(energy_data) else 0.0,
                    'q_value': q_value
                }

        except (OSError, KeyError, ValueError) as e:
            raise RuntimeError(f"Error while reading HDF5 file: {e}") from e

    def _scan_inelastic_mts(self, element: str):
        """
        Helper to scan the HDF5 file ONCE to find which Inelastic MTs (51-91) actually exist.
        """
        if element in self._available_inelastic_mts:
            return

        file_path = os.path.join(self.base_path, f"neutron/{element}.h5")
        found_mts = []
        try:
            with h5py.File(file_path, 'r') as f:
                reactions_root = f"{element}/reactions"
                if reactions_root in f:
                    # Check standard discrete levels 51 through 90
                    for mt in range(51, 92):
                        mt_str = f"reaction_{mt:03}"
                        if mt_str in f[reactions_root]:
                            found_mts.append(mt)
        except:
            pass # Fail silently, will default to empty list

        self._available_inelastic_mts[element] = found_mts

    def get_cross_section_data(self, element: str, mt: int):
        """
        Returns the full data dict (energy, xs, q_value) for a reaction.
        """
        cache_key = (element, mt)
        if cache_key not in self._cache:
            self._load_data_to_cache(element, mt)

        return self._cache[cache_key]

    def get_cross_section(self, element: str, mt: int, energy: float) -> float:
        """
        Get the scalar cross-section value.
        """
        data = self.get_cross_section_data(element, mt)

        if data is None:
            return 0.0

        if energy < data['threshold_energy']:
            return 0.0

        return np.interp(energy, data['energy'], data['xs'])

    def calculate_macroscopic_xs(self, microscopic_xs: float, number_density: float) -> float:
        if microscopic_xs < 0: return 0.0
        return microscopic_xs * 1e-24 * number_density

    def get_inelastic_components(self, element, energy, number_density):
        """
        Returns a breakdown of inelastic scattering components for a given energy.
        Used to determine WHICH inelastic level occurred.

        Returns:
            total_inelastic_xs (float),
            components (list of tuples): [(mt, probability, q_value), ...]
        """
        # Ensure we know which MTs exist
        if element not in self._available_inelastic_mts:
            self._scan_inelastic_mts(element)

        total_xs = 0.0
        components = []

        # Iterate only through known existing MTs
        for mt in self._available_inelastic_mts[element]:
            mic_xs = self.get_cross_section(element, mt, energy)
            if mic_xs > 0:
                mac_xs = self.calculate_macroscopic_xs(mic_xs, number_density)
                if mac_xs > 0:
                    data = self.get_cross_section_data(element, mt)
                    q_val = data['q_value'] if data else 0.0
                    total_xs += mac_xs
                    components.append((mt, mac_xs, q_val))

        return total_xs, components

    def get_cross_sections(self, element, energy, sampler, number_density):
        """
        Get energy-dependent macroscopic cross sections.

        Returns:
            Sigma_elastic:  MT 2
            Sigma_inelastic: Sum of MT 51-91
            Sigma_capture:  MT 102
            Sigma_fission:  MT 18
            Sigma_total:    Sum of above
        """
        # Effective energy for thermal treatment
        energy_cm = calculate_E_cm_prime(energy, 2.5, sampler)

        # 1. Elastic (MT 2)
        mic_el = self.get_cross_section(element, 2, energy_cm)
        Sigma_el = self.calculate_macroscopic_xs(mic_el, number_density)

        # 2. Capture (MT 102)
        mic_cap = self.get_cross_section(element, 102, energy)
        Sigma_cap = self.calculate_macroscopic_xs(mic_cap, number_density)

        # 3. Fission (MT 18)
        mic_fis = self.get_cross_section(element, 18, energy)
        Sigma_fis = self.calculate_macroscopic_xs(mic_fis, number_density)

        # 4. Inelastic (MT 51-91)
        # We just want the sum here for transport probability.
        # Detailed breakdown is only needed if a collision happens.
        Sigma_in, _ = self.get_inelastic_components(element, energy, number_density)

        Sigma_total = Sigma_el + Sigma_in + Sigma_cap + Sigma_fis

        return Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_total
