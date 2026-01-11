from .physics import calculate_E_cm_prime
import os
import h5py
import numpy as np
from .angular_distribution import AngularDistribution

class SecondaryDistribution:
    """
    Handles Tabulated / Kalbach-Mann Secondary Energy Distributions.
    """
    def __init__(self, base_path, element, mt):
        self.base_path = base_path
        self.element = element
        self.mt = mt
        self.loaded = False

        self.incident_energy = None
        self.offsets = None
        self.distribution_data = None
        self.dist_type = None

    def load(self):
        if self.loaded: return

        file_path = os.path.join(self.base_path, f"neutron/{self.element}.h5")
        if not os.path.exists(file_path): return

        try:
            with h5py.File(file_path, 'r') as f:
                # Robust Root Detection (C0 vs C12)
                valid_roots = [k for k in f.keys() if k != "filetype"]
                if not valid_roots: return

                if f"neutron/{self.element}" in f: prefix = f"neutron/{self.element}"
                elif f"{self.element}" in f: prefix = f"{self.element}"
                else: prefix = valid_roots[0]

                base = f"{prefix}/reactions/reaction_{self.mt:03}/product_0/distribution_0"
                if base not in f: return

                grp = f[base]
                if "type" in grp.attrs:
                    self.dist_type = grp.attrs["type"].decode('utf-8')

                if "energy" in grp and "distribution" in grp:
                    self.incident_energy = grp["energy"][:]
                    self.distribution_data = grp["distribution"][:]

                    # Handle Offsets location
                    if "offsets" in grp:
                        self.offsets = grp["offsets"][:]
                    elif "offsets" in grp["energy"].attrs:
                        self.offsets = grp["energy"].attrs["offsets"]
                    elif "offsets" in grp["distribution"].attrs:
                        self.offsets = grp["distribution"].attrs["offsets"]
                    else:
                        return

                    self.loaded = True

        except Exception as e:
            print(f"Warning: Failed to load secondary dist for {self.element} MT={self.mt}: {e}")

    def sample_energy(self, E_in, rng):
        if not self.loaded: return None

        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))

        E_low = self.incident_energy[idx]
        E_high = self.incident_energy[idx+1]

        f = 0.0 if E_high == E_low else (E_in - E_low) / (E_high - E_low)

        E_out_low = self._sample_from_block(idx, rng)
        E_out_high = self._sample_from_block(idx+1, rng)

        return (1 - f) * E_out_low + f * E_out_high

    def _sample_from_block(self, idx, rng):
        start = self.offsets[idx]
        if idx == len(self.incident_energy) - 1:
            end = self.distribution_data.shape[1]
        else:
            end = self.offsets[idx + 1]

        if start >= end: return 0.0

        block_E_out = self.distribution_data[0, start:end]
        block_CDF   = self.distribution_data[2, start:end]

        r = rng.random()
        return np.interp(r, block_CDF, block_E_out)


class CrossSectionReader:
    def __init__(self, base_path: str, temperature: str = "294K"):
        """
        Initialize with a specific temperature.
        Default is 294K (Room Temp).
        """
        self.base_path = base_path
        self.temperature = temperature
        self._cache = {}
        self._available_inelastic_mts = {}
        self.angular_dists = {}
        self.secondary_dists = {}

    def _load_data_to_cache(self, element: str, mt: int):
        if not element.isalnum(): raise ValueError("Invalid element")

        mt_str = f"{mt:03}"
        file_path = os.path.join(self.base_path, f"neutron/{element}.h5")

        if not os.path.exists(file_path):
            self._cache[(element, mt)] = None
            return

        try:
            with h5py.File(file_path, 'r') as f:
                # --- 1. Root Detection ---
                valid_roots = [k for k in f.keys() if k != "filetype"]
                if not valid_roots:
                    self._cache[(element, mt)] = None; return

                if f"neutron/{element}" in f: root = f"neutron/{element}"
                elif f"{element}" in f: root = f"{element}"
                else: root = valid_roots[0]

                reaction_group = f"{root}/reactions/reaction_{mt_str}"
                if reaction_group not in f:
                    self._cache[(element, mt)] = None; return

                # --- 2. Explicit Temperature Selection ---
                # We look for exactly what the user asked for.
                # Only exception: 294K -> 293.6K alias mapping (common in OpenMC)

                available_temps = list(f[reaction_group].keys())
                target_temp = self.temperature

                if target_temp not in available_temps:
                    # Check for alias ONLY if asking for room temp
                    if target_temp == "294K" and "293.6K" in available_temps:
                        target_temp = "293.6K"
                    else:
                        # User asked for X, X is not there. Fail gracefully.
                        self._cache[(element, mt)] = None
                        return

                reaction_temp = f"{reaction_group}/{target_temp}"

                # --- 3. Energy Path ---
                # Standard: {root}/energy/{temp}
                energy_path = f"{root}/energy/{target_temp}"

                # Fallback: Sometimes energy grid is named differently or shared
                if energy_path not in f and f"{root}/energy" in f:
                    # If exact temp missing in energy, try grabbing the first one ONLY if it looks safe
                    # (This is rare, usually energy grids match temp names)
                    e_keys = list(f[f"{root}/energy"].keys())
                    if target_temp in e_keys:
                        energy_path = f"{root}/energy/{target_temp}"
                    elif len(e_keys) > 0:
                        energy_path = f"{root}/energy/{e_keys[0]}"

                if energy_path not in f: raise KeyError("Energy path missing")

                # --- 4. Load Data ---
                energy_data = f[energy_path][:]
                q_value = f[reaction_group].attrs.get('Q_value', 0.0)

                if f"{reaction_temp}/xs" not in f:
                    self._cache[(element, mt)] = None; return

                xs_dataset = f[f"{reaction_temp}/xs"]
                xs_data = xs_dataset[:]
                threshold_idx = xs_dataset.attrs.get('threshold_idx', 0)

                xs_full = np.zeros_like(energy_data)
                if threshold_idx + len(xs_data) <= len(xs_full):
                    xs_full[threshold_idx:threshold_idx + len(xs_data)] = xs_data
                else:
                    available = len(xs_full) - threshold_idx
                    xs_full[threshold_idx:] = xs_data[:available]

                self._cache[(element, mt)] = {
                    'energy': energy_data,
                    'xs': xs_full,
                    'threshold_energy': energy_data[threshold_idx] if threshold_idx < len(energy_data) else 0.0,
                    'q_value': q_value
                }

        except Exception:
            self._cache[(element, mt)] = None

    def _scan_inelastic_mts(self, element: str):
        if element in self._available_inelastic_mts: return

        file_path = os.path.join(self.base_path, f"neutron/{element}.h5")
        found = []
        check_list = [16, 17] + list(range(51, 92))

        try:
            with h5py.File(file_path, 'r') as f:
                valid_roots = [k for k in f.keys() if k != "filetype"]
                if valid_roots:
                    if f"neutron/{element}" in f: root = f"neutron/{element}"
                    elif f"{element}" in f: root = f"{element}"
                    else: root = valid_roots[0]

                    root_reactions = f"{root}/reactions"
                    if root_reactions in f:
                        for mt in check_list:
                            if f"reaction_{mt:03}" in f[root_reactions]: found.append(mt)
        except: pass
        self._available_inelastic_mts[element] = found

    def get_cross_section_data(self, element: str, mt: int):
        cache_key = (element, mt)
        if cache_key not in self._cache:
            self._load_data_to_cache(element, mt)
        return self._cache[cache_key]

    def get_cross_section(self, element: str, mt: int, energy: float) -> float:
        data = self.get_cross_section_data(element, mt)
        if data is None or energy < data['threshold_energy']: return 0.0
        return np.interp(energy, data['energy'], data['xs'])

    def calculate_macroscopic_xs(self, mic, N):
        return mic * 1e-24 * N if mic > 0 else 0.0

    def get_inelastic_components(self, element, energy, N):
        if element not in self._available_inelastic_mts: self._scan_inelastic_mts(element)
        total = 0.0
        comps = []
        for mt in self._available_inelastic_mts[element]:
            mic = self.get_cross_section(element, mt, energy)
            if mic > 0:
                mac = self.calculate_macroscopic_xs(mic, N)
                if mac > 0:
                    data = self.get_cross_section_data(element, mt)
                    q = data['q_value'] if data else 0.0
                    total += mac
                    comps.append((mt, mac, q))
        return total, comps

    def get_elastic_mu(self, element, energy, rng):
        if element not in self.angular_dists:
            dist = AngularDistribution(self.base_path, element, mt=2)
            dist.load()
            self.angular_dists[element] = dist
        return self.angular_dists[element].sample_mu(energy, rng)

    def get_secondary_energy(self, element, mt, energy, rng):
        key = (element, mt)
        if key not in self.secondary_dists:
            dist = SecondaryDistribution(self.base_path, element, mt)
            dist.load()
            self.secondary_dists[key] = dist
        return self.secondary_dists[key].sample_energy(energy, rng)

    def get_cross_sections(self, element, energy, sampler, N, A):
        """
        Get energy-dependent macroscopic cross sections.
        """
        if energy < 10.0:
            lookup_energy = calculate_E_cm_prime(energy, A, sampler)
        else:
            lookup_energy = energy

        mic_el = self.get_cross_section(element, 2, lookup_energy)
        Sigma_el = self.calculate_macroscopic_xs(mic_el, N)

        mic_cap = self.get_cross_section(element, 102, energy)
        Sigma_cap = self.calculate_macroscopic_xs(mic_cap, N)

        mic_fis = self.get_cross_section(element, 18, energy)
        Sigma_fis = self.calculate_macroscopic_xs(mic_fis, N)

        Sigma_in, _ = self.get_inelastic_components(element, energy, N)

        Sigma_total = Sigma_el + Sigma_in + Sigma_cap + Sigma_fis

        return Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_total
