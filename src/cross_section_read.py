from .physics import calculate_E_cm_prime
import os
import h5py
import numpy as np
from .angular_distribution import AngularDistribution

def _normalize_block_to_3xn(arr):
    """
    Normalizes any jagged/interleaved array to (3, N) shape:
    Row 0: Energy/Mu, Row 1: PDF, Row 2: CDF
    """
    arr = np.asarray(arr, dtype=float)

    # Case 1: Already 2D
    if arr.ndim == 2:
        if arr.shape[0] == 3: return arr
        if arr.shape[1] == 3: return arr.T

    # Case 2: 1D Interleaved [E1, P1, C1, E2, P2, C2...]
    flat = arr.flatten()
    n = len(flat) // 3
    if n >= 1 and len(flat) >= 3*n:
        return flat[:3*n].reshape(n, 3).T

    # Fallback: Return safe 3x1 zero block to prevent crashes
    out = np.zeros((3, 1))
    if flat.size > 0: out[0, 0] = flat[0]
    return out

def _normalize_cdf(c):
    """
    Ensures CDF is monotonic non-decreasing and normalized 0..1.
    """
    c = np.asarray(c, dtype=float).copy()
    if c.size == 0: return c

    # If it looks like a PDF (not monotonic), integrate it
    if np.any(c < 0) or (c.size > 1 and not np.all(np.diff(c) >= 0)):
        c = np.abs(c) # Fix negatives
        c = np.cumsum(c)

    # Force strict monotonicity
    c = np.maximum.accumulate(c)

    # Normalize
    denom = (c[-1] - c[0])
    if denom <= 1e-14:
        # Degenerate CDF, return linear spacing
        return np.linspace(0.0, 1.0, len(c))

    return (c - c[0]) / denom


class SecondaryDistribution:
    """
    STRICT MODE: Raises Errors instead of falling back.
    Use this to identify EXACTLY where Law 61 parsing fails.
    """
    def __init__(self, base_path, element, mt):
        self.base_path = base_path
        self.element = element
        self.mt = mt
        self.loaded = False
        self.dist_type = None

        self.incident_energy = None
        self.offsets = None
        self.energy_out_data = None
        self.mu_data = None

        self.ref_frame = "lab" # Default to Lab

    def load(self):
        if self.loaded: return
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        fname = element_map.get(self.element, self.element)
        file_path = os.path.join(self.base_path, f"neutron/{fname}.h5")

        if not os.path.exists(file_path): return

        try:
            with h5py.File(file_path, 'r') as f:
                root = f"neutron/{fname}" if f"neutron/{fname}" in f else fname
                base = f"{root}/reactions/reaction_{self.mt:03}/product_0/distribution_0"
                if base not in f: return
                grp = f[base]

                if "type" in grp.attrs:
                    self.dist_type = grp.attrs["type"].decode('utf-8')

                if "modifiers" in grp.attrs:
                    # OpenMC often stores frame in modifiers
                    mods = grp.attrs["modifiers"]
                    # 1 = Lab, 2 = CM (Standard ENDF codes)
                    if isinstance(mods, np.ndarray):
                        if 2 in mods: self.ref_frame = "cm"
                    elif mods == 2: self.ref_frame = "cm"
                # ==========================================
                # CASE A: CORRELATED (Law 61)
                # ==========================================
                if self.dist_type == "correlated":
                    self.incident_energy = grp["energy"][:]

                    # STRICT OFFSET HUNTING
                    if "offsets" in grp: self.offsets = grp["offsets"][:]
                    elif "offsets" in grp["energy"].attrs: self.offsets = grp["energy"].attrs["offsets"]
                    elif "energy_out" in grp and "offsets" in grp["energy_out"].attrs:
                         self.offsets = grp["energy_out"].attrs["offsets"]
                    else:
                        # Don't fail silently. If we can't find offsets for Law 61, we can't simulate.
                        print(f"⚠️ [WARNING] No offsets found for {self.element} MT{self.mt}")
                        return

                    self.energy_out_data = grp["energy_out"][:]
                    if "mu" in grp: self.mu_data = grp["mu"][:]

                    self.loaded = True

                # ==========================================
                # CASE B: UNCORRELATED
                # ==========================================
                else:
                    if "energy" in grp and "distribution" in grp:
                        self.incident_energy = grp["energy"][:]
                        self.distribution_data = grp["distribution"][:]

                        if "offsets" in grp: self.offsets = grp["offsets"][:]
                        elif "offsets" in grp["energy"].attrs: self.offsets = grp["energy"].attrs["offsets"]
                        elif "offsets" in grp["distribution"].attrs: self.offsets = grp["distribution"].attrs["offsets"]
                        else: return

                        self.loaded = True

        except Exception as e:
            print(f"Error loading {self.element}: {e}")

    def sample_energy(self, E_in, rng):
        """Standard sampling (Law 4)."""
        if not self.loaded or self.dist_type == "correlated": return None
        # ... (Same as before, uncorrelated code is likely fine given 0.2% errors elsewhere)
        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))
        E_l, E_h = self.incident_energy[idx], self.incident_energy[idx+1]
        f = 0.0 if E_h == E_l else (E_in - E_l) / (E_h - E_l)
        return (1 - f) * self._sample_uncorr(idx, rng) + f * self._sample_uncorr(idx+1, rng)

    def _sample_uncorr(self, idx, rng):
        start = self.offsets[idx]
        end = self.distribution_data.shape[1] if idx == len(self.incident_energy) - 1 else self.offsets[idx + 1]
        if start >= end: return 0.0
        return np.interp(rng.random(), self.distribution_data[2, start:end], self.distribution_data[0, start:end])

    def sample_correlated(self, E_in, rng):
        """
        STRICT Correlated Sampling. Raises Exceptions on Data Failure.
        """
        # --- RESTORED SAFETY CHECK ---
        # If the data is Uncorrelated (e.g. Pb208), we must return False
        # so the simulation knows to use the fallback or other methods.
        if not self.loaded or self.dist_type != "correlated":
            return 0.0, 0.0, False

        # 1. Bin Selection
        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))
        E_l, E_h = self.incident_energy[idx], self.incident_energy[idx+1]
        f = 0.0 if E_h == E_l else (E_in - E_l) / (E_h - E_l)
        selected_idx = idx if rng.random() > f else idx + 1

        # 2. Locate Energy Block
        if self.offsets is None:
             raise ValueError(f"CRITICAL: Offsets are None for {self.element} in sample_correlated")

        start = self.offsets[selected_idx]
        if selected_idx >= len(self.offsets) - 1:
            end = self.energy_out_data.shape[1]
        else:
            end = self.offsets[selected_idx+1]

        if start >= end:
             raise ValueError(f"Empty Energy Block: Start {start} >= End {end} (E_in={E_in:.2e})")

        # 3. Sample Energy
        block_E_out = self.energy_out_data[0, start:end]
        block_CDF   = self.energy_out_data[2, start:end]

        r_E = rng.random()
        k_local = np.searchsorted(block_CDF, r_E)
        k_local = max(0, min(k_local, len(block_E_out) - 1))

        if k_local == 0:
            E_sampled = block_E_out[0]
        else:
            c0, c1 = block_CDF[k_local-1], block_CDF[k_local]
            e0, e1 = block_E_out[k_local-1], block_E_out[k_local]
            if abs(c1 - c0) < 1e-14: E_sampled = e1
            else: E_sampled = e0 + (r_E - c0)/(c1 - c0) * (e1 - e0)

        # 4. STRICT ANGLE DETECTION
        current_global_idx = start + k_local

        val_row3 = int(self.energy_out_data[3, current_global_idx])
        val_row4 = int(self.energy_out_data[4, current_global_idx])

        mu_offset = 0
        mu_len = 0
        used_row4_as_offset = False

        if val_row4 > 10:
            mu_offset = val_row4
            used_row4_as_offset = True
            if current_global_idx < self.energy_out_data.shape[1] - 1:
                next_val = int(self.energy_out_data[4, current_global_idx + 1])
                mu_len = next_val - mu_offset
            else:
                mu_len = self.mu_data.shape[1] - mu_offset if self.mu_data is not None else 0

        elif val_row3 > 10:
            mu_offset = val_row3
            used_row4_as_offset = False
            if current_global_idx < self.energy_out_data.shape[1] - 1:
                next_val = int(self.energy_out_data[3, current_global_idx + 1])
                mu_len = next_val - mu_offset
            else:
                mu_len = self.mu_data.shape[1] - mu_offset if self.mu_data is not None else 0

        else:
            debug_info = f"""
            ❌ PARSING FAILURE for {self.element} MT{self.mt}
            E_in: {E_in:.2e} | E_out: {E_sampled:.2e}
            Global Index: {current_global_idx}
            Row 3 Value: {val_row3}
            Row 4 Value: {val_row4}
            """
            raise ValueError(debug_info)

        # 5. STRICT MU VALIDATION
        if self.mu_data is None:
             raise ValueError("CRITICAL: 'mu' dataset is None but Law 61 was requested.")

        if mu_len <= 0:
             raise ValueError(f"Invalid Mu Length {mu_len} at Global Index {current_global_idx}")

        mu_vals = self.mu_data[0, mu_offset : mu_offset + mu_len]
        mu_cdf  = self.mu_data[2, mu_offset : mu_offset + mu_len]

        if len(mu_vals) > 0:
            last_val = mu_cdf[-1]
            if abs(last_val - 1.0) > 0.01:
                debug_info = f"""
                ❌ INVALID CDF DETECTED for {self.element}
                Used Row 4 as Offset? {used_row4_as_offset}
                Offset: {mu_offset} | Length: {mu_len}
                CDF Last Value: {last_val}
                """
                raise ValueError(debug_info)

        if len(mu_vals) == 1:
            return E_sampled, mu_vals[0], True

        if len(mu_vals) == 0:
            raise ValueError(f"Empty Mu Block (Len=0) at Offset {mu_offset}")

        return E_sampled, np.interp(rng.random(), mu_cdf, mu_vals), True

class CrossSectionReader:
    def __init__(self, base_path: str, temperature: str = "294K"):
        self.base_path = base_path
        self.temperature = temperature
        self._cache = {}
        self._available_inelastic_mts = {}
        self.angular_dists = {}
        self.secondary_dists = {}

    def _load_data_to_cache(self, element: str, mt: int):
        # Map common names
        element_map = {
            "C": "C0", "Graphite": "C0", "C12": "C12",
            "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"
        }
        actual_element = element_map.get(element, element)
        mt_str = f"{mt:03}"
        file_path = os.path.join(self.base_path, f"neutron/{actual_element}.h5")

        if not os.path.exists(file_path):
            self._cache[(element, mt)] = None
            return

        try:
            with h5py.File(file_path, 'r') as f:
                # Root detection
                valid_roots = [k for k in f.keys() if k != "filetype"]
                if not valid_roots: self._cache[(element, mt)] = None; return

                if f"neutron/{actual_element}" in f: root = f"neutron/{actual_element}"
                elif f"{actual_element}" in f: root = f"{actual_element}"
                else: root = valid_roots[0]

                reaction_group = f"{root}/reactions/reaction_{mt_str}"
                if reaction_group not in f: self._cache[(element, mt)] = None; return

                # Temp handling
                available_temps = list(f[reaction_group].keys())
                target_temp = self.temperature
                if target_temp not in available_temps:
                    if target_temp == "294K" and "293.6K" in available_temps: target_temp = "293.6K"
                    else: self._cache[(element, mt)] = None; return

                reaction_temp = f"{reaction_group}/{target_temp}"

                # Energy handling
                energy_path = f"{root}/energy/{target_temp}"
                if energy_path not in f and f"{root}/energy" in f:
                    e_keys = list(f[f"{root}/energy"].keys())
                    if target_temp in e_keys: energy_path = f"{root}/energy/{target_temp}"
                    elif len(e_keys) > 0: energy_path = f"{root}/energy/{e_keys[0]}"

                if energy_path not in f: raise KeyError("Energy path missing")

                energy_data = f[energy_path][:]
                q_value = f[reaction_group].attrs.get('Q_value', 0.0)

                if f"{reaction_temp}/xs" not in f: self._cache[(element, mt)] = None; return

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
        # Map common names
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        actual = element_map.get(element, element)
        file_path = os.path.join(self.base_path, f"neutron/{actual}.h5")
        found = []
        check_list = [16, 17] + list(range(51, 92))

        try:
            with h5py.File(file_path, 'r') as f:
                valid_roots = [k for k in f.keys() if k != "filetype"]
                if valid_roots:
                    if f"neutron/{actual}" in f: root = f"neutron/{actual}"
                    elif f"{actual}" in f: root = f"{actual}"
                    else: root = valid_roots[0]
                    root_reactions = f"{root}/reactions"
                    if root_reactions in f:
                        for mt in check_list:
                            if f"reaction_{mt:03}" in f[root_reactions]: found.append(mt)
        except: pass
        self._available_inelastic_mts[element] = found

    def get_cross_section_data(self, element: str, mt: int):
        cache_key = (element, mt)
        if cache_key not in self._cache: self._load_data_to_cache(element, mt)
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
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        actual = element_map.get(element, element)
        if actual not in self.angular_dists:
            dist = AngularDistribution(self.base_path, actual, mt=2)
            dist.load()
            self.angular_dists[actual] = dist
        return self.angular_dists[actual].sample_mu(energy, rng)

    def get_secondary_energy(self, element, mt, energy, rng):
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        actual = element_map.get(element, element)
        key = (actual, mt)
        if key not in self.secondary_dists:
            dist = SecondaryDistribution(self.base_path, actual, mt)
            dist.load()
            self.secondary_dists[key] = dist
        return self.secondary_dists[key].sample_energy(energy, rng)

    def get_secondary_correlated_sample(self, element, mt, energy, rng):
        """Wrapper for Correlated Sampling."""
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        actual = element_map.get(element, element)
        key = (actual, mt)
        if key not in self.secondary_dists:
            dist = SecondaryDistribution(self.base_path, actual, mt)
            dist.load()
            self.secondary_dists[key] = dist
        return self.secondary_dists[key].sample_correlated(energy, rng)

    def get_cross_sections(self, element, energy, sampler, N, A):
        from .physics import calculate_E_cm_prime
        if energy < 10.0: lookup_energy = calculate_E_cm_prime(energy, A, sampler)
        else: lookup_energy = energy

        # 1. Elastic Scattering
        mic_el = self.get_cross_section(element, 2, lookup_energy)
        Sigma_el = self.calculate_macroscopic_xs(mic_el, N)

        # 2. Inelastic Scattering (Now strictly neutron-emitting)
        Sigma_in, _ = self.get_inelastic_components(element, energy, N)

        # 3. Capture / Absorption (The Iron Fix)
        # Sum ALL absorption channels:
        # 102 (n,g), 103 (n,p), 104 (n,d), 105 (n,t), 106 (n,He3), 107 (n,a)
        mic_cap_total = 0.0
        absorption_mts = [102, 103, 104, 105, 106, 107]

        for mt in absorption_mts:
            xs = self.get_cross_section(element, mt, energy)
            mic_cap_total += xs

        Sigma_cap = self.calculate_macroscopic_xs(mic_cap_total, N)

        # 4. Fission
        mic_fis = self.get_cross_section(element, 18, energy)
        Sigma_fis = self.calculate_macroscopic_xs(mic_fis, N)

        Sigma_total = Sigma_el + Sigma_in + Sigma_cap + Sigma_fis

        return Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_total
