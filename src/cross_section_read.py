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
    """Secondary energy-angle reader for (n,xn)/continuum reactions.

    Handles ENDF Law 61 (correlated), Law 4 (uncorrelated) and Kalbach-Mann.
    Any layout it cannot parse leaves the sampler reporting success=False so the
    caller falls back rather than aborting the worker.
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

        self.ref_frame = "lab"

    def load(self):
        if self.loaded: return
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        fname = element_map.get(self.element, self.element)
        file_path = os.path.join(self.base_path, f"neutron/{fname}.h5")

        if not os.path.exists(file_path): return

        try:
            with h5py.File(file_path, 'r') as f:
                root = f"neutron/{fname}" if f"neutron/{fname}" in f else fname
                reaction = f"{root}/reactions/reaction_{self.mt:03}"
                base = f"{reaction}/product_0/distribution_0"
                if base not in f: return
                grp = f[base]

                if "type" in grp.attrs:
                    self.dist_type = grp.attrs["type"].decode('utf-8')

                # frame from the reaction group (center_of_mass: 1 = CM, 0 = lab)
                if reaction in f and int(f[reaction].attrs.get("center_of_mass", 0)) == 1:
                    self.ref_frame = "cm"

                # correlated (Law 61): per-node energy_out blocks plus mu tables
                if self.dist_type == "correlated":
                    self.incident_energy = grp["energy"][:]

                    if "offsets" in grp: self.offsets = grp["offsets"][:]
                    elif "offsets" in grp["energy"].attrs: self.offsets = grp["energy"].attrs["offsets"]
                    elif "energy_out" in grp and "offsets" in grp["energy_out"].attrs:
                         self.offsets = grp["energy_out"].attrs["offsets"]
                    else:
                        return

                    self.energy_out_data = grp["energy_out"][:]
                    if "mu" in grp: self.mu_data = grp["mu"][:]

                    self.loaded = True

                # uncorrelated (Law 4) and Kalbach-Mann: a single distribution block
                else:
                    if "energy" in grp and "distribution" in grp:
                        self.incident_energy = grp["energy"][:]
                        self.distribution_data = grp["distribution"][:]

                        if "offsets" in grp: self.offsets = grp["offsets"][:]
                        elif "offsets" in grp["energy"].attrs: self.offsets = grp["energy"].attrs["offsets"]
                        elif "offsets" in grp["distribution"].attrs: self.offsets = grp["distribution"].attrs["offsets"]
                        else: return

                        self.loaded = True

        except Exception:
            return

    def sample_energy(self, E_in, rng):
        """Law 4 (uncorrelated) outgoing energy, interpolated between incident nodes."""
        if not self.loaded or self.dist_type == "correlated": return None
        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))
        E_l, E_h = self.incident_energy[idx], self.incident_energy[idx+1]
        f = 0.0 if E_h == E_l else (E_in - E_l) / (E_h - E_l)
        f = max(0.0, min(1.0, f))   # clamp to the grid
        return (1 - f) * self._sample_uncorr(idx, rng) + f * self._sample_uncorr(idx+1, rng)

    def _sample_uncorr(self, idx, rng):
        start = self.offsets[idx]
        end = self.distribution_data.shape[1] if idx == len(self.incident_energy) - 1 else self.offsets[idx + 1]
        if start >= end: return 0.0
        return np.interp(rng.random(), self.distribution_data[2, start:end], self.distribution_data[0, start:end])

    def _block_span(self, data, i):
        """[start, end) column slice of the i-th incident-energy block."""
        start = self.offsets[i]
        end = data.shape[1] if i >= len(self.offsets) - 1 else self.offsets[i + 1]
        return start, end

    def _unit_base(self, E_sampled, block_E_out, idx, f, data):
        """Unit-base interpolation: rescale E_sampled onto the [E_min, E_max]
        interpolated between the two bracketing incident nodes, so the tail
        endpoint tracks the incident energy."""
        e_min_sel, e_max_sel = block_E_out[0], block_E_out[-1]
        span_sel = e_max_sel - e_min_sel
        frac = 0.0 if span_sel <= 0 else (E_sampled - e_min_sel) / span_sel

        s_lo, e_lo = self._block_span(data, idx)
        s_hi, e_hi = self._block_span(data, idx + 1)
        lo_min, lo_max = data[0, s_lo], data[0, e_lo - 1]
        hi_min, hi_max = data[0, s_hi], data[0, e_hi - 1]
        e_min = lo_min + f * (hi_min - lo_min)
        e_max = lo_max + f * (hi_max - lo_max)
        return e_min + frac * (e_max - e_min)

    def sample_correlated(self, E_in, rng):
        """Law 61 (E_out, mu, True), or (0, 0, False) when the layout is uninterpretable."""
        if not self.loaded:
            return 0.0, 0.0, False
        if self.dist_type == "kalbach-mann":
            return self._sample_kalbach(E_in, rng)
        if self.dist_type != "correlated":
            return 0.0, 0.0, False
        if self.offsets is None:
            return 0.0, 0.0, False

        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))
        E_l, E_h = self.incident_energy[idx], self.incident_energy[idx+1]
        f = 0.0 if E_h == E_l else (E_in - E_l) / (E_h - E_l)
        f = max(0.0, min(1.0, f))   # clamp to the grid
        selected_idx = idx if rng.random() > f else idx + 1

        # sample the outgoing-energy shape from the chosen incident node
        start, end = self._block_span(self.energy_out_data, selected_idx)
        if start >= end:
            return 0.0, 0.0, False

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

        E_sampled = self._unit_base(E_sampled, block_E_out, idx, f, self.energy_out_data)

        # mu offset lives in row 4 (or row 3 in older layouts); <= 10 means none
        current_global_idx = start + k_local
        val_row3 = int(self.energy_out_data[3, current_global_idx])
        val_row4 = int(self.energy_out_data[4, current_global_idx])

        mu_offset = 0
        mu_len = 0
        if val_row4 > 10:
            mu_offset = val_row4
            if current_global_idx < self.energy_out_data.shape[1] - 1:
                mu_len = int(self.energy_out_data[4, current_global_idx + 1]) - mu_offset
            else:
                mu_len = self.mu_data.shape[1] - mu_offset if self.mu_data is not None else 0
        elif val_row3 > 10:
            mu_offset = val_row3
            if current_global_idx < self.energy_out_data.shape[1] - 1:
                mu_len = int(self.energy_out_data[3, current_global_idx + 1]) - mu_offset
            else:
                mu_len = self.mu_data.shape[1] - mu_offset if self.mu_data is not None else 0
        else:
            return 0.0, 0.0, False

        if self.mu_data is None or mu_len <= 0:
            return 0.0, 0.0, False

        mu_vals = self.mu_data[0, mu_offset : mu_offset + mu_len]
        mu_cdf  = self.mu_data[2, mu_offset : mu_offset + mu_len]

        if len(mu_vals) == 0:
            return 0.0, 0.0, False
        if abs(mu_cdf[-1] - 1.0) > 0.01:   # malformed CDF
            return 0.0, 0.0, False
        if len(mu_vals) == 1:
            return E_sampled, mu_vals[0], True

        return E_sampled, np.interp(rng.random(), mu_cdf, mu_vals), True

    def _sample_kalbach(self, E_in, rng):
        """Kalbach-Mann (E_out, mu). Block rows are [E_out, pdf, cdf, r, a];
        energy from the cdf, cosine from f(mu) ~ cosh(a mu) + r sinh(a mu)."""
        dist = self.distribution_data
        if dist is None or self.offsets is None or dist.shape[0] < 5:
            return 0.0, 0.0, False

        idx = np.searchsorted(self.incident_energy, E_in) - 1
        idx = max(0, min(idx, len(self.incident_energy) - 2))
        E_l, E_h = self.incident_energy[idx], self.incident_energy[idx + 1]
        f = 0.0 if E_h == E_l else (E_in - E_l) / (E_h - E_l)
        f = max(0.0, min(1.0, f))   # clamp to the grid
        sel = idx if rng.random() > f else idx + 1

        start, end = self._block_span(dist, sel)
        if start >= end:
            return 0.0, 0.0, False

        E_out = dist[0, start:end]
        cdf   = dist[2, start:end]
        r_arr = dist[3, start:end]
        a_arr = dist[4, start:end]

        xi = rng.random()
        k = np.searchsorted(cdf, xi)
        k = max(0, min(k, len(E_out) - 1))
        if k == 0:
            E_s, r, a = E_out[0], r_arr[0], a_arr[0]
        else:
            c0, c1 = cdf[k - 1], cdf[k]
            e0, e1 = E_out[k - 1], E_out[k]
            E_s = e1 if abs(c1 - c0) < 1e-14 else e0 + (xi - c0) / (c1 - c0) * (e1 - e0)
            r, a = r_arr[k], a_arr[k]

        # Unit-base interpolation of the outgoing energy between incident grids.
        E_s = self._unit_base(E_s, E_out, idx, f, dist)
        return E_s, self._kalbach_mu(float(a), float(r), rng), True

    @staticmethod
    def _kalbach_mu(a, r, rng):
        """Sample the emission cosine from the Kalbach angular law."""
        a = min(max(a, 1e-6), 700.0)   # guard against overflow / a=0
        if rng.random() <= r:
            xi = rng.random()
            return np.log(xi * np.exp(a) + (1.0 - xi) * np.exp(-a)) / a
        T = (2.0 * rng.random() - 1.0) * np.sinh(a)
        return np.log(T + np.sqrt(T * T + 1.0)) / a


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

    def get_secondary_frame(self, element, mt):
        """Reference frame of an (n,xn) secondary distribution: 'cm' or 'lab'."""
        element_map = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27", "Fe": "Fe56", "Pb": "Pb208"}
        actual = element_map.get(element, element)
        key = (actual, mt)
        if key not in self.secondary_dists:
            dist = SecondaryDistribution(self.base_path, actual, mt)
            dist.load()
            self.secondary_dists[key] = dist
        return self.secondary_dists[key].ref_frame

    def get_cross_sections(self, element, energy, sampler, N, A):
        from .physics import calculate_E_cm_prime
        if energy < 10.0: lookup_energy = calculate_E_cm_prime(energy, A, sampler)
        else: lookup_energy = energy

        # elastic
        mic_el = self.get_cross_section(element, 2, lookup_energy)
        Sigma_el = self.calculate_macroscopic_xs(mic_el, N)

        # neutron-emitting inelastic
        Sigma_in, _ = self.get_inelastic_components(element, energy, N)

        # absorption: 102 (n,g), 103 (n,p), 104 (n,d), 105 (n,t), 106 (n,He3), 107 (n,a)
        mic_cap_total = 0.0
        absorption_mts = [102, 103, 104, 105, 106, 107]

        for mt in absorption_mts:
            xs = self.get_cross_section(element, mt, energy)
            mic_cap_total += xs

        Sigma_cap = self.calculate_macroscopic_xs(mic_cap_total, N)

        # fission
        mic_fis = self.get_cross_section(element, 18, energy)
        Sigma_fis = self.calculate_macroscopic_xs(mic_fis, N)

        Sigma_total = Sigma_el + Sigma_in + Sigma_cap + Sigma_fis

        return Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_total
