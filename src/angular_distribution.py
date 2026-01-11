import numpy as np
import h5py
import os

class AngularDistribution:
    def __init__(self, base_path, element, mt=2):
        self.base_path = base_path
        self.element = element
        self.mt = mt
        self.energy_grid = None
        self.mu_data = None # Shape (3, N) -> [mu, pdf, cdf]
        self.offsets = None
        self.loaded = False

    def load(self):
        """Loads the tabulated angular data from HDF5."""
        if self.loaded: return

        # Try both path conventions
        path_1 = f"neutron/{self.element}"
        path_2 = f"{self.element}"

        file_path = os.path.join(self.base_path, f"neutron/{self.element}.h5")
        if not os.path.exists(file_path):
            return # Fail silently, physics will fallback to isotropic

        try:
            with h5py.File(file_path, 'r') as f:
                base_prefix = path_1 if path_1 in f else path_2
                # Standard ENDFB80 location for elastic scattering (MT=2)
                # .../reactions/reaction_002/product_0/distribution_0/angle
                angle_path = f"{base_prefix}/reactions/reaction_{self.mt:03}/product_0/distribution_0/angle"

                if angle_path not in f:
                    return

                grp = f[angle_path]

                if "energy" in grp and "mu" in grp:
                    self.energy_grid = grp["energy"][:]
                    self.mu_data = grp["mu"][:] # Load giant table (3, 5919)

                    # Extract offsets from ATTRIBUTE
                    if "offsets" in grp["mu"].attrs:
                        self.offsets = grp["mu"].attrs["offsets"]
                    else:
                        # Fallback: assuming uniform? Unlikely.
                        return

                    self.loaded = True

        except Exception as e:
            print(f"Warning: Failed to load angular data: {e}")

    def sample_mu(self, energy_eV, rng):
        """
        Samples a scattering cosine (mu) for a given incident energy.
        Returns: mu_cm (Center of Mass frame)
        """
        if not self.loaded:
            # Fallback to Isotropic if data missing
            return 2 * rng.random() - 1

        # 1. Find Energy Bins (Binary Search)
        # Find i such that E[i] <= E < E[i+1]
        idx = np.searchsorted(self.energy_grid, energy_eV) - 1

        # Clamp to range
        idx = max(0, min(idx, len(self.energy_grid) - 2))

        # 2. Get Interpolation Factor f
        E_low = self.energy_grid[idx]
        E_high = self.energy_grid[idx+1]
        f = (energy_eV - E_low) / (E_high - E_low)

        # 3. Sample mu at E_low
        mu_low = self._sample_from_block(idx, rng)

        # 4. Sample mu at E_high
        mu_high = self._sample_from_block(idx + 1, rng)

        # 5. Interpolate (Standard linear interpolation of the sampled value)
        return (1 - f) * mu_low + f * mu_high

    def _sample_from_block(self, energy_idx, rng):
        """Samples mu from the specific block for energy index i."""
        start = self.offsets[energy_idx]

        # Determine end index
        if energy_idx == len(self.energy_grid) - 1:
            end = self.mu_data.shape[1]
        else:
            end = self.offsets[energy_idx + 1]

        # Extract the block (View)
        # Row 0: mu, Row 2: cdf
        block_mu = self.mu_data[0, start:end]
        block_cdf = self.mu_data[2, start:end]

        # Random sample on CDF [0, 1]
        r = rng.random()

        # Binary search on CDF to find position
        # np.interp performs linear interpolation on the tabulated CDF
        # x = interp(r, xp=cdf, fp=mu)
        return np.interp(r, block_cdf, block_mu)
