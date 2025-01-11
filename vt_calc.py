import numpy as np


class VelocitySampler:
    """A class for sampling velocities based on nuclear physics calculations."""

    def __init__(self, mass, temperature = 294):
        """
        Initialize the VelocitySampler with mass and temperature parameters.

        Parameters:
            mass (float): Mass of the nucleus (in appropriate units)
            temperature (float): Temperature (in appropriate units)
        """
        self.mass = mass
        self.temperature = temperature
        self.k_B = 1.380649e-23  # Boltzmann constant in J/K
        self.beta = self._calculate_beta()

    def _calculate_beta(self):
        """Calculate beta based on the mass of the nucleus and temperature."""
        return np.sqrt(self.mass / (2 * self.k_B * self.temperature))

    def _sample_from_2x3_exp_neg_x2(self):
        """Samples x from the distribution 2x^3 * exp(-x^2)."""
        while True:
            x = np.random.uniform(0, 10)  # Choose an upper limit suitable for the problem
            p = 2 * x**3 * np.exp(-x**2)
            if np.random.uniform(0, 1) < p:
                return x

    def _sample_from_4pi_x2_exp_neg_x2(self):
        """Samples x from the distribution 4πx^2 * exp(-x^2)."""
        while True:
            x = np.random.uniform(0, 10)  # Choose an upper limit suitable for the problem
            p = 4 * x**2 * np.exp(-x**2)/np.sqrt(np.pi)
            if np.random.uniform(0, 1) < p:
                return x

    def sample_velocity(self, vn, max_attempts=1000):
        """
        Samples v_T based on the given process with input vn.

        Parameters:
            vn (float): Input velocity magnitude (from other parts of the code)
            max_attempts (int): Maximum attempts to find an accepted sample

        Returns:
            float: Accepted sample v_T

        Raises:
            ValueError: If no accepted sample is found within max_attempts
        """
        for attempt in range(max_attempts):
            # Step 1: Sample from q(x)
            xi1 = np.random.uniform(0, 1)
            if xi1 < 2 / (np.sqrt(np.pi) * vn + 2):
                # Sample from 2x^3 * exp(-x^2)
                x = self._sample_from_2x3_exp_neg_x2()
            else:
                # Sample from 4πx^2 * exp(-x^2)
                x = self._sample_from_4pi_x2_exp_neg_x2()

            # Convert x to v_T
            v_t = x / self.beta

            # Sample mu
            xi2 = np.random.uniform(0, 1)
            mu = 2 * xi2 - 1

            # Step 2: Acceptance criterion
            xi3 = np.random.uniform(0, 1)
            acceptance_prob = np.sqrt((vn**2 + v_t**2 - 2 * vn * v_t * mu)) / (vn + v_t)

            if xi3 < acceptance_prob:
                return v_t

        raise ValueError("Failed to find an accepted sample within the maximum attempts")
'''
# Example inputs
vn = 3.0  # Velocity magnitude
# Create a sampler instance with specific mass and temperature
sampler = VelocitySampler(mass=1.0e-27, temperature=300)

# Sample a velocity
try:
    v_t = sampler.sample_velocity(vn=3)
    print(f"Sampled velocity: {v_t}")
except ValueError as e:
    print(f"Sampling failed: {e}")
'''
