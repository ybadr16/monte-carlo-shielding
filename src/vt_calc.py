import math

import numpy as np

_SQRT_PI = math.sqrt(math.pi)
_HALF_PI = 0.5 * math.pi


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

    # The two target-speed shapes are sampled directly (inverse-CDF / gamma
    # decomposition) instead of by acceptance-rejection. This is the standard
    # SVT method (as in MCNP/OpenMC) and draws the identical distributions:
    #   x^3 exp(-x^2):  x = sqrt(-ln(r1 r2))                  [t=x^2 ~ Gamma(2)]
    #   x^2 exp(-x^2):  x = sqrt(-ln r1 - ln r2 cos^2(pi/2 r3)) [t=x^2 ~ Gamma(3/2)]
    # The old rejection proposal was uniform on [0,10] with ~10% acceptance, so
    # this removes the dominant per-collision cost in thermal problems.

    def sample_velocity(self, vn, rng=None, max_attempts=1000, return_mu=False):
        """Free-gas (SVT) target velocity by acceptance-rejection on the relative speed.

        With return_mu, returns (v_t, mu); the speed and its cosine must be used
        together or the speed-angle correlation (detailed balance) is lost.

        `rng` is the per-history `RNGHandler`; when supplied the target-velocity
        draws come from that reproducible stream rather than the global NumPy RNG,
        so thermal (sub-cutoff) histories are as seed-reproducible as the rest of
        the transport. It falls back to the global RNG only when no stream is given.
        """
        beta = self.beta
        p_first = 2.0 / (_SQRT_PI * vn + 2.0)
        draw = rng.random_array if rng is not None else np.random.random
        for _attempt in range(max_attempts):
            # one batched draw per attempt: r0 branch, r1-r3 speed shape,
            # r4 cosine, r5 acceptance. 1.0 - r keeps log arguments positive.
            r = draw(6)
            if r[0] < p_first:
                x = math.sqrt(-math.log((1.0 - r[1]) * (1.0 - r[2])))
            else:
                c = math.cos(_HALF_PI * r[3])
                x = math.sqrt(-math.log(1.0 - r[1]) - math.log(1.0 - r[2]) * c * c)

            v_t = x / beta
            mu = 2.0 * r[4] - 1.0

            # accept with probability proportional to the relative speed
            acceptance_prob = math.sqrt(vn * vn + v_t * v_t - 2.0 * vn * v_t * mu) / (vn + v_t)
            if r[5] < acceptance_prob:
                return (v_t, mu) if return_mu else v_t

        raise ValueError("Failed to find an accepted sample within the maximum attempts")

