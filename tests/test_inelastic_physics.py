import pytest
import numpy as np
from src.physics import inelastic_scattering
from src.vt_calc import VelocitySampler
from src.random_number_generator import RNGHandler

class TestInelasticPhysics:

    def test_inelastic_energy_loss_high_mass(self):
        """
        Test hitting a heavy nucleus (like Lead, A=208).
        Recoil should be minimal, so Final Energy approx Initial - Q.
        """
        # Setup
        mass = 208 * 1.67e-27
        sampler = VelocitySampler(mass=mass)
        rng = RNGHandler(seed=42)

        E_initial = 10.0e6  # 10 MeV
        Q_value = 2.0e6     # 2 MeV cost
        A = 208

        # Run Physics
        E_prime, mu_cm, mu_lab = inelastic_scattering(E_initial, A, Q_value, sampler, rng)

        # Logic:
        # Neutron gives 2 MeV to the nucleus.
        # Heavy nucleus barely moves (recoil is tiny).
        # Neutron should have approx 8 MeV left.
        assert E_prime == pytest.approx(8.0e6, rel=0.05) # Allow 5% variance for recoil/thermal

    def test_inelastic_threshold_check(self):
        """
        Test what happens if neutron has LESS energy than Q.
        Should behave like Elastic scattering (no energy loss to Q).
        """
        mass = 208 * 1.67e-27
        sampler = VelocitySampler(mass=mass)
        rng = RNGHandler(seed=42)

        E_initial = 1.0e6   # 1 MeV
        Q_value = 2.0e6     # 2 MeV cost (Too expensive!)
        A = 208

        E_prime, _, _ = inelastic_scattering(E_initial, A, Q_value, sampler, rng)

        # Should behave like elastic scattering off a heavy target -> Energy nearly conserved
        assert E_prime == pytest.approx(E_initial, rel=0.01)

    def test_light_nucleus_recoil(self):
        """
        Test hitting a light nucleus (like Boron, A=10).
        Kinetic energy loss to recoil will be significant + Q value loss.
        """
        mass = 10 * 1.67e-27
        sampler = VelocitySampler(mass=mass)
        rng = RNGHandler(seed=42)

        E_initial = 10.0e6 # 10 MeV
        Q_value = 2.0e6    # 2 MeV
        A = 10

        E_prime, _, _ = inelastic_scattering(E_initial, A, Q_value, sampler, rng)

        # Total Energy Loss = Q + Recoil
        # Recoil for A=10 is significant.
        # Expected is definitely less than 8.0 MeV
        assert E_prime < 7.8e6
        assert E_prime > 5.0e6 # Sanity check
