import pytest
import numpy as np
from src.physics import (
    calculate_E_cm_prime,
    calculate_E_prime,
    elastic_scattering,
    sample_new_direction_cosines
)
from src.vt_calc import VelocitySampler
from src.random_number_generator import RNGHandler


class TestEnergyConversions:
    """Test energy frame transformations"""

    def test_E_cm_prime_high_energy(self):
        """At high energy, thermal motion negligible"""
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)

        E_initial = 1e7  # 10 MeV
        E_cm = calculate_E_cm_prime(E_initial, A=2.5, sampler=sampler)

        # Physics check:
        # v_cm = v_lab / (A + 1)
        # v_n_cm = v_lab - v_cm = v_lab * (A / (A + 1))
        # Energy is proportional to v^2
        # Therefore: E_cm = E_lab * (A / (A + 1))^2
        ratio = 2.5 / 3.5
        expected = E_initial * (ratio ** 2)

        assert E_cm == pytest.approx(expected, rel=0.01)

    def test_E_cm_prime_thermal_energy(self):
        """At thermal energy, target motion matters"""
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)

        E_initial = 0.0253
        E_cm = calculate_E_cm_prime(E_initial, A=2.5, sampler=sampler)

        assert E_cm > 0
        # E_cm should be small, order of thermal energy
        assert 0.001 < E_cm < 1

    def test_E_prime_calculation(self):
        """Test scattered energy calculation"""
        rng = RNGHandler(seed=42)
        E_cm_prime = 1e5
        E_initial = 1e6
        A = 2.5

        E_prime, mu_cm = calculate_E_prime(E_cm_prime, E_initial, A, rng)

        assert E_prime > 0
        assert -1 <= mu_cm <= 1


class TestScatteringPhysics:
    """Test elastic scattering mechanics"""

    def test_elastic_scattering_energy_loss(self):
        """Elastic scattering should reduce energy on average"""
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=42)

        E_initial = 1e6
        A = 2.5

        energies = []
        for _ in range(100):
            E_prime, mu_cm, mu_lab = elastic_scattering(E_initial, A, sampler, rng)
            energies.append(E_prime)

        assert all(E < E_initial for E in energies)

        # Average energy ratio for A=2.5 is approx 0.5-0.6
        avg_energy = np.mean(energies)
        assert 0.4 * E_initial < avg_energy < 0.8 * E_initial

    def test_mu_lab_range(self):
        """Lab frame scattering angle cosine must be in [-1, 1]"""
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=42)

        E_initial = 1e6
        A = 2.5

        for _ in range(100):
            E_prime, mu_cm, mu_lab = elastic_scattering(E_initial, A, sampler, rng)
            assert -1.0001 <= mu_lab <= 1.0001

    def test_forward_scattering_bias_heavy_nucleus(self):
        """Heavy nuclei should show forward scattering preference"""
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=42)

        E_initial = 1e6
        # Use A=10 so the forward bias is statistically obvious in 1000 samples
        A = 10.0

        mu_labs = []
        for _ in range(1000):
            _, _, mu_lab = elastic_scattering(E_initial, A, sampler, rng)
            mu_labs.append(mu_lab)

        # Average cosine should be positive (forward bias)
        assert np.mean(mu_labs) > 0.05


class TestThermalUpscatter:
    """Free-gas elastic scattering must allow upscatter at thermal energies so a
    neutron population can equilibrate with the medium (regression guard for the
    static-target bug that left thermal spectra too cold)."""

    def test_upscatter_occurs_at_thermal(self):
        mass = 12 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=2024)
        E0 = 0.0253  # eV, thermal
        A = 12.0
        outs = [elastic_scattering(E0, A, sampler, rng)[0] for _ in range(2000)]
        # A meaningful fraction must GAIN energy — impossible with a static target.
        frac_up = sum(1 for e in outs if e > E0) / len(outs)
        assert frac_up > 0.1

    def test_no_upscatter_when_fast(self):
        # Well above the thermal regime the target is effectively at rest, so the
        # neutron can only lose energy.
        mass = 12 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=7)
        E0 = 1.0e6
        outs = [elastic_scattering(E0, 12.0, sampler, rng)[0] for _ in range(500)]
        assert all(e <= E0 * 1.0001 for e in outs)


class TestDirectionSampling:
    """Test direction cosine sampling"""

    def test_direction_normalization(self):
        """New direction vector should be unit length"""
        rng = RNGHandler(seed=42)
        u, v, w = 1, 0, 0
        mu_lab = 0.5
        u_new, v_new, w_new, phi = sample_new_direction_cosines(u, v, w, mu_lab, rng)

        magnitude = np.sqrt(u_new**2 + v_new**2 + w_new**2)
        assert magnitude == pytest.approx(1.0, abs=1e-6)

    def test_phi_range(self):
        """Azimuthal angle should be in [0, 2π]"""
        rng = RNGHandler(seed=42)
        u, v, w = 1, 0, 0
        mu_lab = 0.5
        for _ in range(100):
            _, _, _, phi = sample_new_direction_cosines(u, v, w, mu_lab, rng)
            assert 0 <= phi <= 2 * np.pi

    def test_mu_lab_projection(self):
        """Verify scattering angle preserved correctly"""
        rng = RNGHandler(seed=42)

        # Test Case 1: Start with z-direction (Singularity case)
        u, v, w = 0, 0, 1
        mu_lab = 0.7
        u_new, v_new, w_new, _ = sample_new_direction_cosines(u, v, w, mu_lab, rng)
        # If incident is (0,0,1), w_new IS the scattering cosine
        assert w_new == pytest.approx(mu_lab, abs=1e-6)

        # Test Case 2: Start with -z direction
        u, v, w = 0, 0, -1
        mu_lab = 0.7
        u_new, v_new, w_new, _ = sample_new_direction_cosines(u, v, w, mu_lab, rng)
        # If incident is (0,0,-1), w_new is -mu_lab
        assert w_new == pytest.approx(-mu_lab, abs=1e-6)


class TestVelocitySampler:
    """Test target velocity sampling"""

    def test_velocity_sampler_initialization(self):
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        assert sampler.beta > 0
        assert sampler.mass == mass

    def test_velocity_sampling_returns_positive(self):
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        vn = 2000
        for _ in range(50):
            v_t = sampler.sample_velocity(vn)
            assert v_t >= 0

    def test_velocity_distribution_thermal_range(self):
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        vn = 2000
        velocities = [sampler.sample_velocity(vn) for _ in range(100)]
        assert np.mean(velocities) < 1000


class TestConservationLaws:
    """Test physical conservation laws"""

    def test_energy_conservation_on_average(self):
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        rng = RNGHandler(seed=42)
        E_initial = 1e6
        A = 2.5
        energy_ratios = []
        for _ in range(1000):
            E_prime, _, _ = elastic_scattering(E_initial, A, sampler, rng)
            energy_ratios.append(E_prime / E_initial)

        avg_ratio = np.mean(energy_ratios)
        assert 0.4 < avg_ratio < 0.7


class TestRNGHandler:
    """Test random number generator"""
    def test_rng_reproducibility(self):
        rng1 = RNGHandler(seed=12345)
        rng2 = RNGHandler(seed=12345)
        seq1 = [rng1.random() for _ in range(10)]
        seq2 = [rng2.random() for _ in range(10)]
        assert seq1 == seq2

    def test_rng_uniform_range(self):
        rng = RNGHandler(seed=42)
        samples = [rng.uniform(5, 10) for _ in range(1000)]
        assert all(5 <= s <= 10 for s in samples)
        assert 7 < np.mean(samples) < 8

    def test_rng_log_uniform(self):
        rng = RNGHandler(seed=42)
        scale = 0.5
        samples = [rng.log_uniform(scale) for _ in range(100)]
        assert all(s > 0 for s in samples)

# Run all tests if executed directly
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
