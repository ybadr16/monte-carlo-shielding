import pytest
import numpy as np
from src.material import Material
from src.cross_section_read import CrossSectionReader
from src.vt_calc import VelocitySampler


class TestMaterial:
    """Test Material class calculations"""

    def test_material_initialization(self):
        """Test basic material creation"""
        lead = Material(
            name="Lead",
            density=11.35,
            atomic_mass=208,
            atomic_weight_ratio=2.5
        )

        assert lead.name == "Lead"
        assert lead.density == 11.35
        assert lead.atomic_mass == 208
        assert lead.atomic_weight_ratio == 2.5

    def test_number_density_calculation(self):
        """Test number density calculation"""
        # Lead: density = 11.35 g/cm³, atomic mass = 208 g/mol
        lead = Material(name="Lead", density=11.35, atomic_mass=208)

        # N = (ρ × N_A) / A
        # N = (11.35 × 6.022e23) / 208
        expected_N = (11.35 * 6.022e23) / 208

        assert lead.number_density == pytest.approx(expected_N, rel=1e-6)

    def test_atomic_mass_kg(self):
        """Test atomic mass in kg calculation"""
        lead = Material(name="Lead", density=11.35, atomic_mass=208)

        # m_atom = A / (N_A × 1000)  [factor of 1000 for g->kg]
        expected_mass = (208 / 6.022e23) * 1e-3

        assert lead.kg_mass == pytest.approx(expected_mass, rel=1e-6)

    def test_different_materials(self):
        """Test calculations for different materials"""
        materials = [
            Material("Hydrogen", density=0.00009, atomic_mass=1),
            Material("Carbon", density=2.26, atomic_mass=12),
            Material("Iron", density=7.87, atomic_mass=56),
            Material("Uranium", density=19.1, atomic_mass=238),
        ]

        for mat in materials:
            # Number density should be positive
            assert mat.number_density > 0

            # Higher density should give higher number density (for similar A)
            assert mat.kg_mass > 0

    def test_material_repr(self):
        """Test string representation"""
        lead = Material(name="Lead", density=11.35, atomic_mass=208)
        repr_str = repr(lead)

        assert "Lead" in repr_str
        assert "11.35" in repr_str
        assert "208" in repr_str


class TestCrossSectionReader:
    """Test CrossSectionReader functionality"""

    @pytest.fixture
    def reader(self):
        """Fixture to create CrossSectionReader"""
        return CrossSectionReader(base_path="./endfb")

    def test_reader_initialization(self, reader):
        """Test reader initializes correctly"""
        assert reader.base_path == "./endfb"

    def test_get_cross_section_elastic(self, reader):
        """Test retrieving elastic scattering XS (MT=2)"""
        # Test for Pb-208 elastic scattering at 1 MeV
        xs = reader.get_cross_section(
            element="Pb208",
            mt=2,  # Elastic scattering
            energy=1e6  # 1 MeV
        )

        # Cross-section should be positive and reasonable (1-10 barns)
        assert xs > 0
        assert 0.1 < xs < 20  # Typical range for lead

    def test_get_cross_section_capture(self, reader):
        """Test retrieving radiative capture XS (MT=102)"""
        xs = reader.get_cross_section(
            element="Pb208",
            mt=102,  # Radiative capture
            energy=1e6
        )

        assert xs >= 0  # Can be zero at high energies

    def test_cross_section_below_threshold(self, reader):
        """Test XS below reaction threshold returns zero"""
        xs = reader.get_cross_section(
            element="Pb208",
            mt=2,
            energy=1e-5  # Very low energy, below some thresholds
        )

        # Should return zero or very small value
        assert xs >= 0

    def test_macroscopic_xs_calculation(self, reader):
        """Test macroscopic cross-section calculation"""
        microscopic_xs = 5.0  # barns
        N = 3.3e22  # atoms/cm³ (typical for lead)

        Sigma = reader.calculate_macroscopic_xs(microscopic_xs, N)

        # Σ = σ × N, where σ is in cm² and N in atoms/cm³
        # σ = 5 × 10⁻²⁴ cm², so Σ = 5 × 10⁻²⁴ × 3.3e22 = 0.165 cm⁻¹
        expected = 5e-24 * 3.3e22

        assert Sigma == pytest.approx(expected, rel=1e-6)

    def test_macroscopic_xs_zero_density(self, reader):
        """Test macroscopic XS with zero density"""
        Sigma = reader.calculate_macroscopic_xs(5.0, 0.0)
        assert Sigma == 0.0

    def test_macroscopic_xs_negative_inputs(self, reader):
        """A non-physical (negative) microscopic XS is treated as no
        interaction and returns 0.0 — the reader is graceful by design so the
        transport kernel never has to guard every call."""
        assert reader.calculate_macroscopic_xs(-5.0, 1e22) == 0.0
        assert reader.calculate_macroscopic_xs(0.0, 1e22) == 0.0

    def test_get_cross_sections_wrapper(self, reader):
        """Test the convenience wrapper that returns all macroscopic XS.

        Current signature:
            get_cross_sections(element, energy, sampler, N, A)
            -> (Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_total)
        """
        A = 207.2
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass, temperature=294)
        N = 3.3e22

        Sigma_el, Sigma_in, Sigma_cap, Sigma_f, Sigma_t = reader.get_cross_sections(
            "Pb208", 1e6, sampler, N, A
        )

        # All should be non-negative
        for s in (Sigma_el, Sigma_in, Sigma_cap, Sigma_f, Sigma_t):
            assert s >= 0

        # Total should equal the sum of the channels
        assert Sigma_t == pytest.approx(
            Sigma_el + Sigma_in + Sigma_cap + Sigma_f, rel=1e-9)

        # For Pb-208 (non-fissile), fission should be zero
        assert Sigma_f == 0

    def test_invalid_element(self, reader):
        """A missing data file yields 0.0 (graceful), not an exception."""
        xs = reader.get_cross_section(element="Xx999", mt=2, energy=1e6)
        assert xs == 0.0

    def test_invalid_mt_number(self, reader):
        """An MT number with no reaction in the data file returns 0.0 rather
        than raising — absent channels simply contribute nothing."""
        xs = reader.get_cross_section(element="Pb208", mt=9999, energy=1e6)
        assert xs == 0.0


class TestCrossSectionPhysics:
    """Test physical behavior of cross sections"""

    def test_xs_energy_dependence(self):
        """Test that XS changes with energy as expected"""
        reader = CrossSectionReader("./endfb")

        energies = [1e3, 1e4, 1e5, 1e6]  # eV
        xs_values = [
            reader.get_cross_section("Pb208", 2, E)
            for E in energies
        ]

        # All should be positive
        assert all(xs > 0 for xs in xs_values)

        # At low energies, XS typically increases (1/v behavior for some)
        # At high energies, XS typically decreases
        # This is element-dependent, so just check they're different
        assert len(set(xs_values)) > 1  # Not all the same

    def test_absorption_vs_scattering(self):
        """Test that scattering dominates at high energy"""
        reader = CrossSectionReader("./endfb")

        high_energy = 1e6  # 1 MeV

        xs_scatter = reader.get_cross_section("Pb208", 2, high_energy)
        xs_capture = reader.get_cross_section("Pb208", 102, high_energy)

        # At 1 MeV, scattering should dominate
        assert xs_scatter > xs_capture


class TestIntegrationMaterialsAndCrossSections:
    """Integration tests combining materials and cross sections"""

    def test_mean_free_path_calculation(self):
        """Test mean free path calculation: λ = 1/Σ"""
        reader = CrossSectionReader("./endfb")
        lead = Material("Lead", density=11.35, atomic_mass=208,
                        atomic_weight_ratio=207.2)

        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass)

        *_, Sigma_t = reader.get_cross_sections(
            "Pb208", 1e6, sampler, lead.number_density, lead.atomic_weight_ratio
        )

        # Mean free path
        mfp = 1 / Sigma_t

        # For lead at 1 MeV, MFP should be on the order of several cm
        assert 1 < mfp < 20  # Reasonable range

    def test_scattering_probability(self):
        """Test that scattering probability is computed correctly"""
        reader = CrossSectionReader("./endfb")
        lead = Material("Lead", density=11.35, atomic_mass=208,
                        atomic_weight_ratio=207.2)
        mass = 208 * 1.674927471e-27
        sampler = VelocitySampler(mass=mass)

        Sigma_el, Sigma_in, Sigma_cap, Sigma_f, Sigma_t = reader.get_cross_sections(
            "Pb208", 1e6, sampler, lead.number_density, lead.atomic_weight_ratio
        )

        # Scattering probability (elastic + inelastic)
        p_scatter = (Sigma_el + Sigma_in) / Sigma_t

        # Should be between 0 and 1
        assert 0 <= p_scatter <= 1

        # For lead at 1 MeV, scattering dominates
        assert p_scatter > 0.9


# Validation tests (can be marked as slow)
class TestValidationBenchmarks:
    """Validation against known physics benchmarks"""

    @pytest.mark.slow
    def test_1_over_v_law_thermal(self):
        """Test that capture follows 1/v law at thermal energies"""
        reader = CrossSectionReader("./endfb")

        # Test at thermal energies where 1/v applies
        E1 = 0.0253  # eV (thermal)
        E2 = 0.0253 * 4  # 4× thermal

        xs1 = reader.get_cross_section("Pb208", 102, E1)
        xs2 = reader.get_cross_section("Pb208", 102, E2)

        # Should follow: xs1/xs2 ≈ sqrt(E2/E1)
        ratio = xs1 / xs2 if xs2 > 0 else 0
        expected_ratio = np.sqrt(E2 / E1)

        # Allow 20% tolerance (not all isotopes follow this exactly)
        if xs1 > 0 and xs2 > 0:
            assert ratio == pytest.approx(expected_ratio, rel=0.2)


# Run all tests if executed directly
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
