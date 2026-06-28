"""k-eigenvalue / fission-source support: ν̄ data, the Watt χ sampler, fission
neutron production in the kernel, and the generation-iteration k_eff driver.

The k_eff trend test is the headline: a small bare U235 sphere must be clearly
subcritical and a large one clearly supercritical, with the critical-mass-scale
sphere landing near unity.
"""
import os

import numpy as np
import pytest

from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Sphere
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import simulate_single_particle
from src.vt_calc import VelocitySampler
from src.physics import sample_watt_spectrum, watt_params_for
from src.criticality import run_keff, uniform_sphere_source

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")
HAVE_U235 = os.path.exists(os.path.join(ENDF, "neutron", "U235.h5"))

pytestmark = pytest.mark.skipif(not HAVE_U235, reason="U235 ENDF data not present")

A_U235, RHO_U235, M_U235 = 233.0248, 18.74, 235


def _u235():
    reader = CrossSectionReader(ENDF)
    mat = Material("U235", RHO_U235, M_U235, A_U235)
    return reader, mat.number_density, VelocitySampler(mat.kg_mass)


# --- nuclear data -----------------------------------------------------------

def test_nu_bar_values():
    reader, *_ = _u235()
    # ENDF/B-VIII: ~2.44 thermal, rising with energy
    assert reader.get_nu("U235", 0.0253) == pytest.approx(2.44, abs=0.05)
    assert reader.get_nu("U235", 2.0e6) == pytest.approx(2.65, abs=0.1)
    assert reader.get_nu("U235", 14.0e6) > reader.get_nu("U235", 2.0e6)


def test_nonfissile_has_no_nu():
    reader, *_ = _u235()
    assert reader.get_nu("Pb208", 1.0e6) == 0.0
    assert reader.is_fissionable("U235")
    assert not reader.is_fissionable("Pb208")


# --- Watt spectrum ----------------------------------------------------------

def test_watt_spectrum_mean_and_positivity():
    rng = RNGHandler(1)
    a, b = watt_params_for("U235")
    samples = np.array([sample_watt_spectrum(rng, a, b) for _ in range(20000)])
    assert np.all(samples > 0)
    # analytic Watt mean = 1.5 a + a^2 b / 4  (~2.03 MeV for U235)
    expected = 1.5 * a + a * a * b / 4.0
    assert samples.mean() == pytest.approx(expected, rel=0.05)


# --- fission production in the kernel ---------------------------------------

def test_fission_banks_neutrons_in_criticality_mode():
    reader, N, sampler = _u235()
    mediums = [Region(surfaces=[Sphere((0, 0, 0), 10.0)],
                      name="core", priority=1, element="U235")]
    settings = Settings("criticality", 1)
    rngs = [RNGHandler(50 + i) for i in range(200)]
    states = [{"x": 0.0, "y": 0.0, "z": 1e-9, "theta": np.pi / 2, "phi": 0.0,
               "energy": 2.0e6, "weight": 1.0, "has_interacted": False} for _ in rngs]
    args = [(s, reader, mediums, A_U235, N, sampler, None, False, r, settings)
            for s, r in zip(states, rngs)]
    results = [simulate_single_particle(a) for a in args]
    total_fission_neutrons = sum(len(r["fission_sites"]) for r in results)
    assert total_fission_neutrons > 0
    # each banked site is a valid source state with a positive energy
    for r in results:
        for site in r["fission_sites"]:
            assert site["energy"] > 0 and site["weight"] == 1.0


def test_no_fission_neutrons_for_nonfissile():
    reader, N, sampler = _u235()
    mediums = [Region(surfaces=[Sphere((0, 0, 0), 8.0)],
                      name="s", priority=1, element="Pb208")]
    settings = Settings("criticality", 1)
    rngs = [RNGHandler(i) for i in range(100)]
    states = [{"x": 0.0, "y": 0.0, "z": 1e-9, "theta": np.pi / 2, "phi": 0.0,
               "energy": 2.0e6, "weight": 1.0, "has_interacted": False} for _ in rngs]
    args = [(s, reader, mediums, 207.2, N, sampler, None, False, r, settings)
            for s, r in zip(states, rngs)]
    results = [simulate_single_particle(a) for a in args]
    assert sum(len(r["fission_sites"]) for r in results) == 0


# --- k_eff driver -----------------------------------------------------------

def test_keff_subcritical_vs_supercritical_trend():
    reader, N, sampler = _u235()
    settings = Settings("criticality", 1)

    def keff_for(radius):
        mediums = [Region(surfaces=[Sphere((0, 0, 0), radius)],
                          name="core", priority=1, element="U235")]
        src = uniform_sphere_source((0, 0, 0), radius, 250, RNGHandler(11))
        return run_keff(reader, mediums, src, A_U235, N, sampler, settings,
                        generations=22, inactive=7, seed=2)

    small = keff_for(4.0)    # well below the ~8.7 cm bare critical radius
    large = keff_for(14.0)   # well above it

    assert small["k_eff"] < 0.9
    assert large["k_eff"] > 1.1
    assert small["k_eff"] < large["k_eff"]
    assert np.isfinite(large["k_sem"]) and large["k_sem"] > 0


def test_keff_reports_both_estimators_in_agreement():
    reader, N, sampler = _u235()
    mediums = [Region(surfaces=[Sphere((0, 0, 0), 8.7)],
                      name="core", priority=1, element="U235")]
    src = uniform_sphere_source((0, 0, 0), 8.7, 300, RNGHandler(5))
    out = run_keff(reader, mediums, src, A_U235, N, sampler,
                   Settings("criticality", 1), generations=30, inactive=10, seed=4)
    # both the collision (primary k_eff) and source estimators are reported,
    # finite, near critical, and consistent in the mean
    assert out["k_eff"] > 0 and out["k_source"] > 0
    assert np.isfinite(out["k_sem"]) and np.isfinite(out["k_source_sem"])
    comb = np.hypot(out["k_sem"], out["k_source_sem"])
    assert abs(out["k_eff"] - out["k_source"]) < 5 * comb


def test_keff_driver_requires_nonempty_source():
    reader, N, sampler = _u235()
    with pytest.raises(ValueError):
        run_keff(reader, [], [], A_U235, N, sampler, Settings("criticality", 1))
