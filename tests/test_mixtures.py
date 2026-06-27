"""Material-mixture support: the Nuclide/composition data model and the
multi-isotope transport path in the kernel.

The cornerstone is backward compatibility: a region carrying a single-isotope
`composition` must be byte-identical to the legacy `element=` path (the kernel
draws no isotope-selection random number when there is only one isotope), and a
mixture must reproduce the same physics as the equivalent single isotope at the
same total number density.
"""
import os

import numpy as np
import pytest

from src.cross_section_read import CrossSectionReader
from src.material import Material, Nuclide
from src.medium import Region, Sphere
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import simulate_single_particle
from src.vt_calc import VelocitySampler

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")
HAVE_DATA = os.path.exists(os.path.join(ENDF, "neutron", "Pb208.h5"))

pytestmark = pytest.mark.skipif(not HAVE_DATA, reason="ENDF data not present")


def _run(mediums, *, energy=2.0e6, n=200, mode="shielding",
         A=207.2, N=None, sampler=None):
    """Track `n` particles from an isotropic point source at the origin."""
    reader = CrossSectionReader(ENDF)
    if N is None:
        N = Material("Pb208", 11.35, 208, A).number_density
    if sampler is None:
        sampler = VelocitySampler(Material("Pb208", 11.35, 208, A).kg_mass)
    settings = Settings(mode, n)
    rngs = [RNGHandler(100 + i) for i in range(n)]
    states = [{
        "x": 0.0, "y": 0.0, "z": 1e-9,
        "theta": np.pi / 2, "phi": 0.0,
        "energy": energy, "weight": 1.0, "has_interacted": False,
    } for _ in rngs]
    args = [(s, reader, mediums, A, N, sampler, None, False, r, settings)
            for s, r in zip(states, rngs)]
    return [simulate_single_particle(a) for a in args]


def _leakage(results, n):
    return sum(r["final_weight"] for r in results) / n


# --- data model -------------------------------------------------------------

def test_nuclide_carries_density_and_sampler():
    nuc = Nuclide("U235", number_density=4.8e22, atomic_weight_ratio=233.0,
                  atomic_mass=235.04)
    assert nuc.element == "U235"
    assert nuc.number_density == pytest.approx(4.8e22)
    assert nuc.atomic_weight_ratio == 233.0
    # the sampler mass is built from the molar mass (kg per atom)
    assert nuc.sampler.mass == pytest.approx((235.04 / 6.022e23) * 1e-3)


def test_mixture_single_component_matches_material_density():
    # A one-component mixture must reproduce the plain Material number density.
    mat = Material("Lead", 11.35, 208, 207.2)
    comp = Material.mixture("Lead", 11.35, [("Pb208", 207.2, 1.0, 208)])
    assert len(comp) == 1
    assert comp[0].number_density == pytest.approx(mat.number_density, rel=1e-9)


def test_mixture_atom_fractions_partition_total_density():
    # Number densities split by normalised atom fraction and sum to the total.
    comp = Material.mixture("PbMix", 11.35,
                            [("Pb207", 205.0, 0.25, 207),
                             ("Pb208", 206.0, 0.75, 208)])
    total = sum(c.number_density for c in comp)
    assert comp[0].number_density == pytest.approx(0.25 * total, rel=1e-9)
    assert comp[1].number_density == pytest.approx(0.75 * total, rel=1e-9)


def test_mixture_rejects_nonpositive_fractions():
    with pytest.raises(ValueError):
        Material.mixture("bad", 1.0, [("Pb208", 207.2, 0.0)])


# --- transport equivalence --------------------------------------------------

def test_single_isotope_composition_is_byte_identical_to_legacy():
    # The whole backward-compatibility guarantee: a one-isotope composition must
    # give the exact same per-history results as the legacy element= path,
    # because no isotope-selection random number is drawn for a single isotope.
    A, N = 207.2, Material("Pb208", 11.35, 208, 207.2).number_density
    legacy = [Region(surfaces=[Sphere((0, 0, 0), 8.0)],
                     name="s", priority=1, element="Pb208")]
    mixed = [Region(surfaces=[Sphere((0, 0, 0), 8.0)], name="s", priority=1,
                    composition=[Nuclide("Pb208", N, A, atomic_mass=208)])]

    r_legacy = _run(legacy, n=150)
    r_mixed = _run(mixed, n=150)

    for a, b in zip(r_legacy, r_mixed):
        assert a["final_weight"] == b["final_weight"]
        assert a["absorbed_weight"] == b["absorbed_weight"]
        assert a["final_energy"] == b["final_energy"]


def test_split_density_mixture_reproduces_single_isotope_physics():
    # Splitting one isotope into two half-density entries of the *same* isotope
    # must reproduce the single-isotope result (statistically — the extra
    # isotope-selection draw perturbs the RNG stream but not the physics).
    A = 207.2
    N = Material("Pb208", 11.35, 208, A).number_density
    n = 600

    pure = [Region(surfaces=[Sphere((0, 0, 0), 8.0)],
                   name="s", priority=1, element="Pb208")]
    split = [Region(surfaces=[Sphere((0, 0, 0), 8.0)], name="s", priority=1,
                    composition=[Nuclide("Pb208", N / 2, A, atomic_mass=208),
                                 Nuclide("Pb208", N / 2, A, atomic_mass=208)])]

    leak_pure = _leakage(_run(pure, n=n), n)
    leak_split = _leakage(_run(split, n=n), n)
    assert leak_split == pytest.approx(leak_pure, abs=0.03)


def test_two_isotope_mixture_runs_and_leaks_physically():
    # A genuine two-isotope mixture (heavy scatterer + light moderator) must
    # transport without error and leak a physical fraction in (0, 1].
    comp = Material.mixture("PbC", 11.0,
                            [("Pb208", 207.2, 0.5, 208),
                             ("C12", 11.9, 0.5, 12)])
    mediums = [Region(surfaces=[Sphere((0, 0, 0), 8.0)],
                      name="s", priority=1, composition=comp)]
    res = _run(mediums, energy=2.0e6, n=200)
    leak = _leakage(res, 200)
    assert 0.0 < leak <= 1.0
    assert all(np.isfinite(r["final_energy"]) for r in res)
