"""Smoke tests for the per-history diagnostic fields the validation suite relies
on (secondary-sampling counters and the uncollided-escape flag)."""
import os

import numpy as np
import pytest

from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Sphere, Box
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import simulate_single_particle
from src.vt_calc import VelocitySampler

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")
HAVE_DATA = os.path.exists(os.path.join(ENDF, "neutron", "Pb208.h5"))

pytestmark = pytest.mark.skipif(not HAVE_DATA, reason="ENDF data not present")


def _run(case_energy, geo, A=207.2, el="Pb208", rho=11.35, n=200, mode="shielding"):
    reader = CrossSectionReader(ENDF)
    mat = Material(el, rho, 208, A)
    N = mat.number_density
    settings = Settings(mode, n)
    rngs = [RNGHandler(100 + i) for i in range(n)]
    states = [{
        "x": 0.0, "y": 0.0, "z": 1e-9,
        "theta": geo["theta"], "phi": 0.0,
        "energy": case_energy, "weight": 1.0, "has_interacted": False,
    } for _ in rngs]
    args = [(s, reader, geo["mediums"], A, N, VelocitySampler(mat.kg_mass),
             None, False, r, settings) for s, r in zip(states, rngs)]
    return [simulate_single_particle(a) for a in args]


def test_result_has_diagnostic_keys():
    geo = {"theta": np.pi / 2,
           "mediums": [Region(surfaces=[Sphere((0, 0, 0), 8.0)],
                              name="s", priority=1, element="Pb208")]}
    res = _run(2.0e6, geo, n=50)
    for r in res:
        assert "physics_counts" in r and isinstance(r["physics_counts"], dict)
        assert "source_uncollided" in r and isinstance(r["source_uncollided"], bool)


def test_nxn_counters_fire_at_14mev():
    # A 14 MeV Pb208 sphere should produce (n,2n) events and log parent counters.
    geo = {"theta": np.pi / 2,
           "mediums": [Region(surfaces=[Sphere((0, 0, 0), 10.0)],
                              name="s", priority=1, element="Pb208")]}
    res = _run(14.0e6, geo, n=300)
    total = {}
    for r in res:
        for k, v in r["physics_counts"].items():
            total[k] = total.get(k, 0) + v
    assert total.get("nxn_events", 0) > 0
    # Parent source must be accounted for on every (n,xn) event.
    parent = sum(total.get(f"parent_{s}", 0) for s in ("law61", "law4", "maxwellian"))
    assert parent == total.get("nxn_events", 0)


def test_uncollided_escape_thin_vs_thick():
    # A pencil beam through a thin slab should leave many uncollided; a thick
    # slab essentially none.
    thin = {"theta": 0.0,
            "mediums": [Region(surfaces=[Box(-1e4, 1e4, -1e4, 1e4, 0.0, 0.2)],
                              name="s", priority=1, element="Pb208")]}
    thick = {"theta": 0.0,
             "mediums": [Region(surfaces=[Box(-1e4, 1e4, -1e4, 1e4, 0.0, 60.0)],
                               name="s", priority=1, element="Pb208")]}
    f_thin = np.mean([r["source_uncollided"] for r in _run(2.0e6, thin, n=300, mode="criticality")])
    f_thick = np.mean([r["source_uncollided"] for r in _run(2.0e6, thick, n=300, mode="criticality")])
    assert f_thin > 0.7
    assert f_thick < 0.05
