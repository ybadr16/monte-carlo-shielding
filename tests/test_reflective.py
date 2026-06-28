"""Reflective boundary condition: specular-reflection math, geometry crossing,
and an end-to-end kernel check that a reflective box leaks nothing.

These tests need no nuclear data: the kernel test uses a pure-absorber mock
cross section, so the only physics exercised is streaming, boundary reflection,
and analog absorption.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.geometry import reflect_direction, calculate_nearest_crossing
from src.medium import Box, Plane, Region, Sphere
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import run_particle_kernel
from src.vt_calc import VelocitySampler


class _AbsorberReader:
    """Constant pure-absorber cross section (no scattering, no fission)."""
    def __init__(self, sigma_a):
        self.sigma_a = sigma_a

    def get_cross_sections(self, element, energy, sampler, number_density, awr):
        # (sigma_elastic, sigma_inelastic, sigma_capture, sigma_fission, extra)
        return 0.0, 0.0, self.sigma_a, 0.0, None


# ---------------------------------------------------------------- reflection math
def test_reflect_normal_incidence_reverses():
    assert reflect_direction(1.0, 0.0, 0.0, 1.0, 0.0, 0.0) == (-1.0, 0.0, 0.0)


def test_reflect_flips_only_normal_component():
    d = (1 / math.sqrt(2), -1 / math.sqrt(2), 0.0)
    out = reflect_direction(*d, 0.0, 1.0, 0.0)  # plane with normal +y
    assert out[0] == pytest.approx(1 / math.sqrt(2))
    assert out[1] == pytest.approx(1 / math.sqrt(2))   # v flipped
    assert out[2] == pytest.approx(0.0)


def test_reflect_preserves_unit_norm():
    d = np.array([0.3, -0.6, 0.7421])
    d /= np.linalg.norm(d)
    n = np.array([0.0, 0.0, 1.0])
    out = np.array(reflect_direction(*d, *n))
    assert np.linalg.norm(out) == pytest.approx(1.0, abs=1e-9)


def test_normal_orientation_irrelevant():
    d = (0.5, 0.5, math.sqrt(0.5))
    a = reflect_direction(*d, 0.0, 0.0, 1.0)
    b = reflect_direction(*d, 0.0, 0.0, -1.0)   # flipped normal
    assert a == pytest.approx(b)


# --------------------------------------------------------------- geometry crossing
def test_crossing_reports_reflective_surface():
    box = Box(-1, 1, -1, 1, -1, 1, boundary_type="reflective")
    state = {"x": 0.0, "y": 0.0, "z": 0.0}
    pt, reg, dist, surf = calculate_nearest_crossing(state, [box], 1.0, 0.0, 0.0)
    assert dist == pytest.approx(1.0)
    assert surf.boundary_type == "reflective"
    # reflecting off this face turns +x around; nudged point stays inside the box
    nx, ny, nz = surf.normal(*pt)
    u, v, w = reflect_direction(1.0, 0.0, 0.0, nx, ny, nz)
    assert u == pytest.approx(-1.0)


def test_default_boundary_is_transmission():
    assert Plane(1, 0, 0, 0).boundary_type == "transmission"
    assert Box(-1, 1, -1, 1, -1, 1).surfaces[0].boundary_type == "transmission"
    assert Sphere((0, 0, 0), 1).boundary_type == "transmission"


# ----------------------------------------------------------------- kernel end-to-end
def _run_box(boundary_type, n=400, sigma_a=0.05, half=1.0, seed=7):
    """Fire n neutrons from the centre of a cube with the given wall BC and
    return (escaped_count, statuses, all_inside)."""
    box = Box(-half, half, -half, half, -half, half, boundary_type=boundary_type)
    box.priority = 1
    reader = _AbsorberReader(sigma_a)
    settings = Settings("criticality", 1)        # analog: clean absorption
    sampler = VelocitySampler(1.0)
    rng = RNGHandler(seed)

    escaped = 0
    statuses = []
    all_inside = True
    tol = 1e-3
    for i in range(n):
        state = {
            "x": 0.0, "y": 0.0, "z": 0.0,
            "theta": np.arccos(2 * rng.random() - 1),
            "phi": 2 * np.pi * rng.random(),
            "energy": 1.0e6, "weight": 1.0, "has_interacted": False,
        }
        status, _children, traj, _aw, _ac = run_particle_kernel(
            state, reader, [box], 1.0, 1.0, sampler,
            None, True, RNGHandler(seed + 1 + i), settings)
        statuses.append(status)
        if status == "escaped":
            escaped += 1
        if traj:
            for (x, y, z) in traj:
                if max(abs(x), abs(y), abs(z)) > half + tol:
                    all_inside = False
    return escaped, statuses, all_inside


def test_reflective_box_leaks_nothing():
    escaped, statuses, all_inside = _run_box("reflective")
    assert escaped == 0, f"reflective box should not leak; got {escaped} escapes"
    assert all(s == "absorbed" for s in statuses)
    assert all_inside, "every tracked point should remain inside the reflective box"


def test_transmission_box_leaks_some():
    # Control: with vacuum walls and a weak absorber, a sizeable fraction leaks.
    escaped, _statuses, _inside = _run_box("transmission")
    assert escaped > 0, "vacuum-walled box should leak at least some neutrons"
