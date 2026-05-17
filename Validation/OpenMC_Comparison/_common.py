"""
Shared utilities for PyNeut vs OpenMC validation scripts.

All validate_*.py scripts import from here. Handles path setup,
PyNeut simulation runner, statepoint reading, and result formatting.
"""
import os
import sys
import numpy as np
from multiprocessing import Pool

# ---------------------------------------------------------------------------
# Path resolution — works when scripts are run from anywhere
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
ENDF_PATH = os.path.join(PROJECT_ROOT, 'endfb')
STATEPOINT_DIR = os.path.join(HERE, 'benchmark_diverse')

sys.path.insert(0, PROJECT_ROOT)

from src.cross_section_read import CrossSectionReader
from src.vt_calc import VelocitySampler
from src.simulation import simulate_single_particle
from src.material import Material
from src.medium import Region, Sphere, Cylinder, Box, Plane
from src.random_number_generator import RNGHandler
from src.settings import Settings

try:
    import openmc
    HAS_OPENMC = True
except ImportError:
    HAS_OPENMC = False

# Lazy-loaded global reader — shared across all cases in a single run
_reader = None

def get_reader():
    global _reader
    if _reader is None:
        _reader = CrossSectionReader(ENDF_PATH)
    return _reader


# ---------------------------------------------------------------------------
# PyNeut runner
# ---------------------------------------------------------------------------

def run_pyneut(case, n_particles=10000):
    """
    Run a PyNeut shielding simulation for a test case dict.

    Case dict keys:
        el   : isotope name matching ENDF filename, e.g. 'Pb208'
        rho  : density [g/cm3]
        A    : atomic weight ratio
        E    : source energy [eV]
        geo  : 'sphere' | 'slab' | 'cyl'
        dims : [radius] for sphere, [half-thickness] for slab,
               [radius, half-height] for cyl (all in cm)

    Returns (leakage_fraction, avg_escape_energy_eV).
    """
    reader = get_reader()
    mat = Material(case['el'], case['rho'], case['A'], case['A'])
    N = mat.number_density
    A = case['A']

    geo = case['geo']
    if geo == 'sphere':
        R = case['dims'][0]
        mediums = [Region(
            surfaces=[Sphere((0, 0, 0), R)],
            operation='intersection', name='Sphere', priority=1, element=case['el'],
        )]
    elif geo == 'slab':
        d = case['dims'][0]
        mediums = [Region(
            surfaces=[Box(-d, d, -50, 50, -50, 50)],
            operation='intersection', name='Slab', priority=1, element=case['el'],
        )]
    elif geo == 'cyl':
        R, H = case['dims']
        mediums = [Region(
            surfaces=[Cylinder('z', R, (0, 0, 0)), Plane(0, 0, 1, H), Plane(0, 0, -1, H)],
            operation='intersection', name='Cyl', priority=1, element=case['el'],
        )]
    else:
        raise ValueError(f"Unknown geometry: {geo!r}")

    settings = Settings('shielding', n_particles)
    rngs = [RNGHandler(12345 + i) for i in range(n_particles)]
    states = [{
        'x': 0.0, 'y': 0.0, 'z': 0.0,
        'theta': np.arccos(1 - 2 * r.random()),
        'phi': r.uniform(0, 2 * np.pi),
        'has_interacted': False,
        'energy': case['E'],
        'weight': 1.0,
    } for r in rngs]

    try:
        reader.get_inelastic_components(case['el'], case['E'], N)
    except Exception:
        pass

    args = [
        (s, reader, mediums, A, N, VelocitySampler(mat.kg_mass), None, False, r, settings)
        for s, r in zip(states, rngs)
    ]

    with Pool() as pool:
        results = pool.map(simulate_single_particle, args)

    escaped = [(r['final_energy'], r['final_weight']) for r in results if r['final_weight'] > 0]
    total_w = sum(w for _, w in escaped)
    leakage = total_w / n_particles
    avg_e = sum(e * w for e, w in escaped) / total_w if total_w > 0 else 0.0
    return leakage, avg_e


# ---------------------------------------------------------------------------
# OpenMC statepoint reader
# ---------------------------------------------------------------------------

def read_statepoint(sp_path):
    """
    Read leakage fraction and average escape energy from a stored OpenMC
    statepoint.  Requires the openmc Python package.

    Expects tallies named 'leakage' and 'spectrum' (standard in all
    benchmark_diverse runs).

    Returns (leakage_fraction, avg_energy_eV) or (None, None) if openmc
    is not installed or reading fails.
    """
    if not HAS_OPENMC:
        return None, None
    try:
        with openmc.StatePoint(sp_path) as sp:
            leak = sp.get_tally(name='leakage')
            leak_val = float(np.sum(np.abs(leak.mean.flatten())))

            spec = sp.get_tally(name='spectrum')
            df = spec.get_pandas_dataframe()
            df['E_mid'] = (df['energy low [eV]'] + df['energy high [eV]']) / 2
            total_w = df['mean'].abs().sum()
            avg_e = float((df['E_mid'] * df['mean'].abs()).sum() / total_w) if total_w > 0 else 0.0
            return leak_val, avg_e
    except Exception as e:
        print(f'  [warn] statepoint read failed ({sp_path}): {e}')
        return None, None


# ---------------------------------------------------------------------------
# Result formatting
# ---------------------------------------------------------------------------

# Thresholds for pass/fail
_LEAKAGE_WARN = 1.0   # %
_LEAKAGE_FAIL = 3.0   # %
_ENERGY_WARN  = 2.0   # %
_ENERGY_FAIL  = 5.0   # %


def make_result(case_name, isotope, geometry, energy_mev,
                omc_leak, omc_energy, pyn_leak, pyn_energy,
                openmc_source='notebook'):
    """
    Build a result dict suitable for CSV output and terminal printing.

    openmc_source values:
        'notebook'   — hardcoded from prior notebook run
        'statepoint' — read live from benchmark_diverse statepoint
        'live'       — OpenMC was run on the fly
        'N/A'        — no OpenMC reference available
    """
    if omc_leak is not None and omc_energy is not None and omc_energy > 0:
        l_diff = abs(omc_leak - pyn_leak) / omc_leak * 100
        e_diff = abs(omc_energy - pyn_energy) / omc_energy * 100
        if l_diff <= _LEAKAGE_WARN and e_diff <= _ENERGY_WARN:
            status = 'OK'
        elif l_diff <= _LEAKAGE_FAIL and e_diff <= _ENERGY_FAIL:
            status = 'WARNING'
        else:
            status = 'FAIL'
    else:
        l_diff = e_diff = None
        status = 'NO_REF'

    return {
        'case':               case_name,
        'isotope':            isotope,
        'geometry':           geometry,
        'energy_MeV':         energy_mev,
        'openmc_leakage':     round(omc_leak, 5)  if omc_leak  is not None else None,
        'pyneut_leakage':     round(pyn_leak, 5),
        'leakage_diff_pct':   round(l_diff, 2)    if l_diff    is not None else None,
        'openmc_energy_eV':   round(omc_energy)   if omc_energy is not None else None,
        'pyneut_energy_eV':   round(pyn_energy),
        'energy_diff_pct':    round(e_diff, 2)    if e_diff    is not None else None,
        'openmc_source':      openmc_source,
        'status':             status,
    }


def print_result(r):
    tag = {'OK': '✓', 'WARNING': '!', 'FAIL': '✗', 'NO_REF': '?'}.get(r['status'], '?')
    print(f"\n  [{tag}] {r['case']}")
    print(f"      {'Metric':<22} {'OpenMC':>12}  {'PyNeut':>12}  {'Diff':>8}")
    print(f"      {'-'*58}")

    omc_l = f"{r['openmc_leakage']:.5f}"  if r['openmc_leakage']  is not None else 'N/A'
    omc_e = f"{r['openmc_energy_eV']:,}"  if r['openmc_energy_eV'] is not None else 'N/A'
    l_d   = f"{r['leakage_diff_pct']:.2f}%"  if r['leakage_diff_pct'] is not None else 'N/A'
    e_d   = f"{r['energy_diff_pct']:.2f}%"   if r['energy_diff_pct']  is not None else 'N/A'

    print(f"      {'Leakage fraction':<22} {omc_l:>12}  {r['pyneut_leakage']:>12.5f}  {l_d:>8}")
    print(f"      {'Avg escape energy (eV)':<22} {omc_e:>12}  {r['pyneut_energy_eV']:>12,}  {e_d:>8}")
    print(f"      source: {r['openmc_source']}")
