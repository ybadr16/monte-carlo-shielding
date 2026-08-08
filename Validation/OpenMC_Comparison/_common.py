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

sys.path.insert(0, PROJECT_ROOT)

from src.cross_section_read import CrossSectionReader
from src.vt_calc import VelocitySampler
from src.simulation import simulate_single_particle
from src.material import Material, NEUTRON_MASS_AMU
from src.medium import Region, Sphere, Cylinder, Box, Plane
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.statistics import escape_statistics, ratio_mean_uncertainty

try:
    import openmc
    HAS_OPENMC = True
except ImportError:
    HAS_OPENMC = False

# Number of batches used for PyNeut batch-means uncertainties
N_BATCHES = 20

# Shared leakage-spectrum energy grid. Both codes estimate the average escape
# energy from this identical histogram (bin midpoints), so any histogram bias
# cancels and the comparison reflects genuine physics differences, not the
# choice of binning.
E_NBINS = 500


def energy_grid(E):
    """Common leakage-spectrum bin edges for a source energy E (eV).

    The top is at least 1 eV so that thermal *upscatter* (neutrons gaining
    energy from the motion of the target nuclei, which carries them above the
    source energy) is captured rather than truncated.
    """
    top = max(1.01 * E, 1.0)
    return np.linspace(0.0, top, E_NBINS + 1)


# ---------------------------------------------------------------------------
# OpenMC nuclear-data library setup
# ---------------------------------------------------------------------------
# The endfb/neutron/*.h5 files are already in OpenMC HDF5 format.  We build a
# cross_sections.xml on the fly so OpenMC can run live without any external
# data library, then point OPENMC_CROSS_SECTIONS at it.  This makes the whole
# validation suite self-contained and reproducible on any machine that has the
# repo + openmc installed.
XS_XML = os.path.join(ENDF_PATH, 'cross_sections.xml')


def ensure_openmc_data():
    """Build endfb/cross_sections.xml from the local HDF5 files if needed and
    export OPENMC_CROSS_SECTIONS so live OpenMC runs find the data."""
    if not HAS_OPENMC:
        return False
    if not os.path.exists(XS_XML):
        import glob
        import openmc.data
        lib = openmc.data.DataLibrary()
        for fp in sorted(glob.glob(os.path.join(ENDF_PATH, 'neutron', '*.h5'))):
            try:
                lib.register_file(fp)
            except Exception as e:
                print(f'  [warn] could not register {os.path.basename(fp)}: {e}')
        lib.export_to_xml(XS_XML)
    os.environ['OPENMC_CROSS_SECTIONS'] = XS_XML
    return True

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

def data_manifest(endf_path=None):
    """Fingerprint of the nuclear-data directory a result was measured against.

    A cross-code comparison is a statement about a code pair AND a snapshot of
    the evaluation files: refreshing endfb/neutron moved the BEAVRS pin cell by
    a few hundred pcm in both codes, which once masqueraded as a transport bias.
    Every recorded reference therefore carries this digest, so two numbers can
    be checked for comparability before they are differenced.

    Returns {'n_files', 'total_bytes', 'sha256'}, the digest being over the
    sorted (name, size) pairs of endfb/neutron/*.h5 -- cheap, and sensitive to
    any file being added, removed or replaced.
    """
    import glob
    import hashlib
    base = endf_path or ENDF_PATH
    files = sorted(glob.glob(os.path.join(base, 'neutron', '*.h5')))
    h = hashlib.sha256()
    total = 0
    for f in files:
        size = os.path.getsize(f)
        total += size
        h.update(f"{os.path.basename(f)}:{size}\n".encode())
    return {'n_files': len(files), 'total_bytes': total,
            'sha256': h.hexdigest()[:16]}


def _assert_awr_matches_data(element, A, tol=1e-3):
    """Fail loudly if a case's A is not the data file's atomic-weight ratio.

    Both codes must see the same target mass; see the note in run_pyneut.
    """
    import h5py
    path = os.path.join(ENDF_PATH, 'neutron', f'{element}.h5')
    if not os.path.exists(path):
        return
    with h5py.File(path, 'r') as f:
        grp = f[element] if element in f else f[f'neutron/{element}']
        awr = float(grp.attrs['atomic_weight_ratio'])
    if abs(A / awr - 1.0) > tol:
        raise ValueError(
            f"{element}: case A={A} is {100 * (A / awr - 1):+.2f}% off the data "
            f"file's atomic_weight_ratio={awr:.5f}. A must be the AWR "
            f"(mass/neutron mass), not the molar mass in g/mol — otherwise "
            f"PyNeut and OpenMC run different masses and number densities.")


def run_pyneut(case, n_particles=10000, n_batches=N_BATCHES):
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

    Returns a dict with batch-means statistics:
        {leakage, leakage_sem, avg_energy, avg_energy_sem}.
    """
    reader = get_reader()
    # case['A'] is the atomic-weight ratio (mass / neutron mass). The molar mass
    # in g/mol is A * m_neutron, NOT A: passing A as the mass makes the free-gas
    # sampler's target ~0.9% too light and its number density ~0.9% too high,
    # which for a light moderator biases the thermal escape spectrum hot. OpenMC
    # uses the true isotopic mass, so the mass must be A * NEUTRON_MASS_AMU here.
    #
    # For the same reason case['A'] must be the file's AWR and not the element's
    # molar mass in g/mol: OpenMC builds its material from rho plus its own
    # atomic mass, so a case that supplies (say) 55.845 for Fe56 instead of
    # 55.4544 hands the two codes different kinematic masses AND different
    # number densities, which is a model difference rather than the kernel
    # difference the suite is meant to measure. Six case files carried the molar
    # mass until 2026-08-02; the check below stops that from recurring.
    _assert_awr_matches_data(case['el'], case['A'])
    mat = Material(case['el'], case['rho'],
                   case['A'] * NEUTRON_MASS_AMU, case['A'])
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

    # Batch-means uncertainties: one (leaked weight, escape energy) per history,
    # including non-leakers (weight 0) so the leakage denominator is correct.
    leaked_w = np.array([r['final_weight'] for r in results], dtype=float)
    esc_e = np.array([r['final_energy'] for r in results], dtype=float)
    stats = escape_statistics(leaked_w, esc_e, n_particles, n_batches,
                              e_bins=energy_grid(case['E']))

    # Aggregate the secondary-sampling diagnostic counters across histories.
    counts = {}
    for r in results:
        for k, v in r.get('physics_counts', {}).items():
            counts[k] = counts.get(k, 0) + v
    stats['physics_counts'] = counts

    # Per-history escape spectrum (same grid as OpenMC) for shape comparison,
    # with batch-means per-bin uncertainties.
    stats['spectrum'], stats['spectrum_sem'] = _pyneut_spectrum(
        leaked_w, esc_e, energy_grid(case['E']), n_batches)
    return stats


def _pyneut_spectrum(leaked_w, esc_e, e_bins, n_batches):
    """Batch-means leakage spectrum (current per energy bin) and its SEM."""
    from src.statistics import batch_edges, mean_sem
    edges = batch_edges(len(leaked_w), n_batches)
    per_batch = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        hist, _ = np.histogram(esc_e[lo:hi], bins=e_bins, weights=leaked_w[lo:hi])
        per_batch.append(hist / (hi - lo))      # normalise per source history
    per_batch = np.array(per_batch)
    nb = per_batch.shape[0]
    mean = per_batch.mean(axis=0)
    sem = per_batch.std(axis=0, ddof=1) / np.sqrt(nb) if nb > 1 else np.full_like(mean, np.nan)
    return mean, sem


# ---------------------------------------------------------------------------
# OpenMC live runner — geometry mirrors run_pyneut() exactly
# ---------------------------------------------------------------------------

def run_openmc(case, n_particles=10000, batches=20):
    """
    Run a live OpenMC fixed-source simulation for the same case dict used by
    run_pyneut().  Source, geometry and the leakage definition (current through
    every bounding surface) are matched so the two codes are directly
    comparable.

    Returns a dict {leakage, leakage_sem, avg_energy, avg_energy_sem} with
    OpenMC's own batch statistics (the spectrum SEM is obtained by ratio error
    propagation over the energy bins), or None if OpenMC/its data are
    unavailable.
    """
    if not ensure_openmc_data():
        return None

    el, rho, A, E, geo, dims = (
        case['el'], case['rho'], case['A'], case['E'], case['geo'], case['dims'])

    mat = openmc.Material()
    try:
        mat.add_nuclide(el, 1.0)
    except Exception:
        mat.add_element(el.rstrip('0123456789'), 1.0)
    mat.set_density('g/cm3', rho)
    materials = openmc.Materials([mat])

    if geo == 'sphere':
        s = openmc.Sphere(r=dims[0], boundary_type='vacuum')
        cell = openmc.Cell(region=-s, fill=mat)
        surfs = [s]
    elif geo == 'slab':
        d = dims[0]
        xl = openmc.XPlane(-d, boundary_type='vacuum')
        xr = openmc.XPlane(+d, boundary_type='vacuum')
        yl = openmc.YPlane(-50, boundary_type='vacuum')
        yr = openmc.YPlane(+50, boundary_type='vacuum')
        zl = openmc.ZPlane(-50, boundary_type='vacuum')
        zr = openmc.ZPlane(+50, boundary_type='vacuum')
        cell = openmc.Cell(region=+xl & -xr & +yl & -yr & +zl & -zr, fill=mat)
        surfs = [xl, xr, yl, yr, zl, zr]
    elif geo == 'cyl':
        R, H = dims
        cyl = openmc.ZCylinder(r=R, boundary_type='vacuum')
        bot = openmc.ZPlane(-H, boundary_type='vacuum')
        top = openmc.ZPlane(+H, boundary_type='vacuum')
        cell = openmc.Cell(region=-cyl & +bot & -top, fill=mat)
        surfs = [cyl, bot, top]
    else:
        raise ValueError(f'Unknown geometry: {geo!r}')

    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = batches
    settings.particles = max(1, n_particles // batches)
    settings.output = {'summary': False, 'tallies': False}
    settings.verbosity = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([E], [1.0]),
    )

    t_leak = openmc.Tally(name='leakage')
    t_leak.filters = [openmc.SurfaceFilter(surfs)]
    t_leak.scores = ['current']

    e_bins = energy_grid(E)
    e_mid = 0.5 * (e_bins[:-1] + e_bins[1:])
    t_spec = openmc.Tally(name='spectrum')
    t_spec.filters = [openmc.SurfaceFilter(surfs), openmc.EnergyFilter(e_bins)]
    t_spec.scores = ['current']
    tallies = openmc.Tallies([t_leak, t_spec])

    run_dir = os.path.join(HERE, 'openmc_runs', f"omc_{case['name']}")
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geometry, materials, settings, tallies)
        sp_path = model.run(output=False)
        if sp_path is None:
            # model.run() occasionally returns None even though the statepoint
            # was written; fall back to the newest statepoint in the run dir.
            import glob
            sps = sorted(glob.glob(os.path.join(run_dir, 'statepoint.*.h5')),
                         key=os.path.getmtime)
            if not sps:
                raise RuntimeError('OpenMC produced no statepoint')
            sp_path = sps[-1]
        with openmc.StatePoint(sp_path) as sp:
            # --- Leakage: sum over surfaces; SEMs add in quadrature ---
            lt = sp.get_tally(name='leakage')
            leak = float(np.abs(lt.mean.flatten()).sum())
            leak_sem = float(np.sqrt((lt.std_dev.flatten() ** 2).sum()))

            # --- Spectrum: ordered [surface, energy]; collapse the surface axis
            # (means add, SEMs in quadrature), then propagate the ratio
            # uncertainty of the weighted-average escape energy.
            # (get_pandas_dataframe is avoided — it breaks on SurfaceFilter with
            # the installed pandas/pyarrow backend.)
            st = sp.get_tally(name='spectrum')
            mean = np.abs(st.mean.flatten()).reshape(len(surfs), len(e_mid))
            sd = st.std_dev.flatten().reshape(len(surfs), len(e_mid))
            bin_w = mean.sum(axis=0)
            bin_sd = np.sqrt((sd ** 2).sum(axis=0))
            avg_e, avg_e_sem = ratio_mean_uncertainty(e_mid, bin_w, bin_sd)
    except Exception as e:
        print(f'  [warn] OpenMC run failed for {case["name"]}: {e}')
        return None
    finally:
        os.chdir(cwd)
    return {
        'leakage': leak,
        'leakage_sem': leak_sem,
        'avg_energy': avg_e,
        'avg_energy_sem': avg_e_sem,
        'spectrum': bin_w,           # leakage current per source neutron, per bin
        'spectrum_sem': bin_sd,
    }


# ---------------------------------------------------------------------------
# Result formatting
# ---------------------------------------------------------------------------
#
# Grading is statistical: the difference between codes is expressed in units of
# the combined standard error,  z = |OMC - PyN| / sqrt(sem_OMC^2 + sem_PyN^2).
# Two codes are "consistent" when z is small regardless of the raw % difference.
#   OK       : both metrics within 2 sigma
#   WARNING  : both within 4 sigma
#   FAIL     : either exceeds 4 sigma   (a statistically significant bias)
_Z_OK   = 2.0
_Z_WARN = 4.0


def compare_spectra(omc, pyn, n_groups=20):
    """
    Compare two leakage spectra (fine bins) by coarse-grouping and computing a
    chi-square goodness-of-fit of the *shape*:

        chi2 = Σ_g (P_g - O_g)^2 / (σ_Pg^2 + σ_Og^2)

    over groups that carry signal.  Returns (chi2, ndof, chi2_per_dof) or
    (None, None, None) if either spectrum is missing.
    """
    if omc is None or 'spectrum' not in omc or 'spectrum' not in pyn:
        return None, None, None
    O, sO = np.asarray(omc['spectrum']), np.asarray(omc['spectrum_sem'])
    P, sP = np.asarray(pyn['spectrum']), np.asarray(pyn['spectrum_sem'])
    L = len(O)
    edges = np.linspace(0, L, n_groups + 1).astype(int)

    chi2 = 0.0
    ndof = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        Og = O[lo:hi].sum()
        Pg = P[lo:hi].sum()
        var = np.nansum(sO[lo:hi] ** 2) + np.nansum(sP[lo:hi] ** 2)
        # Only score groups that carry signal and have a usable variance.
        if var <= 0 or (Og + Pg) <= 0:
            continue
        chi2 += (Pg - Og) ** 2 / var
        ndof += 1
    if ndof == 0:
        return None, None, None
    return chi2, ndof, chi2 / ndof


def _fallback_summary(physics_counts):
    """Derive per-(n,xn)-event fallback fractions from raw counters."""
    c = physics_counts or {}
    nxn = c.get('nxn_events', 0)
    out = {
        'nxn_events': nxn,
        'n3n_events': c.get('n3n_events', 0),
        'parent_maxwellian_pct': None,
        'child_maxwellian_pct': None,
    }
    if nxn > 0:
        out['parent_maxwellian_pct'] = round(100 * c.get('parent_maxwellian', 0) / nxn, 1)
    n2n_children = (c.get('child_law61', 0) + c.get('child_law4', 0)
                    + c.get('child_maxwellian', 0))
    if n2n_children > 0:
        out['child_maxwellian_pct'] = round(100 * c.get('child_maxwellian', 0) / n2n_children, 1)
    return out


def _zscore(omc, omc_sem, pyn, pyn_sem):
    """|difference| in units of the combined standard error, or None if the
    combined SEM is not finite/positive (cannot judge significance)."""
    diff = abs(omc - pyn)
    comb = np.hypot(omc_sem if np.isfinite(omc_sem) else 0.0,
                    pyn_sem if np.isfinite(pyn_sem) else 0.0)
    if not np.isfinite(comb) or comb <= 0:
        return None
    return diff / comb


def make_result(case_name, isotope, geometry, energy_mev,
                omc, pyn, openmc_source='live'):
    """
    Build a result dict (for CSV + terminal) from two statistics dicts
    {leakage, leakage_sem, avg_energy, avg_energy_sem}.  ``omc`` may be None
    when no OpenMC reference is available.
    """
    def pct(a, b):
        return abs(a - b) / a * 100 if a else None

    def e_round(x):
        # Keep sub-eV thermal energies meaningful instead of rounding to 0.
        if x is None:
            return None
        return round(x, 4) if abs(x) < 100 else round(x)

    if omc is not None and omc['avg_energy'] > 0:
        l_diff = pct(omc['leakage'], pyn['leakage'])
        e_diff = pct(omc['avg_energy'], pyn['avg_energy'])
        z_leak = _zscore(omc['leakage'], omc['leakage_sem'],
                         pyn['leakage'], pyn['leakage_sem'])
        z_energy = _zscore(omc['avg_energy'], omc['avg_energy_sem'],
                           pyn['avg_energy'], pyn['avg_energy_sem'])
        zs = [z for z in (z_leak, z_energy) if z is not None]
        z_max = max(zs) if zs else None
        if z_max is None:                       # no usable uncertainty
            status = 'NO_SIGMA'
        elif z_max <= _Z_OK:
            status = 'OK'
        elif z_max <= _Z_WARN:
            status = 'WARNING'
        else:
            status = 'FAIL'
    else:
        l_diff = e_diff = z_leak = z_energy = None
        status = 'NO_REF'

    _, _, chi2_dof = compare_spectra(omc, pyn)
    fb = _fallback_summary(pyn.get('physics_counts'))

    # Separate verdict on the *shape* of the leakage spectrum (the integral
    # status above only judges leakage and mean energy).
    if chi2_dof is None:
        spectrum_status = None
    elif chi2_dof <= 2.0:
        spectrum_status = 'GOOD'
    elif chi2_dof <= 5.0:
        spectrum_status = 'FAIR'
    else:
        spectrum_status = 'POOR'

    return {
        'case':               case_name,
        'isotope':            isotope,
        'geometry':           geometry,
        'energy_MeV':         energy_mev,
        'openmc_leakage':     round(omc['leakage'], 5) if omc else None,
        'openmc_leakage_sem': round(omc['leakage_sem'], 5) if omc else None,
        'pyneut_leakage':     round(pyn['leakage'], 5),
        'pyneut_leakage_sem': round(pyn['leakage_sem'], 5),
        'leakage_diff_pct':   round(l_diff, 2) if l_diff is not None else None,
        'leakage_zscore':     round(z_leak, 2) if z_leak is not None else None,
        'openmc_energy_eV':   e_round(omc['avg_energy']) if omc else None,
        'openmc_energy_sem':  e_round(omc['avg_energy_sem']) if omc else None,
        'pyneut_energy_eV':   e_round(pyn['avg_energy']),
        'pyneut_energy_sem':  e_round(pyn['avg_energy_sem']),
        'energy_diff_pct':    round(e_diff, 2) if e_diff is not None else None,
        'energy_zscore':      round(z_energy, 2) if z_energy is not None else None,
        'spectrum_chi2_dof':  round(chi2_dof, 2) if chi2_dof is not None else None,
        'spectrum_status':    spectrum_status,
        'nxn_events':         fb['nxn_events'],
        'parent_maxwellian_pct': fb['parent_maxwellian_pct'],
        'child_maxwellian_pct':  fb['child_maxwellian_pct'],
        'n3n_events':         fb['n3n_events'],
        'openmc_source':      openmc_source,
        'status':             status,
    }


def validate_case(case, n_particles=10000, isotope=None, n_batches=N_BATCHES):
    """
    Run both codes for a case dict and return a result row.  OpenMC is run
    live against the local data library, so every reference is reproducible
    and uses identical geometry/source to PyNeut.
    """
    pyn = run_pyneut(case, n_particles, n_batches)
    omc = run_openmc(case, n_particles, batches=n_batches)
    return make_result(
        case_name     = case['name'],
        isotope       = isotope or case['el'],
        geometry      = case['geo'],
        energy_mev    = case['E'] / 1e6,
        omc           = omc,
        pyn           = pyn,
        openmc_source = 'live' if omc is not None else 'N/A',
    )


def print_result(r):
    tag = {'OK': '✓', 'WARNING': '!', 'FAIL': '✗',
           'NO_REF': '?', 'NO_SIGMA': '?'}.get(r['status'], '?')
    print(f"\n  [{tag}] {r['case']}   ({r['status']})")
    print(f"      {'Metric':<22} {'OpenMC':>18}  {'PyNeut':>18}  {'Diff':>7}  {'z':>5}")
    print(f"      {'-'*74}")

    if r['openmc_leakage'] is not None:
        omc_l = f"{r['openmc_leakage']:.5f} ± {r['openmc_leakage_sem']:.5f}"
        omc_e = f"{r['openmc_energy_eV']:,} ± {r['openmc_energy_sem']:,}"
    else:
        omc_l = omc_e = 'N/A'
    pyn_l = f"{r['pyneut_leakage']:.5f} ± {r['pyneut_leakage_sem']:.5f}"
    pyn_e = f"{r['pyneut_energy_eV']:,} ± {r['pyneut_energy_sem']:,}"
    l_d = f"{r['leakage_diff_pct']:.2f}%" if r['leakage_diff_pct'] is not None else 'N/A'
    e_d = f"{r['energy_diff_pct']:.2f}%" if r['energy_diff_pct'] is not None else 'N/A'
    l_z = f"{r['leakage_zscore']:.1f}" if r['leakage_zscore'] is not None else '-'
    e_z = f"{r['energy_zscore']:.1f}" if r['energy_zscore'] is not None else '-'

    print(f"      {'Leakage fraction':<22} {omc_l:>18}  {pyn_l:>18}  {l_d:>7}  {l_z:>5}")
    print(f"      {'Avg escape energy (eV)':<22} {omc_e:>18}  {pyn_e:>18}  {e_d:>7}  {e_z:>5}")
    if r.get('spectrum_chi2_dof') is not None:
        print(f"      Spectrum shape: χ²/dof = {r['spectrum_chi2_dof']:.2f}")
    if r.get('nxn_events'):
        pm = r.get('parent_maxwellian_pct')
        cm = r.get('child_maxwellian_pct')
        msg = f"      (n,xn) events: {r['nxn_events']:,}"
        if pm is not None:
            msg += f"  |  parent Maxwellian fallback: {pm:.1f}%"
        if cm is not None:
            msg += f"  |  (n,2n) child: {cm:.1f}%"
        print(msg)
        if r.get('n3n_events'):
            print(f"      (n,3n) events: {r['n3n_events']:,}")
    print(f"      source: {r['openmc_source']}")
