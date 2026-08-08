"""PyNeut vs OpenMC k-infinity for a BEAVRS HZP pin cell (reflective lattice).

This is the suite's thermal/moderated criticality check, complementing the
fast-metal ICSBEP spheres in validate_keff.py. The geometry is a single
Westinghouse 17x17 fuel pin (UO2 / He gap / Zircaloy-4 / borated water) in an
all-reflective cube, i.e. an infinite lattice, so the eigenvalue is k_infinity.

PyNeut has no S(alpha,beta) bound-thermal-scattering kernel and treats the
moderator with free-gas scattering at 294 K, whereas the realistic BEAVRS HZP
state is 600 K with c_H_in_H2O. To keep the z-score a fair test of the transport
kernel rather than of the thermal-scattering data, the validation runs *both*
codes in matched physics (free-gas hydrogen, 294 K, no S(a,b)). The realistic
600 K + S(a,b) OpenMC value is also reported, so the gap between it and the
matched value documents exactly the physics PyNeut omits.

Compositions are identical by construction: the BEAVRS materials are built once
in OpenMC and their per-nuclide atom densities are fed unchanged to PyNeut.

Performance note: an all-reflective lattice has *no* leakage, so a neutron only
terminates on absorption. Thermal neutrons therefore random-walk through the
moderator for hundreds of collisions per history, which makes the PyNeut side
far more expensive than the leakage-dominated fast-metal spheres in
validate_keff.py (order 100x more collisions per history). The OpenMC reference
is quick; budget accordingly for the pure-Python run, or use a smaller --n.

The defaults are the schedule reported in the paper: 3,000 neutrons/generation
in PyNeut against 8,000 in OpenMC, 140 generations discarding the first 40, and
five independent seeds per code. Because an all-reflective lattice has a
dominance ratio near unity, the single-run standard error understates the true
spread, so sigma is the spread across seeds rather than the batch-means error
of one run. At those settings the reference result is PyNeut 1.29297 +/-
0.00069 against OpenMC 1.29225 +/- 0.00069, a gap of +72 +/- 98 pcm (z = 0.74);
see pincell_reference.json, whose data_note records the nuclear-data snapshot
those numbers belong to.

Budget: the full default run is hours of pure-Python transport (see the note
above). Use --quick for a smoke test, or --seeds 1 --n 1000 for a spot check.

Run:  python validate_pincell.py --endf /path/to/endfb
      python validate_pincell.py --quick          # fast smoke test
"""
import argparse
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF_JSON = os.path.join(HERE, "pincell_reference.json")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))

from _common import ENDF_PATH as DEFAULT_ENDF, HAS_OPENMC  # noqa: E402
from src.cross_section_read import CrossSectionReader  # noqa: E402
from src.medium import Region, Cylinder, Box  # noqa: E402
from src.material import Nuclide  # noqa: E402
from src.vt_calc import VelocitySampler  # noqa: E402
from src.settings import Settings  # noqa: E402
from src.random_number_generator import RNGHandler  # noqa: E402
from src.criticality import run_keff  # noqa: E402
from src.physics import sample_watt_spectrum, watt_params_for  # noqa: E402

NEUTRON_MASS_AMU = 1.00866491588

# BEAVRS HZP pin-cell dimensions (cm) and state.
PITCH = 1.25984
R_PELLET = 0.39218
R_INNER_CLAD = 0.40005
R_OUTER_CLAD = 0.45720
T_HZP = 600.0


# ----------------------------------------------------------------- BEAVRS materials
def build_openmc_materials(temperature, thermal):
    """The four BEAVRS HZP pin-cell materials. `thermal=True` adds the
    c_H_in_H2O bound-scattering law (realistic); False leaves free-gas H."""
    import openmc

    fuel = openmc.Material(name='UO2_3.1')
    fuel.set_density('g/cm3', 10.24)
    fuel.add_element('U', 1.0, enrichment=3.1)
    fuel.add_element('O', 2.0)
    fuel.temperature = temperature

    gap = openmc.Material(name='Helium_Gap')
    gap.set_density('g/cm3', 0.0015)
    gap.add_element('He', 1.0)
    gap.temperature = temperature

    clad = openmc.Material(name='Zircaloy')
    clad.set_density('g/cm3', 6.55)
    clad.add_element('Zr', 0.98115, percent_type='wo')
    clad.add_element('Sn', 0.01450, percent_type='wo')
    clad.add_element('Fe', 0.00210, percent_type='wo')
    clad.add_element('Cr', 0.00100, percent_type='wo')
    clad.temperature = temperature

    water = openmc.Material(name='Borated_Water')
    water.set_density('g/cm3', 0.743)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_element('B', 9.75e-4)
    if thermal:
        water.add_s_alpha_beta('c_H_in_H2O')
    water.temperature = temperature

    return {'fuel': fuel, 'gap': gap, 'clad': clad, 'water': water}


def pyneut_composition(omc_material, temperature):
    """Convert an OpenMC material to a list of PyNeut Nuclides with identical
    atom densities (atoms/b-cm -> atoms/cm3)."""
    import openmc.data
    nuclides = []
    for name, dens_barn in omc_material.get_nuclide_atom_densities().items():
        mass = openmc.data.atomic_mass(name)            # amu
        awr = mass / NEUTRON_MASS_AMU
        nuclides.append(Nuclide(name, float(dens_barn) * 1e24, awr,
                                atomic_mass=mass, temperature=temperature))
    return nuclides


# ----------------------------------------------------------------- geometry / source
def pyneut_mediums(comps):
    """Concentric pin in an all-reflective cube (infinite lattice).

    The pin cylinders are infinite in z, so every region is clipped to the cube;
    this puts the reflective box planes on every region's boundary and prevents
    a neutron streaming along the axis from leaking out the (unbounded) cylinder.
    """
    half = PITCH / 2.0
    pellet = Cylinder('z', R_PELLET, (0, 0, 0))
    inner = Cylinder('z', R_INNER_CLAD, (0, 0, 0))
    outer = Cylinder('z', R_OUTER_CLAD, (0, 0, 0))
    cube = Box(-half, half, -half, half, -half, half, boundary_type='reflective')

    inside_inner = Region(surfaces=[inner, cube], operation='intersection')
    inside_outer = Region(surfaces=[outer, cube], operation='intersection')

    fuel = Region(surfaces=[pellet, cube], operation='intersection',
                  name='fuel', priority=4, composition=comps['fuel'])
    gap = Region(surfaces=[inside_inner, pellet], operation='difference',
                 name='gap', priority=3, composition=comps['gap'])
    clad = Region(surfaces=[inside_outer, inner], operation='difference',
                  name='clad', priority=2, composition=comps['clad'])
    water = Region(surfaces=[cube, outer], operation='difference', name='water',
                   priority=1, composition=comps['water'])
    return [fuel, gap, clad, water]


def pellet_source(n, rng, element='U235'):
    """Uniform fission source inside the fuel pellet (infinite axial extent)."""
    half = PITCH / 2.0
    a_w, b_w = watt_params_for(element)
    states = []
    for _ in range(n):
        r = R_PELLET * np.sqrt(rng.random())
        ang = 2 * np.pi * rng.random()
        states.append({
            "x": r * np.cos(ang), "y": r * np.sin(ang),
            "z": rng.uniform(-half, half),
            "theta": np.arccos(2 * rng.random() - 1),
            "phi": 2 * np.pi * rng.random(),
            "energy": sample_watt_spectrum(rng, a_w, b_w),
            "weight": 1.0, "has_interacted": False,
        })
    return states


# ----------------------------------------------------------------------- the two runs
def run_pyneut(endf, comps, n, gens, inactive, seed=1):
    reader = CrossSectionReader(endf, temperature="294K")
    mediums = pyneut_mediums(comps)
    src = pellet_source(n, RNGHandler(seed))
    settings = Settings("criticality", 1)
    from multiprocessing import Pool
    with Pool() as pool:
        out = run_keff(reader, mediums, src, 1.0, 1.0, VelocitySampler(1.0),
                       settings, generations=gens, inactive=inactive,
                       seed=seed, pool=pool)
    return out["k_eff"], out["k_sem"]


def run_openmc(materials, n, gens, inactive, temperature, tag, seed=None):
    import openmc
    half = PITCH / 2.0

    mats = openmc.Materials(list(materials.values()))

    pellet = openmc.ZCylinder(r=R_PELLET)
    inner = openmc.ZCylinder(r=R_INNER_CLAD)
    outer = openmc.ZCylinder(r=R_OUTER_CLAD)
    xl = openmc.XPlane(-half, boundary_type='reflective')
    xr = openmc.XPlane(half, boundary_type='reflective')
    yb = openmc.YPlane(-half, boundary_type='reflective')
    yt = openmc.YPlane(half, boundary_type='reflective')
    zb = openmc.ZPlane(-half, boundary_type='reflective')
    zt = openmc.ZPlane(half, boundary_type='reflective')
    box = +xl & -xr & +yb & -yt & +zb & -zt

    # clip every cell to the reflective cube (the cylinders are infinite in z)
    fuel_c = openmc.Cell(fill=materials['fuel'], region=-pellet & box)
    gap_c = openmc.Cell(fill=materials['gap'], region=+pellet & -inner & box)
    clad_c = openmc.Cell(fill=materials['clad'], region=+inner & -outer & box)
    water_c = openmc.Cell(fill=materials['water'], region=+outer & box)
    geom = openmc.Geometry([fuel_c, gap_c, clad_c, water_c])

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = gens
    settings.inactive = inactive
    settings.particles = n
    settings.temperature = {'method': 'nearest', 'default': temperature}
    settings.output = {'summary': False, 'tallies': False}
    settings.verbosity = 1
    if seed is not None:
        settings.seed = int(seed)
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-R_PELLET, -R_PELLET, -half),
                               (R_PELLET, R_PELLET, half)),
        constraints={'fissionable': True})

    run_dir = os.path.join(HERE, 'openmc_runs', f'omc_pincell_{tag}')
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(run_dir)
        model = openmc.Model(geom, mats, settings)
        sp_path = model.run(output=False)
        if sp_path is None:
            sp_path = sorted(glob.glob(os.path.join(run_dir, 'statepoint.*.h5')),
                             key=os.path.getmtime)[-1]
        with openmc.StatePoint(sp_path) as sp:
            k = sp.keff
        return float(k.nominal_value), float(k.std_dev)
    finally:
        os.chdir(cwd)


def main():
    ap = argparse.ArgumentParser(description="PyNeut vs OpenMC k-inf (BEAVRS pin cell)")
    ap.add_argument("--endf", default=os.environ.get("PYNEUT_PINCELL_ENDF", DEFAULT_ENDF),
                    help="ENDF/B HDF5 base dir (must hold the 36 pin-cell nuclides)")
    ap.add_argument("--n", type=int, default=3000,
                    help="PyNeut neutrons per generation (paper: 3000)")
    ap.add_argument("--n-omc", type=int, default=8000,
                    help="OpenMC neutrons per generation (paper: 8000)")
    ap.add_argument("--gens", type=int, default=140, help="total generations")
    ap.add_argument("--inactive", type=int, default=40, help="inactive generations")
    ap.add_argument("--seeds", type=int, default=5,
                    help="independent replications per code (paper: 5); the "
                         "reported sigma is the spread across seeds, not the "
                         "single-run standard error")
    ap.add_argument("--save-reference", action="store_true",
                    help="write the result to pincell_reference.json (the "
                         "artifact the paper cites)")
    ap.add_argument("--quick", action="store_true", help="fast smoke test")
    ap.add_argument("--no-realistic", action="store_true",
                    help="skip the 600 K + S(a,b) OpenMC reference run")
    args = ap.parse_args()

    if args.quick:
        args.n, args.n_omc, args.gens, args.inactive, args.seeds = (
            300, 800, 40, 10, 2)

    if not HAS_OPENMC:
        print("OpenMC not importable; cannot run the live comparison.")
        return

    xs_xml = os.path.join(args.endf, "cross_sections.xml")
    if os.path.exists(xs_xml):
        os.environ["OPENMC_CROSS_SECTIONS"] = xs_xml
    else:
        from _common import ensure_openmc_data
        ensure_openmc_data()

    print(f"\nPyNeut vs OpenMC k-infinity — BEAVRS HZP pin cell "
          f"(N_pyn={args.n}/gen, N_omc={args.n_omc}/gen, {args.gens} gens, "
          f"{args.inactive} inactive, {args.seeds} seeds/code)")
    print(f"data: {args.endf}\n")

    # --- matched physics: free-gas H, 294 K, no S(a,b), in BOTH codes ---
    # Independent replications: an all-reflective lattice has a dominance ratio
    # near unity, so the single-run standard error understates the true spread.
    # Both codes are therefore replicated over independent seeds and reported as
    # the seed mean +/- the standard error of that mean.
    matched_mats = build_openmc_materials(temperature=294.0, thermal=False)
    comps = {k: pyneut_composition(m, 294.0) for k, m in matched_mats.items()}

    pyn_k, omc_k = [], []
    for i in range(args.seeds):
        seed = i + 1
        k, _ = run_pyneut(args.endf, comps, args.n, args.gens, args.inactive,
                          seed=seed)
        pyn_k.append(k)
        # rebuild materials (the matched set is consumed for atom densities)
        mats = build_openmc_materials(temperature=294.0, thermal=False)
        ko_i, _ = run_openmc(mats, args.n_omc, args.gens, args.inactive, 294.0,
                             f"matched_seed{seed}", seed=seed)
        omc_k.append(ko_i)
        print(f"  seed {seed}: PyNeut {k:.5f}   OpenMC {ko_i:.5f}", flush=True)

    def mean_sem(v):
        a = np.asarray(v, dtype=float)
        return float(a.mean()), float(a.std(ddof=1) / np.sqrt(a.size)) if a.size > 1 else 0.0

    kp, sp = mean_sem(pyn_k)
    ko, so = mean_sem(omc_k)

    denom = np.hypot(sp, so)
    z = abs(ko - kp) / denom if denom > 0 else float("inf")
    status = "OK" if z <= 2 else ("WARNING" if z <= 4 else "FAIL")

    print("\nMatched physics (free-gas H, 294 K, no S(a,b)) — the validation result")
    print(f"  PyNeut k_inf : {kp:.5f} ± {sp:.5f}   (spread {np.std(pyn_k, ddof=1):.5f})")
    print(f"  OpenMC k_inf : {ko:.5f} ± {so:.5f}   (spread {np.std(omc_k, ddof=1):.5f})")
    print(f"  gap  : {1e5 * (kp - ko):+.0f} ± {1e5 * denom:.0f} pcm")
    print(f"  z = {z:.2f}   [{status}]   (OK <= 2, WARNING <= 4, FAIL > 4)\n")

    if args.save_reference:
        import json
        from datetime import date
        from _common import data_manifest
        ref = {
            "data_manifest": data_manifest(args.endf),
            "pyneut_5seeds": {str(i + 1): round(v, 5) for i, v in enumerate(pyn_k)},
            "pyneut_mean": round(kp, 5), "pyneut_sem": round(sp, 5),
            "openmc_i40_5seeds": [round(v, 5) for v in omc_k],
            "openmc_mean": round(ko, 5), "openmc_sem": round(so, 5),
            "gap_pcm": round(1e5 * (kp - ko), 1), "z": round(z, 2),
            "config": (f"{args.seeds}v{args.seeds} seeds, N_pyn={args.n} "
                       f"N_omc={args.n_omc}, gens={args.gens}, "
                       f"inactive={args.inactive}, matched physics"),
            "data_note": (
                f"Regenerated {date.today().isoformat()} by validate_pincell.py "
                f"--save-reference against the endfb snapshot then on disk. A "
                f"cross-code k_inf is a statement about a code pair AND a data "
                f"snapshot: this lattice moves by ~330 pcm across an evaluation "
                f"refresh, so values from different snapshots are not comparable."),
        }
        with open(REF_JSON, "w") as f:
            json.dump(ref, f, indent=1)
        print(f"  wrote {REF_JSON}\n")

    # --- realistic BEAVRS HZP (OpenMC only): documents the omitted physics ---
    # Optional, and skipped cleanly when the library has no S(a,b) sublibrary:
    # ensure_openmc_data() builds cross_sections.xml from endfb/neutron/*.h5,
    # which are neutron sublibrary files only, so c_H_in_H2O is absent unless
    # the user has pointed --endf at a full OpenMC data library. This must not
    # take down the validation result computed above.
    if not args.no_realistic:
        try:
            real_mats = build_openmc_materials(temperature=T_HZP, thermal=True)
            kr, sr = run_openmc(real_mats, args.n_omc, args.gens, args.inactive,
                                T_HZP, "realistic")
        except Exception as exc:
            first = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
            print("Realistic BEAVRS HZP (600 K + c_H_in_H2O) — SKIPPED")
            print(f"  {first}")
            print("  This run needs a data library carrying the thermal-scattering\n"
                  "  sublibrary (c_H_in_H2O). It is optional: the matched-physics\n"
                  "  comparison above is the validation result. Pass --no-realistic\n"
                  "  to skip it silently.\n")
        else:
            print("Realistic BEAVRS HZP (600 K + c_H_in_H2O, OpenMC) — documents the gap")
            print(f"  OpenMC k_inf : {kr:.5f} ± {sr:.5f}")
            print(f"  realistic - matched = {kr - ko:+.5f} dk  "
                  f"({1e5 * (kr - ko):+.0f} pcm): the bound-H + temperature physics "
                  f"PyNeut does not model.\n")


if __name__ == "__main__":
    main()
