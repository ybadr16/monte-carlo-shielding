<p align="center">
  <img src="docs/assets/pyneut-logo.jpg" alt="PyNeut logo" width="420">
</p>

# PyNeut — Python Neutronics

A Python Monte Carlo code for continuous-energy neutron transport through
user-defined materials and geometries, driven by ENDF/B-VIII nuclear data.
**Intended for educational and research-prototyping use only.**

[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Validated vs OpenMC](https://img.shields.io/badge/validated%20vs-OpenMC-success.svg)](#validation-status)

<p align="center">
  <img src="docs/images/validation_summary.png"
       alt="Every PyNeut validation case plotted against OpenMC as a z-score; twelve of thirteen fall inside the two-sigma agreement band"
       width="88%">
</p>

<p align="center">
  <img src="docs/images/xs_verification.png" alt="U235 cross sections read by PyNeut overlaid on OpenMC" width="49%">
  <img src="docs/images/leakage_spectra.png" alt="Leakage spectra from PyNeut and OpenMC" width="49%">
</p>

## Overview

PyNeut implements continuous-energy neutron physics, including:
- Elastic scattering with tabulated ENDF angular distributions and free-gas thermal motion
- Inelastic scattering — discrete levels and continuum / multiplication, (n,2n), (n,3n)
- Energy-dependent cross sections read directly from ENDF/B-VIII HDF5 files
- Multi-region constructive-solid-geometry (planes, spheres, cylinders, boxes, voids, priorities)
- Single-isotope or multi-isotope material mixtures per region (reacting isotope sampled per collision)
- Analog and survival-biased (implicit-capture + roulette) transport
- k-eigenvalue (criticality) via a fission-source power iteration with ν and a Watt χ spectrum
- Parallel particle tracking via `multiprocessing`, mesh tallies (VTK) and trajectory recording

<p align="center">
  <img src="docs/images/geometry_pincell.png"
       alt="Material map of the BEAVRS pin cell: UO2 fuel, helium gap, Zircaloy-4 clad and borated water"
       width="66%">
</p>

Geometry is sampled with the same containment test the tracker uses, so a
material map shows exactly what the transport sees — priorities included.
Render one for any model with `tools/plot_geometry.py`.

## Requirements

- Python 3.x
- Runtime: `numpy`, `h5py` (see `requirements.txt`)
- Validation suite: `openmc`; geometry tests: `trimesh`, `rtree`

## Installation

```bash
pip install -e .
# or
pip install -r requirements.txt
```

## Quick Start

```python
import numpy as np
from multiprocessing import Pool

from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Sphere
from src.simulation import simulate_single_particle
from src.vt_calc import VelocitySampler
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.tally import Tally

reader = CrossSectionReader("./endfb")                       # ENDF/B-VIII data
lead = Material("Lead", 11.35, 208, atomic_weight_ratio=207.2)
sampler = VelocitySampler(lead.kg_mass)

regions = [Region(surfaces=[Sphere((0, 0, 0), 10.0)],
                  name="Lead", priority=1, element="Pb208")]

settings = Settings(mode="shielding", particles=10_000)      # implicit capture
N = settings.num_particles
rngs = [RNGHandler(seed=12345 + i) for i in range(N)]

states = [{
    "x": 0.0, "y": 0.0, "z": 0.0,
    "theta": np.arccos(1 - 2 * r.random()),                  # isotropic
    "phi": r.uniform(0, 2 * np.pi),
    "energy": 2.0e6,                                          # eV
    "weight": 1.0,                                            # required
    "has_interacted": False,
} for r in rngs]

# NOTE the trailing `settings` element in each args tuple.
args = [(s, reader, regions, lead.atomic_weight_ratio, lead.number_density,
         sampler, None, False, r, settings) for s, r in zip(states, rngs)]

with Pool() as pool:
    results = pool.map(simulate_single_particle, args)

tally = Tally()
for res in results:
    tally.merge_partial_results(res)
tally.print_summary(N)
```

See `main.py` for a complete example with a mesh tally, and `docs/QUICK_START.md`
for slab/cylinder/detector variations.

## ENDF/B-VIII Data

PyNeut reads OpenMC-format HDF5 files from `./endfb/neutron/`. The data files are
**not shipped in the repository** (they are large); download ENDF/B-VIII.0 from
the [OpenMC data libraries](https://openmc.org/data-libraries/) or
[NNDC ENDF/B-VIII](https://www.nndc.bnl.gov/endf/) and place the `*.h5` files
there. The validation suite uses ten isotopes — Al27, B10, Be9, C12, Cd112, Fe56,
O16, Pb207, Pb208, U235 — plus U233/U234/U238 and Pu239–241/Ga for the ICSBEP
benchmarks.

## Code Structure

```
src/
├── cross_section_read.py   # ENDF reader: XS, angular & secondary-energy distributions
├── angular_distribution.py # Tabulated elastic (MT=2) angular sampling
├── material.py             # Material number density, atom mass, Nuclide & mixtures
├── medium.py               # CSG primitives (Region, Plane, Cylinder, Sphere, Box)
├── geometry.py             # Ray–surface boundary search
├── physics.py              # Elastic/inelastic kinematics, evaporation, Watt χ
├── vt_calc.py              # Free-gas target velocity sampling
├── settings.py             # Transport mode (analog vs implicit capture)
├── simulation.py           # Transport loop + secondary-neutron & fission bank
├── criticality.py          # k-eigenvalue fission-source power iteration
├── tally.py                # Scalar results accumulation
├── mesh.py                 # Cartesian mesh tally → VTK
└── random_number_generator.py
```

## Validation Status

PyNeut is validated against **OpenMC 0.15** (live, same ENDF/B data, matched
geometry/source), against **exact analytic results**, and at the
**cross-section-reader** level.

Agreement is graded statistically by the z-score
`z = |OpenMC − PyNeut| / √(σ²_OMC + σ²_PyN)` (**OK ≤ 2σ, WARNING ≤ 4σ,
FAIL > 4σ**). The average escape energy uses the same energy histogram in both
codes (binning bias cancels), and the full leakage **spectrum shape** is graded
separately by a coarse-group χ²/dof (GOOD ≤ 2, FAIR ≤ 5, POOR > 5).

| Case | E | Regime | Leak z | E z | Integral | Spectrum |
|------|--:|--------|------:|----:|:--------:|:--------:|
| Pb208 elastic (sph R=10) | 2 MeV | fast | 0.9 | 0.3 | OK | GOOD (0.90) |
| Pb208 (sph R=10) | 5 MeV | fast | 1.1 | 0.1 | OK | GOOD (0.65) |
| Pb208 (n,2n) (sph R=10) | 14 MeV | fast (n,2n) | 0.4 | 0.6 | OK | POOR (29) |
| Fe56 slab (d=2 cm) | 14 MeV | fast (n,xn) | 2.5 | 1.4 | WARN | POOR (5.5) |
| Fe56 sphere (R=10) | 1 MeV | fast | 0.8 | 0.4 | OK | GOOD (1.0) |
| Fe56 sphere (R=15) | 25 keV | intermediate | 0.1 | 0.6 | OK | GOOD (0.85) |
| Be9 (n,2n) (sph R=10) | 14 MeV | fast (n,2n) | 0.9 | 0.1 | OK | POOR (46) |
| Al27 slab (d=0.5 cm) | 1 MeV | fast | 0.0 | 0.2 | OK | GOOD (0.76) |
| C12 cyl (R=10,H=20) | 2 MeV | fast | 0.0 | 0.4 | OK | GOOD (0.78) |
| O16 sphere (R=20) | 10 MeV | fast | 0.4 | 1.6 | OK | GOOD (1.0) |
| B10 sphere (R=1) | 1 keV | epithermal | 0.6 | 1.2 | OK | POOR (15) |
| Cd112 sphere (R=2) | 0.0253 eV | thermal | 0.7 | 0.3 | OK | GOOD (0.08) |
| C12 sphere (R=5) | 0.0253 eV | thermal | 1.2 | 0.4 | OK | GOOD (0.47) |

**12 OK, 1 warning, 0 failures.** Fast-neutron transport (the intended domain)
agrees with OpenMC within ~2σ on leakage, with GOOD spectrum shape on the
elastic, discrete-inelastic and resonance cases. **Thermal scattering is
validated** — free-gas vector kinematics (below 10 eV) let neutrons upscatter
and equilibrate to the 294 K Maxwellian (before that fix c12 was 86σ off), and
the earlier ~3 % hot bias in the `c12_thermal` mean escape energy has been
closed: it came from three compounding errors (re-broadening an already
Doppler-broadened flight cross section, the SVT branch probability using the
raw neutron speed instead of the dimensionless ratio βv, and the sampler mass
taken as the AWR rather than A·mₙ). Both thermal cases now agree within 0.4σ.

Every case takes the target mass from the data file's atomic-weight ratio, and
`_common.py` refuses a case whose declared `A` disagrees with it — six case
files (ten of the thirteen cases) previously supplied the element's molar mass
in g/mol, which handed the
two codes different masses *and* number densities (0.5–1.0 %) and produced a
spurious 2.1σ escape-energy warning on `c12_cyl_2mev`.

<p align="center">
  <img src="docs/images/cases/spectrum_pb208_n2n.png"
       alt="Pb208 14 MeV leakage spectrum from both codes with their ratio below, showing PyNeut low between 4 and 6 MeV"
       width="62%">
</p>

`python run_all.py --plots` writes one of these per case into
`docs/images/cases/`. The panel above is what a POOR spectrum verdict looks
like: integrals and mean escape energy both pass, while the ratio panel shows
where the shape actually departs.

The (n,2n) cases pass on integrals and mean
escape energy but their **leakage spectrum shape** is POOR. The single-collision
(n,xn) emission is verified to match OpenMC node-for-node (`validate_secondary.py`,
χ²/dof ≈ 1) — with the CM→lab boost, unit-base outgoing-energy interpolation and
correlated emission angle all applied — and the Maxwellian fallback fires **0 %**
of the time; the residual is dominated by the emitted neutrons being sampled
independently, with no constraint conserving the total available energy across
the (n,xn) multiplicity.

Also: **analytic benchmarks 13/13 PASS** (Beer–Lambert attenuation; elastic
[α,1] / ξ within 0.3σ), the **cross-section reader matches OpenMC to 0.000 %**
on every channel, and the `pytest` suite has **165 passing** tests.

### Material mixtures vs OpenMC

`Validation/OpenMC_Comparison/validate_mixture.py` validates the per-collision
isotope-selection machinery on a B4C-like **B10 + C12** sphere (strong absorber +
moderator), with **identical number densities** set in both codes:

| Metric | OpenMC | PyNeut | z |
|--------|-------:|-------:|--:|
| leakage | 0.69090 ± 0.00280 | 0.69077 ± 0.00157 | 0.0 |
| avg escape energy (eV) | 772 558 | 771 890 | 0.4 |

Leakage — the integral the mixture machinery most directly drives — agrees within
**0.1σ**, the mean escape energy within **0.4σ**, and the leakage spectrum shape
is **GOOD** (χ²/dof = 1.4). Run at the default `--n 40000`.

### Criticality (k-eigenvalue) vs OpenMC

`Validation/OpenMC_Comparison/validate_keff.py` compares PyNeut's fission-source
power iteration (lower-variance **collision estimator**) against OpenMC's
`eigenvalue` mode on bare U235 metal spheres (same density and data, matched
geometry, graded by the same k_eff z-score):

| Sphere R (cm) | PyNeut k_eff | OpenMC k_eff | z | Status |
|--------------:|-------------:|-------------:|--:|:------:|
| 6.0 (subcritical) | 0.74347 ± 0.00220 | 0.7396 ± 0.0012 | 1.5 | OK |
| 8.7 (supercritical) | 1.02999 ± 0.00250 | 1.0242 ± 0.0015 | 2.0 | OK |
| 11.0 (supercritical) | 1.23619 ± 0.00311 | 1.2283 ± 0.0019 | 2.2 | WARN |

<p align="center">
  <img src="docs/images/keff_vs_radius.png" alt="k_eff versus bare U235 sphere radius for PyNeut and OpenMC" width="49%">
  <img src="docs/images/keff_convergence.png" alt="Fission-source convergence of k_eff over successive generations" width="49%">
</p>

At 2 000 neutrons/generation, 120 generations, 20 inactive. The sweep in
`paper_template/make_figures.py` crosses k = 1 at **R ≈ 8.4 cm**
(k = 1.0013 ± 0.0027), which is the critical radius of *pure* U235 metal at
18.74 g/cm³ — Godiva's 8.7407 cm belongs to the 93.7 %-enriched alloy, which
carries U238/U234 this bare model omits. At R = 8.7 cm both codes agree the
sphere is ~3 % supercritical. PyNeut tracks OpenMC across the range within
**z ≤ 2.2**, the largest departure being at the most supercritical point, and
is consistent with PyNeut's **Watt χ** approximation vs OpenMC's tabulated
ENDF fission spectrum.

**ICSBEP fast-metal benchmarks** (`validate_keff.py --benchmarks`) — the three
canonical bare-sphere critical assemblies, each run through *both* PyNeut and
OpenMC and compared to the published k_eff ≡ 1.0000:

| Benchmark | Fissile | PyNeut k_eff | OpenMC k_eff | z (OMC) | z (bench) |
|-----------|:-------:|-------------:|-------------:|:-------:|:---------:|
| Godiva (HEU-MET-FAST-001) | U235 | 1.0014 ± 0.0021 | 1.0000 ± 0.0013 | 0.6 | 0.6 |
| Jezebel-23 (U233-MET-FAST-001) | U233 | 0.9981 ± 0.0020 | 1.0011 ± 0.0013 | 1.3 | 0.8 |
| Jezebel (PU-MET-FAST-001) | Pu239 | 1.0000 ± 0.0022 | 0.9991 ± 0.0012 | 0.4 | 0.0 |

<p align="center">
  <img src="docs/images/icsbep_benchmarks.png"
       alt="PyNeut and OpenMC k_eff for the Godiva, Jezebel-23 and Jezebel benchmarks against the published unity value"
       width="72%">
</p>

Across **three different fissile isotopes**, PyNeut reproduces each recognized
benchmark within **z ≤ 0.8 of the evaluated k_eff** and **z ≤ 1.3 of OpenMC** on
the identical model; OpenMC itself lands within 0.7σ of every evaluated value. (Requires the ENDF/B-VIII U233/U234/U238 and
Pu239/Pu240/Pu241/Ga data in `endfb/neutron/`; benchmarks whose isotopes are
absent are skipped.)

## Physics Models

- **Elastic scattering** — above 10 eV: static-target kinematics + tabulated
  ENDF angular data; below 10 eV: free-gas vector kinematics (target nucleus
  thermal motion sampled → upscatter, equilibrates to the medium temperature).
- **Inelastic (discrete)** — MT 51–90, two-body kinematics, isotropic in CM.
- **Inelastic (continuum / multiplication)** — (n,2n), (n,3n), continuum via
  ENDF Law 61 / Kalbach–Mann correlated sampling, with Law 4 and Maxwellian
  fallbacks. Outgoing energies use unit-base interpolation between incident
  grids, and CM-frame reactions are boosted to the lab frame; the emitted
  neutrons are sampled independently and banked/tracked.
- **Absorption** — MT 102/103/104/105/106/107 summed.
- **Unresolved-resonance self-shielding** — for heavy nuclides carrying URR
  probability tables (OpenMC HDF5 `urr` group), a cross-section band is sampled
  per collision in the unresolved range and drives both the flight and the
  reaction selection; mean-unbiased outside the self-shielding effect.
- **Fission** — MT 18; in analog `criticality` mode it emits ν (tabulated total
  ν̄) prompt neutrons per fission with a Watt χ spectrum into the next
  generation's fission bank, enabling k_eff. In `shielding` mode the fission
  cross section is still scored as absorption.
- **Cross sections** — ENDF/B-VIII, lin–lin interpolation (verified vs OpenMC).

## Limitations

- **Leakage spectrum shape in multiplying/strongly-moderating cases** — integrals
  and mean escape energies agree, but the leakage *spectrum shape* differs in the
  (n,2n) sphere cases and the B10 epithermal case (χ²/dof ≈ 15–46). The
  single-collision (n,xn) emission is **verified to match OpenMC**
  (`validate_secondary.py`: outgoing-energy χ²/dof ≈ 1 at every incident grid
  node) now that the CM→lab frame transform and unit-base outgoing-energy
  interpolation are applied and the correlated emission angle (Law 61 / Kalbach–
  Mann) is used, and it is **not** the Maxwellian fallback (0 % usage). The
  residual is dominated by the **independent sampling of the emitted
  multiplicity** — no constraint conserves the total available energy across the
  (n,xn) neutrons — which distorts the emission spectrum while leaving the
  neutron balance (and hence the integral leakage) correct.
- **Thermal scattering is free-gas only (no S(α,β))** — the 294 K free-gas model
  is validated (with upscatter, and the earlier hot bias closed), but
  molecular/crystalline binding in real moderators (water, graphite) is not
  modelled. The BEAVRS pin-cell comparison below is therefore a matched-physics
  code-to-code test, not a prediction of the physical k∞.
- The neutrons emitted by (n,2n)/(n,3n) are each sampled independently from the
  single-particle emission distribution, with no constraint conserving the total
  available energy across the multiplicity — the main driver of the POOR
  (n,2n) leakage-spectrum shapes.
- Material mixtures are supported (multiple isotopes per region via
  `Region(composition=[Nuclide, …])` or `Material.mixture`); there is no built-in
  natural-abundance database, so isotopes are specified explicitly.
- k-eigenvalue uses a Watt χ spectrum (not the tabulated ENDF fission spectrum)
  and total ν̄ (prompt + delayed lumped, no delayed-neutron time kinetics).
- Charged-particle and photon products are not transported; neutrons only.
- Source states are built by hand (no built-in spectrum/spatial sampler).
- Geometry limited to planes, spheres and axis-aligned cylinders/boxes.

## Testing

```bash
pytest -q                                              # 165 unit tests
cd Validation/OpenMC_Comparison && python run_all.py   # PyNeut vs OpenMC
python Validation/OpenMC_Comparison/validate_keff.py   # k-eigenvalue vs OpenMC
python Validation/analytic_benchmarks.py               # exact analytic checks
```

## Documentation

- `readme.md` — this file
- `docs/index.md` — project overview & validation results
- `docs/api-reference.md` — class & function reference
- `docs/QUICK_START.md` — worked examples
- `Validation/OpenMC_Comparison/README.md` — validation suite guide

## Future Work

Ordered by engineering leverage (see `docs/index.md` for the full rationale):

1. **Seamless animation pipeline** — replace the bare-coordinate JSON export with
   a self-describing run bundle (serialised geometry plus per-vertex
   energy/event/weight tracks tagged by `parent_id`) and a generic
   three-dimensional Manim scene, so geometry changes no longer require
   hand-editing the animation script; move the export to a binary container as
   particle counts grow.
2. **Criticality refinements** — k-eigenvalue is implemented and validated
   (`criticality.py`: ν-weighted fission bank, Watt χ, generation power iteration,
   lower-variance collision estimator) against OpenMC and the **Godiva /
   Jezebel-23 / Jezebel ICSBEP benchmarks** (z ≤ 1.4 of the evaluated k_eff across
   three fissile isotopes). Remaining work: read the tabulated ENDF χ instead of
   the Watt approximation, separate prompt and delayed ν (time kinetics / β_eff),
   add reflected and thermal-solution benchmarks, and ship a small fissile-isotope
   data subset so the benchmarks are reproducible without a full library download.
3. **Performance** — the profile-guided pass is done (unionized per-nuclide
   energy grids, scalar-float geometry hot path, direct SVT free-gas sampling:
   ~5–6× faster per history, all exact/verified). Remaining headroom is the
   cross-nuclide energy-grid union and Numba/`njit` or an event-vectorised
   rewrite, both of which demand significant array-flattening work.
4. **Physics breadth** — S(α,β) thermal data, multi-temperature and
   Doppler-broadened libraries, a track-length flux estimator, built-in source
   samplers (Watt, volumetric), a natural-abundance database on top of the
   existing material-mixture support, and additional validated isotopes.
5. **Software hygiene** — replace the ten-element positional argument tuple to
   `simulate_single_particle` with a dataclass, and consolidate the duplicated
   isotope `element_map`.

## References

- ENDF/B-VIII Nuclear Data: https://www.nndc.bnl.gov/endf/
- OpenMC Documentation: https://docs.openmc.org/en/stable/methods/index.html
- Duderstadt & Hamilton, *Nuclear Reactor Analysis*
- J. Wood, *Computational Methods in Reactor Shielding*

## License

[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Contact

Youssef — Nuclear and Radiation Engineering, Alexandria University

---

**Disclaimer**: For educational / research-prototyping use only. Not validated
for production, regulatory or safety-critical applications.
