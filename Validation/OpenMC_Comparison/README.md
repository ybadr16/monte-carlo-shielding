# PyNeut vs OpenMC validation suite

A fully scripted comparison of PyNeut against
[OpenMC](https://docs.openmc.org) on a set of single-isotope shielding
problems. Every OpenMC reference is generated **live** from the same ENDF/B
HDF5 data PyNeut reads, with identical geometry and source, so the two codes
are directly comparable and the numbers are reproducible on any machine with
`openmc` installed.

## Layout

| File | Purpose |
|------|---------|
| `_common.py` | Shared infrastructure: builds `endfb/cross_sections.xml` on the fly; runs PyNeut (`run_pyneut`) and OpenMC live (`run_openmc`); computes the spectrum χ²/dof (`compare_spectra`) and fallback summary; formats results (`validate_case`, `make_result`, `print_result`). |
| `validate_<iso>.py` | Per-isotope case definitions (Pb208, Fe56, Be9, Al27, C12, O16). Each is runnable standalone. |
| `validate_thermal.py` | Thermal / epithermal regime (B10, Cd112, C12) — probes the part of phase space PyNeut is *least* validated in. |
| `validate_xs.py` | Cross-section interpolation check — PyNeut's reader vs `openmc.data.IncidentNeutron.from_hdf5`, channel by channel. No `cross_sections.xml` needed. |
| `run_all.py` | Master runner. Executes every isotope module + the XS check, prints a summary table, and writes `validation_results.csv`. |
| `main_benchmark.py` | Minimal standalone PyNeut-only Pb208 sphere benchmark. |
| `../analytic_benchmarks.py` | PyNeut vs **exact analytic** results (no OpenMC): Beer–Lambert attenuation and elastic [α,1] / ξ kinematics. |

## Running

```bash
cd Validation/OpenMC_Comparison

# Full suite (N = 10 000 particles/case)
python run_all.py

# Faster smoke test, subset of isotopes
python run_all.py --n 2000 --isotopes Pb208,Fe56

# Cross-section reader check only
python validate_xs.py

# Exact analytic benchmarks (no OpenMC needed)
python ../analytic_benchmarks.py
```

Requirements: `openmc`, `numpy`, `h5py` and the ENDF/B HDF5 files in
`endfb/neutron/`. The OpenMC data library (`endfb/cross_sections.xml`) is built
automatically from those files the first time a live run is requested.

## What is compared

1. **Integral metrics** — leakage fraction and average escape energy, each with a
   standard error (batch means for PyNeut; OpenMC's native batch stats, with
   ratio error propagation for the energy). The average energy uses the **same
   energy histogram** in both codes (`E_NBINS` bins, `energy_grid()`), so binning
   bias cancels.
2. **Spectrum shape** — the full leakage spectrum, via a coarse-group χ²/dof
   goodness-of-fit (`compare_spectra`), reported *separately* from the integral
   status (GOOD ≤ 2, FAIR ≤ 5, POOR > 5).
3. **Fallback diagnostics** — the transport kernel counts how each (n,xn)
   secondary energy was sampled (Law 61 / Law 4 / Maxwellian evaporation); the
   per-case "Max%" column is the fraction that hit the evaporation fallback.

## Grading (integral status)

```
z = |OpenMC - PyNeut| / sqrt(sem_OMC^2 + sem_PyN^2)
```

| Status | Condition (worse of leakage / energy) |
|--------|----------------------------------------|
| OK | z ≤ 2 (consistent within 2σ) |
| WARNING | 2 < z ≤ 4 |
| FAIL | z > 4 (statistically significant bias) |
| NO_SIGMA | uncertainty unavailable (e.g. a single batch) |

See `../../docs/index.md` for the latest committed results table and the
interpretation by energy regime.
