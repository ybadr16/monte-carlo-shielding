---
title: API Reference
layout: default
---

# API Reference

Complete reference documentation for all classes and functions in the Monte Carlo neutron transport code.

---

## Module: `cross_section_read`

### Class: `CrossSectionReader`

Reads and processes nuclear cross-section data from ENDF/B-VIII HDF5 files.

#### Constructor

```python
CrossSectionReader(base_path: str)
```

**Parameters:**
- `base_path` (str): Base directory containing ENDF/B-VIII HDF5 files in `neutron/` subdirectory

**Example:**
```python
reader = CrossSectionReader("./endfb")
```

#### Methods

##### `get_cross_section(element, mt, energy)`

Get microscopic cross-section for a specific reaction and energy.

**Parameters:**
- `element` (str): Nuclide identifier (e.g., "Pb208", "U235")
- `mt` (int): ENDF MT reaction number (1-999)
  - 2: Elastic scattering
  - 18: Fission
  - 102: Radiative capture
- `energy` (float): Neutron energy in eV

**Returns:**
- `float`: Microscopic cross-section in barns (10⁻²⁴ cm²)

**Behavior:** Graceful by design — returns `0.0` (rather than raising) when the
data file is missing, the reaction/MT is absent, or the energy is below the
reaction threshold. This lets the transport kernel sum channels without guarding
every call.

**Example:**
```python
sigma = reader.get_cross_section("Pb208", 2, 1e6)  # 1 MeV elastic scatter
sigma = reader.get_cross_section("Pb208", 9999, 1e6)  # -> 0.0 (no such MT)
```

##### `calculate_macroscopic_xs(microscopic_xs, number_density)`

Convert microscopic to macroscopic cross-section.

**Parameters:**
- `microscopic_xs` (float): Microscopic cross-section in barns
- `number_density` (float): Number density in atoms/cm³

**Returns:**
- `float`: Macroscopic cross-section in cm⁻¹

**Formula:**
```
Σ = σ × N × 10⁻²⁴
```

**Behavior:** Returns `0.0` for a non-positive microscopic input (treated as no
interaction).

**Example:**
```python
N = 3.3e22  # atoms/cm³ for lead
Sigma = reader.calculate_macroscopic_xs(10.0, N)  # cm⁻¹
```

##### `get_cross_sections(element, energy, sampler, N, A)`

Get all macroscopic transport cross-sections at once. Below 10 eV the elastic
lookup energy is shifted to the CM frame to approximate thermal motion.

**Parameters:**
- `element` (str): Nuclide identifier
- `energy` (float): Neutron energy in eV
- `sampler` (VelocitySampler): For thermal-motion energy shift
- `N` (float): Number density in atoms/cm³
- `A` (float): Atomic weight ratio

**Returns (5-tuple, all macroscopic cm⁻¹):**
- `Sigma_el` — elastic scattering
- `Sigma_in` — neutron-emitting inelastic (MT 16, 17, 51–91)
- `Sigma_cap` — absorption (MT 102/103/104/105/106/107)
- `Sigma_fis` — fission (MT 18), 0 if not fissile
- `Sigma_total` — sum of the above

**Example:**
```python
Sigma_el, Sigma_in, Sigma_cap, Sigma_fis, Sigma_t = reader.get_cross_sections(
    "U235", 1e6, sampler, 4.8e22, 233.0
)
```

##### Other reader methods

- `get_inelastic_components(element, energy, N)` → `(total, [(mt, Σ, Q), …])`
  for the available neutron-emitting inelastic channels.
- `get_elastic_mu(element, energy, rng)` → CM scattering cosine sampled from the
  tabulated MT=2 angular distribution.
- `get_secondary_correlated_sample(element, mt, energy, rng)` →
  `(E_out, mu, success)` from ENDF Law 61 correlated data.
- `get_secondary_energy(element, mt, energy, rng)` → secondary energy from
  uncorrelated (Law 4) data, or `None`.

---

## Module: `material`

### Class: `Material`

Represents a material with its physical properties.

#### Constructor

```python
Material(name, density, atomic_mass, atomic_weight_ratio=None)
```

**Parameters:**
- `name` (str): Material name (e.g., "Lead", "Water")
- `density` (float): Density in g/cm³
- `atomic_mass` (float): Atomic mass in g/mol
- `atomic_weight_ratio` (float, optional): Mass ratio A (target/neutron)

**Attributes:**
- `name`: Material name
- `density`: Density (g/cm³)
- `atomic_mass`: Atomic mass (g/mol)
- `atomic_weight_ratio`: A value
- `number_density`: Calculated number density (atoms/cm³)
- `kg_mass`: Single atom mass (kg)

**Example:**
```python
lead = Material(
    name="Lead",
    density=11.35,
    atomic_mass=208,
    atomic_weight_ratio=2.5
)
print(f"Number density: {lead.number_density:.2e} atoms/cm³")
```

#### Methods

##### `calculate_number_density()`

Calculate number density from material properties.

**Returns:**
- `float`: Number density in atoms/cm³

**Formula:**
```
N = (ρ × N_A) / A
```
where ρ = density (g/cm³), N_A = Avogadro's number, A = atomic mass (g/mol)

##### `calculate_atomic_mass_kg()`

Calculate single atom mass in kilograms.

**Returns:**
- `float`: Atomic mass in kg

---

## Module: `medium`

### Class: `Region`

Defines a spatial region using surface boundaries.

#### Constructor

```python
Region(surfaces=None, operation="intersection", name=None, 
       priority=0, is_void=False, element=None)
```

**Parameters:**
- `surfaces` (list): List of surface objects defining boundaries
- `operation` (str): Boolean operation - "intersection", "union", "complement", "difference"
- `name` (str): Region name for identification
- `priority` (int): Priority level (higher = higher priority at overlaps)
- `is_void` (bool): If True, no particle interactions occur
- `element` (str): Nuclear data identifier (e.g., "Pb208")

**Example:**
```python
cylinder_region = Region(
    surfaces=[
        Cylinder("z", 10, (0, 0, 0)),
        Plane(0, 0, -1, 10),
        Plane(0, 0, 1, 10)
    ],
    operation="intersection",
    name="Lead Shield",
    priority=1,
    element="Pb208"
)
```

#### Methods

##### `contains(x, y, z)`

Check if a point is inside the region.

**Parameters:**
- `x, y, z` (float): Coordinates in cm

**Returns:**
- `bool`: True if point is inside region

##### `add_surface(surface)`

Add a surface to the region's surface list.

---

### Class: `Plane`

Defines a half-space. **Convention:** `Plane(A, B, C, D)` represents the plane
`A·x + B·y + C·z = D`, and a point is *inside* (``evaluate`` ≤ 0) when
`A·x + B·y + C·z ≤ D`. Internally the normal is normalized and the stored
offset is `-D/‖n‖`, so to describe the plane `z = 5` you write `Plane(0, 0, 1, 5)`.

#### Constructor

```python
Plane(A, B, C, D)
```

**Parameters:**
- `A, B, C` (float): Normal vector components (automatically normalized)
- `D` (float): Plane offset (in the `A·x+B·y+C·z = D` sense)

**Example:**
```python
plane = Plane(0, 0, 1, 5)    # the plane z = 5  (inside: z <= 5)
plane = Plane(-1, 0, 0, 0)   # the plane x = 0  (inside: x >= 0)
```

#### Methods

##### `evaluate(x, y, z)`

Evaluate plane equation at point. Returns ≤ 0 if point satisfies inequality.

##### `nearest_surface_method(x, y, z, u, v, w)`

Calculate distance to plane intersection along ray.

**Parameters:**
- `x, y, z`: Starting point
- `u, v, w`: Direction vector (not necessarily normalized)

**Returns:**
- `float`: Distance to intersection, or None if parallel/no intersection

##### `normal(x, y, z)`

Get outward normal vector at point (normalized).

**Returns:**
- `tuple`: (nx, ny, nz)

---

### Class: `Cylinder`

Infinite cylinder along coordinate axis.

#### Constructor

```python
Cylinder(axis, radius, center)
```

**Parameters:**
- `axis` (str): "x", "y", or "z" - cylinder axis direction
- `radius` (float): Cylinder radius in cm
- `center` (tuple): (x0, y0, z0) - center point coordinates

**Example:**
```python
# Cylinder along z-axis, radius 5 cm, centered at origin
cyl = Cylinder("z", 5.0, (0, 0, 0))
```

#### Methods

##### `evaluate(x, y, z)`

Evaluate cylinder equation. Returns ≤ 0 if inside.

##### `nearest_surface_method(x, y, z, u, v, w)`

Calculate distance to cylinder surface along ray.

**Returns:**
- `float`: Distance to intersection
- Special cases:
  - Returns None if ray parallel to axis with no intersection
  - Returns exit distance if inside
  - Returns entry distance if outside

##### `normal(x, y, z)`

Get radial normal vector at point on surface.

---

### Class: `Sphere`

Spherical surface.

#### Constructor

```python
Sphere(center, radius)
```

**Parameters:**
- `center` (tuple): (x0, y0, z0)
- `radius` (float): Sphere radius in cm

**Example:**
```python
sphere = Sphere((0, 0, 0), 10.0)
```

---

### Class: `Box`

Rectangular box region (convenience class, inherits from Region).

#### Constructor

```python
Box(x_min, x_max, y_min, y_max, z_min, z_max)
```

**Parameters:**
- `x_min, x_max`: X bounds
- `y_min, y_max`: Y bounds
- `z_min, z_max`: Z bounds

**Example:**
```python
# 10×20×30 cm box centered at origin
box = Box(-5, 5, -10, 10, -15, 15)
```

---

## Module: `physics`

### Function: `elastic_scattering`

Simulate elastic scattering event.

```python
elastic_scattering(initial_energy, A, sampler, rng)
```

**Parameters:**
- `initial_energy` (float): Incident neutron energy (eV)
- `A` (float): Target atomic weight ratio
- `sampler` (VelocitySampler): For thermal motion
- `rng` (RNGHandler): Random number generator

**Returns:**
- `E_prime` (float): Scattered neutron energy (eV)
- `mu_cm` (float): Scattering cosine in CM frame
- `mu_lab` (float): Scattering cosine in lab frame

**Physics:**
Accounts for:
- Target thermal motion (free-gas sampling for E < 10 eV)
- Center-of-mass ↔ lab frame vector transformation
- Isotropic scattering in the CM frame (the lab-frame forward bias it produces
  is the genuine kinematic 2/(3A) effect, not an imposed anisotropy)

> The transport loop uses this free-gas kernel **below 10 eV** (where target
> motion and upscatter matter); **above 10 eV** it uses static-target kinematics
> with the **tabulated ENDF angular distribution** (`get_elastic_mu`). The free
> gas requires the target speed *and* its cosine relative to the neutron from the
> same accepted sample — `VelocitySampler.sample_velocity(..., return_mu=True)`.

---

### Function: `inelastic_scattering(initial_energy, A, Q_value, sampler, rng)`

Discrete-level inelastic scattering (two-body kinematics, isotropic in CM).
Falls back to elastic if the available CM energy is below `Q_value`.

**Returns:** `(E_prime, mu_cm, mu_lab)`.

### Function: `sample_maxwellian(temperature_ev, rng)`

Sample an energy from a Maxwellian (evaporation/fission) spectrum
`f(E) ∝ √E·exp(-E/T)`. Used as the secondary-energy fallback for (n,xn)/continuum.

### Function: `get_nuclear_temperature(energy, A)`

Fermi-gas nuclear temperature `T = √(U/a)`, `a = A/8 MeV⁻¹`, for the evaporation
model. Returns eV.

---

### Function: `calculate_E_cm_prime`

Calculate neutron energy in center-of-mass frame before collision.

```python
calculate_E_cm_prime(initial_energy, A, sampler)
```

**Parameters:**
- `initial_energy` (float): Lab energy (eV)
- `A` (float): Mass ratio
- `sampler` (VelocitySampler): For velocity sampling

**Returns:**
- `float`: Energy in CM frame (eV)

---

### Function: `calculate_E_prime`

Calculate post-collision energy in lab frame.

```python
calculate_E_prime(E_cm_prime, initial_energy, A, rng)
```

**Parameters:**
- `E_cm_prime` (float): CM energy before collision
- `initial_energy` (float): Lab energy before collision
- `A` (float): Mass ratio
- `rng` (RNGHandler): RNG for sampling scattering angle

**Returns:**
- `E_prime` (float): Lab energy after collision
- `mu_cm` (float): Scattering cosine in CM

---

### Function: `sample_new_direction_cosines`

Sample new flight direction after scattering.

```python
sample_new_direction_cosines(u, v, w, mu_lab, rng)
```

**Parameters:**
- `u, v, w` (float): Incident direction cosines
- `mu_lab` (float): Scattering cosine in lab frame
- `rng` (RNGHandler): Random number generator

**Returns:**
- `u_new, v_new, w_new` (float): New direction cosines (normalized)
- `phi` (float): Azimuthal scattering angle

**Note:** Assumes azimuthal symmetry (isotropic in φ).

---

## Module: `vt_calc`

### Class: `VelocitySampler`

Samples target nucleus velocities from Maxwell-Boltzmann distribution.

#### Constructor

```python
VelocitySampler(mass, temperature=294)
```

**Parameters:**
- `mass` (float): Target nucleus mass (kg)
- `temperature` (float): Temperature (K), default 294K

**Example:**
```python
m_Pb = 208 * 1.66e-27  # kg
sampler = VelocitySampler(mass=m_Pb, temperature=300)
```

#### Methods

##### `sample_velocity(vn, max_attempts=1000)`

Sample target velocity magnitude accounting for relative motion.

**Parameters:**
- `vn` (float): Neutron velocity magnitude (m/s)
- `max_attempts` (int): Maximum rejection sampling attempts

**Returns:**
- `float`: Target velocity magnitude (m/s)

**Raises:**
- `ValueError`: If sampling fails after max_attempts

**Algorithm:** Uses acceptance-rejection sampling with weighted proposal distribution.

---

## Module: `geometry`

### Function: `calculate_nearest_boundary`

Find nearest surface intersection along particle trajectory.

```python
calculate_nearest_boundary(state, regions, u, v, w, epsilon=1e-10)
```

**Parameters:**
- `state` (dict): Particle state with x, y, z coordinates
- `regions` (list): List of Region objects
- `u, v, w` (float): Direction vector components
- `epsilon` (float): Tolerance for coincident surfaces

**Returns:**
- `nearest_point` (tuple): (x, y, z) of intersection
- `nearest_region` (Region): Region beyond boundary
- `nearest_distance` (float): Distance to boundary

**Returns:**
- `(None, None, inf)` if no boundary found (particle escapes)

---

### Function: `calculate_direction_cosines`

Calculate direction cosines from two points.

```python
calculate_direction_cosines(x, y, z, x_prev, y_prev, z_prev)
```

**Returns:**
- `u, v, w` (float): Normalized direction cosines

---

### Function: `count_coordinates_in_boundary`

Count coordinates within specified box.

```python
count_coordinates_in_boundary(coordinates, x_bounds, y_bounds, z_bounds)
```

**Parameters:**
- `coordinates` (list): List of (x, y, z) tuples
- `x_bounds, y_bounds, z_bounds` (tuple): (min, max) for each dimension

**Returns:**
- `int`: Number of coordinates in box

---

## Module: `simulation`

### Function: `simulate_single_particle(args)`

Entry point mapped over a multiprocessing `Pool`. Simulates one source
particle **and all of its secondary neutrons** (a bank/stack handles (n,2n),
(n,3n) descendants).

**`args` is a 10-tuple, in order:**

```python
(state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings)
```

| # | Item | Notes |
|---|------|-------|
| 1 | `state` (dict) | `x, y, z` (cm); `theta, phi` (rad); `energy` (eV); **`weight`** (float, e.g. 1.0); `has_interacted` (bool) |
| 2 | `reader` | `CrossSectionReader` |
| 3 | `mediums` | list of `Region` |
| 4 | `A` | atomic weight ratio |
| 5 | `N` | number density (atoms/cm³) |
| 6 | `sampler` | `VelocitySampler` |
| 7 | `region_bounds` | optional `(x_min,x_max,y_min,y_max,z_min,z_max)` detector box, or `None` |
| 8 | `track_coordinates` | `bool` — record trajectory |
| 9 | `rng` | `RNGHandler` |
| 10 | `settings` | `Settings` — selects analog vs implicit-capture transport |

> The `weight` field and the trailing `settings` element are **required** — they
> were added when survival biasing and secondary-neutron banking were
> introduced.

**Returns** a result `dict`:

| Key | Meaning |
|-----|---------|
| `result` | `"simulated"` |
| `absorbed_weight` | total weight removed by capture |
| `absorbed_coords` | list of `(x, y, z, weight_lost)` for mesh scoring |
| `fissioned` | number of fission events |
| `final_energy` | weighted-average escape energy (eV) |
| `final_weight` | total weight that leaked out |
| `region_detected` | `bool` — entered the detector box |
| `trajectory` | `[(x,y,z), …]` if tracking, else `None` |

The core loop lives in `run_particle_kernel(...)` (same arguments); it samples
distance to collision, handles void/boundary crossings, applies the analog or
implicit-capture collision kernel, and dispatches elastic/inelastic kinematics.

---

## Module: `settings`

### Class: `Settings`

```python
Settings(mode="shielding", particles=1000, batches=10)
```

| `mode` | `use_implicit_capture` | `weight_cutoff` | `roulette_survival_prob` |
|--------|------------------------|-----------------|--------------------------|
| `"shielding"` | True | 1e-4 | 0.1 |
| `"criticality"` | False (analog) | 0.0 | 1.0 |

**Attributes:** `mode`, `num_particles`, `batches`, `use_implicit_capture`,
`weight_cutoff`, `roulette_survival_prob`.

---

## Module: `mesh`

### Class: `MeshTally`

A regular Cartesian mesh that accumulates scored weight (e.g. energy/dose
deposition) and exports to a legacy VTK file for ParaView.

```python
MeshTally(x_bounds, y_bounds, z_bounds, dims)
```

**Parameters:** `x_bounds`/`y_bounds`/`z_bounds` are `(min, max)` tuples;
`dims` is `(nx, ny, nz)`.

**Methods:**
- `score(x, y, z, weight)` — add `weight` to the voxel containing the point.
- `write_vtk(filename)` — write `STRUCTURED_POINTS` cell data to `filename`.

---

## Module: `tally`

### Class: `Tally`

Accumulates simulation statistics.

#### Constructor

```python
Tally()
```

#### Attributes

- `results` (dict): `{"absorbed": float, "escaped": int, "killed": int}`
  (`absorbed` accumulates weight; `killed` counts roulette deaths)
- `absorbed_coordinates` (list): `(x, y, z, weight)` absorption events
- `energy_spectrum` (list): Per-particle escape energies
- `region_count` (int): Particles detected in the optional detector box

#### Methods

##### `update(result, absorbed=None, final_energy=None, region_detected=False)`

Update tallies with single particle result.

##### `merge_partial_results(partial_results)`

Merge results from multiprocessing worker.

##### `print_summary(num_particles)`

Print formatted summary of simulation results.

##### `get_results()`

Return dictionary of all tally data.

---

## Module: `random_number_generator`

### Class: `RNGHandler`

Wrapper for NumPy random number generation.

#### Constructor

```python
RNGHandler(seed=None)
```

**Parameters:**
- `seed` (int, optional): Random seed for reproducibility

#### Methods

##### `random()`
Returns float in [0, 1)

##### `uniform(low, high)`
Returns float in [low, high)

##### `log_uniform(scale)`
Returns -ln(1-ξ)/scale for exponential sampling

##### `choice(a, size=None, replace=True, p=None)`
Random selection from sequence

---

## Usage Example

Complete workflow:

```python
from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Cylinder, Plane
from src.vt_calc import VelocitySampler
from src.simulation import simulate_single_particle
from src.tally import Tally
from src.random_number_generator import RNGHandler
from src.settings import Settings
from multiprocessing import Pool
import numpy as np

# 1. Initialize cross-section reader
reader = CrossSectionReader("./endfb")

# 2. Define material (single isotope per region)
lead = Material("Lead", 11.35, 208, atomic_weight_ratio=207.2)

# 3. Create velocity sampler
sampler = VelocitySampler(mass=lead.kg_mass)

# 4. Define geometry
regions = [
    Region(
        surfaces=[
            Cylinder("z", 10, (0, 0, 0)),
            Plane(0, 0, -1, 10),   # z >= -10
            Plane(0, 0, 1, 10),    # z <= 10
        ],
        name="Shield",
        priority=1,
        element="Pb208"
    )
]

# 5. Settings + particles (note the `weight` field)
settings = Settings(mode="shielding", particles=100)
rngs = [RNGHandler(seed=12345 + i) for i in range(settings.num_particles)]
particle_states = [
    {
        "x": -15.0, "y": 0.0, "z": 0.0,
        "theta": 0.0, "phi": 0.0,
        "energy": 1e6,
        "weight": 1.0,
        "has_interacted": False,
    }
    for _ in rngs
]

# 6. Prepare simulation arguments (10-tuple, trailing `settings`)
args = [
    (state, reader, regions, lead.atomic_weight_ratio,
     lead.number_density, sampler, None, False, rng, settings)
    for state, rng in zip(particle_states, rngs)
]

# 7. Run simulation
with Pool() as pool:
    results = pool.map(simulate_single_particle, args)

# 8. Collect results
tally = Tally()
for result in results:
    tally.merge_partial_results(result)

tally.print_summary(settings.num_particles)
```