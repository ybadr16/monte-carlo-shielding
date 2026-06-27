# Quick Start Guide

Step-by-step examples for common use cases.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/ybadr16/PyNeut
cd PyNeut
```

2. Install dependencies:
```bash
pip install -e .          # or: pip install numpy h5py
```

3. ENDF/B-VIII data:
- HDF5 files live in `./endfb/neutron/`. The repo ships Al27, B10, Be9, C12,
  Cd112, Fe56, O16, Pb207, Pb208, U235.

> ### ⚠️ API note (read this first)
> Every example below passes a **10-element** args tuple to
> `simulate_single_particle`, ending in a `Settings` object, and every particle
> `state` includes a **`weight`** field:
>
> ```python
> from src.settings import Settings
> settings = Settings(mode="shielding", particles=num_particles)
>
> state = {"x":…, "y":…, "z":…, "theta":…, "phi":…,
>          "energy": 1e6, "weight": 1.0, "has_interacted": False}
>
> args = (state, reader, regions, A, N, sampler,
>         region_bounds, track_coordinates, rng, settings)   # <- trailing settings
> ```
>
> `mode="shielding"` enables implicit capture (weighted leakage); `mode="criticality"`
> runs analog transport. Energies are in **eV**, distances in **cm**.

## Example 1: Simple Slab Shielding

Simulate neutron transmission through a lead slab.

```python
from src.cross_section_read import CrossSectionReader
from src.material import Material
from src.medium import Region, Plane
from src.vt_calc import VelocitySampler
from src.simulation import simulate_single_particle
from src.tally import Tally
from src.random_number_generator import RNGHandler
from src.settings import Settings
from multiprocessing import Pool
import numpy as np

# Setup
reader = CrossSectionReader("./endfb")
lead = Material("Lead", 11.35, 208, atomic_weight_ratio=207.2)
sampler = VelocitySampler(lead.kg_mass)
settings = Settings(mode="shielding", particles=1000)

# Geometry: 10 cm lead slab from x=0 to x=10
regions = [
    Region(
        surfaces=[
            Plane(-1, 0, 0, 0),    # x >= 0
            Plane(1, 0, 0, 10),     # x <= 10
            Plane(0, -1, 0, 50),    # y >= -50
            Plane(0, 1, 0, 50),     # y <= 50
            Plane(0, 0, -1, 50),    # z >= -50
            Plane(0, 0, 1, 50)      # z <= 50
        ],
        name="Lead Slab",
        priority=1,
        element="Pb208"
    ),
    Region(
        surfaces=[
            Plane(-1, 0, 0, 100),
            Plane(1, 0, 0, 100),
            Plane(0, -1, 0, 100),
            Plane(0, 1, 0, 100),
            Plane(0, 0, -1, 100),
            Plane(0, 0, 1, 100)
        ],
        name="Void",
        priority=0,
        is_void=True
    )
]

# Source: 1000 neutrons at 1 MeV, starting at x=-5
num_particles = settings.num_particles
rngs = [RNGHandler(seed=12345 + i) for i in range(num_particles)]
particle_states = [
    {
        "x": -5, "y": 0, "z": 0,
        "theta": 0, "phi": 0,  # Moving in +x direction
        "energy": 1e6,
        "weight": 1.0,
        "has_interacted": False
    }
    for _ in rngs
]

# Simulate (note trailing `settings`)
args = [
    (state, reader, regions, lead.atomic_weight_ratio,
     lead.number_density, sampler, None, False, rng, settings)
    for state, rng in zip(particle_states, rngs)
]

with Pool() as pool:
    results = pool.map(simulate_single_particle, args)

# Analyze
tally = Tally()
for result in results:
    tally.merge_partial_results(result)

tally.print_summary(num_particles)

# Calculate transmission
transmission = tally.results["escaped"] / num_particles
print(f"Transmission coefficient: {transmission:.3f}")
```

## Example 2: Cylindrical Detector Response

Simulate neutron flux in a cylindrical detector.

```python
from src.medium import Cylinder

# Geometry: Small detector cylinder at x=15
detector_region = (14.9, 15.1, -5, 5, -5, 5)  # x, y, z bounds

regions = [
    # Lead shield
    Region(
        surfaces=[
            Cylinder("z", 10, (0, 0, 0)),
            Plane(0, 0, -1, 10),
            Plane(0, 0, 1, 10)
        ],
        name="Shield",
        priority=1,
        element="Pb208"
    ),
    # Void
    Region(
        surfaces=[
            Plane(-1, 0, 0, 20),
            Plane(1, 0, 0, 20),
            Plane(0, -1, 0, 20),
            Plane(0, 1, 0, 20),
            Plane(0, 0, -1, 20),
            Plane(0, 0, 1, 20)
        ],
        name="Void",
        priority=0,
        is_void=True
    )
]

# Source: Isotropic point source
particle_states = [
    {
        "x": -15, "y": 0, "z": 0,
        "theta": rng.uniform(0, np.pi),
        "phi": rng.uniform(0, 2*np.pi),
        "energy": 1e6,
        "weight": 1.0,
        "has_interacted": False
    }
    for rng in rngs
]

# Run simulation with detector region (note trailing `settings`)
args = [
    (state, reader, regions, lead.atomic_weight_ratio,
     lead.number_density, sampler, detector_region, False, rng, settings)
    for state, rng in zip(particle_states, rngs)
]

with Pool() as pool:
    results = pool.map(simulate_single_particle, args)

# Count detections
tally = Tally()
for result in results:
    tally.merge_partial_results(result)

detection_efficiency = tally.region_count / num_particles
print(f"Detection efficiency: {detection_efficiency:.4f}")
```

## Example 3: Energy Spectrum Analysis

Track energy degradation through multiple scatterings.

```python
import matplotlib.pyplot as plt

# Run simulation (same as Example 1)
# ...

# Extract energy spectrum
final_energies = tally.energy_spectrum

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(final_energies, bins=50, edgecolor='black')
plt.xlabel('Final Energy (eV)')
plt.ylabel('Count')
plt.title('Neutron Energy Spectrum After Shield')
plt.yscale('log')
plt.xscale('log')
plt.grid(True, alpha=0.3)
plt.savefig('energy_spectrum.png', dpi=300)
plt.show()

# Statistics
print(f"Mean final energy: {np.mean(final_energies):.2e} eV")
print(f"Median final energy: {np.median(final_energies):.2e} eV")
print(f"Min final energy: {np.min(final_energies):.2e} eV")
print(f"Max final energy: {np.max(final_energies):.2e} eV")
```

## Example 4: Trajectory Visualization

Record and visualize particle paths.

```python
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Enable trajectory tracking
track_coordinates = True

args = [
    (state, reader, regions, lead.atomic_weight_ratio,
     lead.number_density, sampler, None, track_coordinates, rng, settings)
    for state, rng in zip(particle_states[:10], rngs[:10])  # Only 10 particles
]

# Run simulation
with Pool() as pool:
    results = pool.map(simulate_single_particle, args)

# Save trajectories
trajectories = {}
for idx, result in enumerate(results):
    if result["trajectory"]:
        trajectories[idx] = result["trajectory"]

with open("trajectories.json", "w") as f:
    json.dump(trajectories, f, indent=2)

# Plot 3D trajectories
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

for idx, traj in trajectories.items():
    coords = np.array(traj)
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
            alpha=0.6, linewidth=1, label=f"Particle {idx}")

ax.set_xlabel('X (cm)')
ax.set_ylabel('Y (cm)')
ax.set_zlabel('Z (cm)')
ax.set_title('Neutron Trajectories')
ax.legend(fontsize=8)
plt.savefig('trajectories_3d.png', dpi=300)
plt.show()
```

## Example 5: Multi-Energy Source

Simulate with energy distribution.

```python
# Energy distribution (e.g., fission spectrum)
def watt_spectrum(rng):
    """Sample from Watt fission spectrum"""
    a = 0.988  # MeV
    b = 2.249  # MeV^-1
    
    # Simplified sampling (not exact)
    while True:
        E = -a * np.log(rng.random())
        if rng.random() < np.exp(-E / a) * np.sinh(np.sqrt(b * E)):
            return E * 1e6  # Convert to eV

# Create particles with energy distribution
particle_states = [
    {
        "x": -15, "y": 0, "z": 0,
        "theta": rng.uniform(0, np.pi),
        "phi": rng.uniform(0, 2*np.pi),
        "energy": watt_spectrum(rng),
        "weight": 1.0,
        "has_interacted": False
    }
    for rng in rngs
]

# Run simulation
# ... (same as before)

# Analyze initial vs final energies
initial_energies = [state["energy"] for state in particle_states]

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(initial_energies, bins=50, alpha=0.7, label='Initial')
plt.xlabel('Energy (eV)')
plt.ylabel('Count')
plt.title('Initial Energy Distribution')
plt.xscale('log')
plt.legend()

plt.subplot(1, 2, 2)
plt.hist(tally.energy_spectrum, bins=50, alpha=0.7, 
         label='Final', color='orange')
plt.xlabel('Energy (eV)')
plt.ylabel('Count')
plt.title('Final Energy Distribution')
plt.xscale('log')
plt.legend()

plt.tight_layout()
plt.savefig('energy_comparison.png', dpi=300)
plt.show()
```

## Example 6: Parametric Study

Study transmission vs shield thickness.

```python
thicknesses = [1, 2, 5, 10, 15, 20]  # cm
transmissions = []

for thickness in thicknesses:
    print(f"Simulating {thickness} cm shield...")
    
    # Update geometry
    regions = [
        Region(
            surfaces=[
                Plane(-1, 0, 0, 0),
                Plane(1, 0, 0, thickness),
                Plane(0, -1, 0, 50),
                Plane(0, 1, 0, 50),
                Plane(0, 0, -1, 50),
                Plane(0, 0, 1, 50)
            ],
            name="Lead Slab",
            priority=1,
            element="Pb208"
        ),
        Region(
            surfaces=[
                Plane(-1, 0, 0, 100),
                Plane(1, 0, 0, 100),
                Plane(0, -1, 0, 100),
                Plane(0, 1, 0, 100),
                Plane(0, 0, -1, 100),
                Plane(0, 0, 1, 100)
            ],
            name="Void",
            priority=0,
            is_void=True
        )
    ]
    
    # Run simulation (note trailing `settings`)
    args = [
        (state, reader, regions, lead.atomic_weight_ratio,
         lead.number_density, sampler, None, False, rng, settings)
        for state, rng in zip(particle_states, rngs)
    ]
    
    with Pool() as pool:
        results = pool.map(simulate_single_particle, args)
    
    # Calculate transmission
    tally = Tally()
    for result in results:
        tally.merge_partial_results(result)
    
    transmission = tally.results["escaped"] / num_particles
    transmissions.append(transmission)
    print(f"  Transmission: {transmission:.4f}")

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(thicknesses, transmissions, 'o-', linewidth=2, markersize=8)
plt.xlabel('Shield Thickness (cm)')
plt.ylabel('Transmission Coefficient')
plt.title('Neutron Transmission vs Lead Shield Thickness')
plt.grid(True, alpha=0.3)
plt.yscale('log')
plt.savefig('transmission_vs_thickness.png', dpi=300)
plt.show()

# Fit exponential attenuation
from scipy.optimize import curve_fit

def exponential(x, a, b):
    return a * np.exp(-b * x)

popt, _ = curve_fit(exponential, thicknesses, transmissions)
print(f"Fitted attenuation coefficient: {popt[1]:.4f} cm^-1")
print(f"Mean free path: {1/popt[1]:.2f} cm")
```

## Example 7: Complex Geometry

Multi-region geometry with different materials.

```python
from src.medium import Sphere

# Three-layer spherical geometry
regions = [
    # Inner sphere (void - represents source region)
    Region(
        surfaces=[Sphere((0, 0, 0), 2)],
        name="Source",
        priority=3,
        is_void=True
    ),
    # First layer (lead)
    Region(
        surfaces=[
            Sphere((0, 0, 0), 7),
            Region(
                surfaces=[Sphere((0, 0, 0), 2)],
                operation="complement"
            )
        ],
        operation="intersection",
        name="Lead Layer",
        priority=2,
        element="Pb208"
    ),
    # Outer void
    Region(
        surfaces=[
            Plane(-1, 0, 0, 50),
            Plane(1, 0, 0, 50),
            Plane(0, -1, 0, 50),
            Plane(0, 1, 0, 50),
            Plane(0, 0, -1, 50),
            Plane(0, 0, 1, 50)
        ],
        name="Outer Void",
        priority=0,
        is_void=True
    )
]

# Isotropic point source at origin
particle_states = [
    {
        "x": 0, "y": 0, "z": 0,
        "theta": np.arccos(2*rng.random() - 1),  # Isotropic
        "phi": 2*np.pi*rng.random(),
        "energy": 1e6,
        "weight": 1.0,
        "has_interacted": False
    }
    for rng in rngs
]

# Run simulation
# ...
```

## Common Pitfalls

### 1. Particles Starting Outside Geometry
```python
# BAD: Particle outside all regions
state = {"x": 1000, "y": 0, "z": 0, ...}

# GOOD: Verify particle is in a region
def validate_position(x, y, z, regions):
    for region in regions:
        if region.contains(x, y, z):
            return True
    return False

if not validate_position(state["x"], state["y"], state["z"], regions):
    print("Warning: Particle outside geometry!")
```

### 2. Inconsistent Units
```python
# BAD: Mixing units
energy = 1  # MeV - WRONG!

# GOOD: Always use eV for energy, cm for distance
energy = 1e6  # eV
position = 10  # cm
```

### 3. Overlapping Regions Without Priority
```python
# BAD: Two regions overlap with same priority
region1 = Region(..., priority=1)
region2 = Region(..., priority=1)  # Ambiguous if they overlap

# GOOD: Use different priorities
region1 = Region(..., priority=2)  # Higher priority
region2 = Region(..., priority=1)
```

### 4. Forgetting Void Region
```python
# BAD: No void region - particles have nowhere to escape
regions = [shield_region]

# GOOD: Always include void region
regions = [shield_region, void_region]
```

## Performance Tips

1. **Use multiprocessing**: Scales linearly with CPU cores
2. **Disable trajectory tracking** for production runs
3. **Reduce number of particles** during development
4. **Use larger epsilon** for faster (but less accurate) simulations
5. **Cache cross-section data** if running many similar simulations

## Next Steps

- See [`api-reference.md`](api-reference.md) for detailed function documentation
- See [`index.md`](index.md) for the validation results and current limitations
- See `main.py` for a complete runnable example with a mesh tally
- See `Validation/OpenMC_Comparison/` for the PyNeut-vs-OpenMC suite
