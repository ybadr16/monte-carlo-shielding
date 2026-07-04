import numpy as np
from src.simulation import simulate_single_particle
from src.medium import Region, Plane
from src.settings import Settings
from src.random_number_generator import RNGHandler

# --- 1. MOCK READER ---
class MockReader:
    def get_cross_sections(self, element, energy, sampler, N):
        # PURE SCATTERER: Sigma_s = 1.0, Sigma_a = 0.0
        # Sigma_t = 1.0
        return 1.0, 0.0, 0.0, 1.0

# --- 2. MOCK SAMPLER (THE FIX) ---
class MockSampler:
    """Simulates a target at 0 Kelvin (Rest)."""
    def sample_velocity(self, vn=None, rng=None, return_mu=False):
        return (0.0, 0.0) if return_mu else 0.0

def run_reflection_test():
    print("Running Reflection (Albedo) Validation...")
    print("firing neutrons into a thick wall of pure scatterer.")

    # Very thick slab (20 mean free paths)
    T = 20.0
    num_particles = 2000

    # Use Criticality (Analog) mode
    settings = Settings(mode="criticality")

    # Geometry
    slab = Region(
        surfaces=[
            Plane(0,0,-1, 0),   # Front (z>=0)
            Plane(0,0,1, T),    # Back (z<=T)
            Plane(-1,0,0, 100), # Wide boundaries
            Plane(1,0,0, 100),
            Plane(0,-1,0, 100),
            Plane(0,1,0, 100)
        ],
        name="Scatterium", priority=1, element="Scatterium"
    )

    mediums = [slab]

    # Counters
    reflected = 0
    transmitted = 0
    absorbed = 0

    print("-" * 60)
    print(f"{'Metric':<20} | {'Expected':<15} | {'Result':<15}")
    print("-" * 60)

    for i in range(num_particles):
        rng = RNGHandler(seed=i)

        # Start just inside front face
        state = {
            "x": 0.0, "y": 0.0, "z": 0.0001,
            "theta": 0.0, "phi": 0.0,
            "has_interacted": False, "energy": 1e6, "weight": 1.0
        }

        track = True

        # Pass the MockSampler instead of None
        reader = MockReader()
        sampler = MockSampler()

        args = (state, reader, mediums, 1.0, 1.0, sampler, None, track, rng, settings)
        data = simulate_single_particle(args)

        if data["result"] == "escaped":
            # Check exit location
            # Note: trajectory[-1] is the last recorded point.
            # In your simulation loop, you append coordinate at step start.
            # To be precise, we check if the LAST recorded Z is near 0 or near T.
            last_x, last_y, last_z = data["trajectory"][-1]

            # Simple check: Which boundary is closer?
            if abs(last_z - 0) < abs(last_z - T):
                reflected += 1
            else:
                transmitted += 1

        elif data["result"] == "absorbed":
            absorbed += 1

    # ANALYZE
    total_escaped = reflected + transmitted
    reflection_rate = reflected / num_particles * 100
    transmission_rate = transmitted / num_particles * 100

    print(f"{'Absorbed':<20} | {'0':<15} | {absorbed:<15}")
    print(f"{'Conservation':<20} | {num_particles:<15} | {total_escaped:<15}")
    print(f"{'Reflection %':<20} | {'~80-90%':<15} | {reflection_rate:.2f}%")
    print(f"{'Transmission %':<20} | {'~10-20%':<15} | {transmission_rate:.2f}%")

    if absorbed == 0 and total_escaped == num_particles:
        print("\n[PASS] Conservation of mass holds. Physics engine is consistent.")
    else:
        print("\n[FAIL] Particles were lost or absorbed in a pure scatterer.")

if __name__ == "__main__":
    run_reflection_test()
