import numpy as np
from src.simulation import simulate_single_particle
from src.medium import Region, Plane
from src.settings import Settings
from src.random_number_generator import RNGHandler

# --- 1. MOCK READER (50% Scatter, 50% Absorb) ---
class MockReader:
    def get_cross_sections(self, element, energy, sampler, N):
        # Sigma_t = 1.0
        # Sigma_s = 0.5 (Scattering)
        # Sigma_a = 0.5 (Absorption)
        return 0.5, 0.5, 0.0, 1.0

# Mock Sampler (Target at rest)
class MockSampler:
    def sample_velocity(self, vn=None, rng=None, return_mu=False):
        return (0.0, 0.0) if return_mu else 0.0

def run_consistency_test():
    print("Running Implicit vs Analog Consistency Validation...")
    print("Comparing transmission through a 50% absorbing slab.")

    T = 2.0  # 2 Mean Free Paths
    N_particles = 100000

    # --- SETUP GEOMETRY ---
    slab = Region(
        surfaces=[
            Plane(0,0,-1, 0), Plane(0,0,1, T),
            Plane(-1,0,0, 50), Plane(1,0,0, 50),
            Plane(0,-1,0, 50), Plane(0,1,0, 50)
        ], name="Half-Absorbium", priority=1, element="Mix"
    )
    mediums = [slab]

    # ---------------------------------------------------------
    # RUN 1: ANALOG (CRITICALITY)
    # ---------------------------------------------------------
    print(f"\n[Run 1] Analog Mode (N={N_particles})...")
    settings_analog = Settings(mode="criticality")

    count_escaped_analog = 0

    for i in range(N_particles):
        rng = RNGHandler(seed=i)
        # Start inside, pointing forward
        state = {
            "x": 0.0, "y": 0.0, "z": 0.0001,
            "theta": 0.0, "phi": 0.0,
            "has_interacted": False, "energy": 1e6, "weight": 1.0
        }

        args = (state, MockReader(), mediums, 1.0, 1.0, MockSampler(), None, False, rng, settings_analog)
        data = simulate_single_particle(args)

        if data["result"] == "escaped":
            # In Analog, every survivor counts as 1.0
            count_escaped_analog += 1

    transmission_analog = count_escaped_analog / N_particles

    # ---------------------------------------------------------
    # RUN 2: IMPLICIT CAPTURE (SHIELDING)
    # ---------------------------------------------------------
    print(f"[Run 2] Implicit Mode (N={N_particles})...")
    settings_implicit = Settings(mode="shielding")
    # Disable roulette for this pure comparison to avoid noise
    settings_implicit.weight_cutoff = 0.0

    weight_escaped_implicit = 0.0

    for i in range(N_particles):
        rng = RNGHandler(seed=i + 10000) # Different seed sequence
        state = {
            "x": 0.0, "y": 0.0, "z": 0.0001,
            "theta": 0.0, "phi": 0.0,
            "has_interacted": False, "energy": 1e6, "weight": 1.0
        }

        args = (state, MockReader(), mediums, 1.0, 1.0, MockSampler(), None, False, rng, settings_implicit)
        data = simulate_single_particle(args)

        if data["result"] == "escaped":
            # In Implicit, we sum the WEIGHT of the survivor
            weight_escaped_implicit += data["final_weight"]

    transmission_implicit = weight_escaped_implicit / N_particles

    # ---------------------------------------------------------
    # COMPARISON
    # ---------------------------------------------------------
    print("-" * 50)
    print(f"{'Metric':<20} | {'Analog':<12} | {'Implicit':<12}")
    print("-" * 50)
    print(f"{'Transmission':<20} | {transmission_analog:<12.4f} | {transmission_implicit:<12.4f}")

    diff = abs(transmission_analog - transmission_implicit)
    print(f"\nDifference: {diff:.4f}")

    # Statistical tolerance (approx 2-3 sigma)
    # For N=5000, error is roughly 1-2%
    if diff < 0.02:
        print("[PASS] Results match within statistical error.")
        print("       Implicit Capture logic is valid.")
    else:
        print("[FAIL] Results diverge. Check weighting logic.")

if __name__ == "__main__":
    run_consistency_test()
