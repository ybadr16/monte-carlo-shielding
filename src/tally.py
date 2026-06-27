# src/tally.py
import numpy as np

from .statistics import escape_statistics


class Tally:
    def __init__(self):
        self.results = {"absorbed": 0.0, "escaped": 0, "killed": 0}
        self.absorbed_coordinates = []
        self.energy_spectrum = []        # per-history escape energy (eV)
        self.leaked_weights = []         # per-history leaked weight (0 if none)
        self.region_count = 0            # particles detected in a specific region

    def update(self, result, absorbed_weight=0.0, absorbed_coords=None,
               final_energy=None, final_weight=0.0, region_detected=False):
        # Discrete-event counters (analog mode)
        if result in ["escaped", "killed"]:
            self.results[result] += 1

        # Weight accounting
        self.results["absorbed"] += absorbed_weight
        self.leaked_weights.append(final_weight)

        if absorbed_coords:
            self.absorbed_coordinates.extend(absorbed_coords)

        if final_energy is not None:
            self.energy_spectrum.append(final_energy)

        if region_detected:
            self.region_count += 1

    def merge_partial_results(self, partial_results):
        """Merge one source-history result dict from a worker process."""
        res_type = partial_results["result"]
        if res_type in ["escaped", "killed"]:
            self.results[res_type] += 1

        self.results["absorbed"] += partial_results["absorbed_weight"]
        # kept in source order for the batch-means uncertainties
        self.leaked_weights.append(partial_results.get("final_weight", 0.0))

        if partial_results["absorbed_coords"]:
            self.absorbed_coordinates.extend(partial_results["absorbed_coords"])

        if partial_results["final_energy"] is not None:
            self.energy_spectrum.append(partial_results["final_energy"])

        if partial_results["region_detected"]:
            self.region_count += 1

    def statistics(self, num_particles, n_batches=10):
        """Batch-means leakage fraction and average escape energy, each with a
        standard error.  Returns the dict from ``escape_statistics``."""
        n = len(self.leaked_weights)
        if n == 0:
            return {'leakage': 0.0, 'leakage_sem': float('nan'),
                    'avg_energy': 0.0, 'avg_energy_sem': float('nan')}
        weights = np.asarray(self.leaked_weights, dtype=float)
        # Energies align with leaked weights when both are recorded per history.
        if len(self.energy_spectrum) == n:
            energies = np.asarray(self.energy_spectrum, dtype=float)
        else:
            energies = np.zeros(n)
        return escape_statistics(weights, energies, num_particles, n_batches)

    def print_summary(self, num_particles, n_batches=10):
        stats = self.statistics(num_particles, n_batches)

        def pm(value, sem, fmt):
            sem_str = format(sem, fmt) if np.isfinite(sem) else 'n/a'
            return f"{format(value, fmt)} ± {sem_str}"

        print("--- Simulation Results ---")
        print(f"  Total Particles Started:    {num_particles}")
        # Leakage is weight-based (the meaningful quantity under implicit
        # capture, where histories are not simply 'escaped'/'absorbed').
        print(f"  Leakage Fraction (weight):  {pm(stats['leakage'], stats['leakage_sem'], '.5f')}")
        print(f"  Total Weight Absorbed:      {self.results['absorbed']:.4f}")
        print(f"  Absorption Events Recorded: {len(self.absorbed_coordinates)}")
        if self.results['escaped'] or self.results['killed']:
            print(f"  Analog escaped / killed:    "
                  f"{self.results['escaped']} / {self.results['killed']}")
        print(f"  Detected in region:         {self.region_count}")

        if stats['avg_energy'] > 0:
            print(f"  Average escape energy:      "
                  f"{pm(stats['avg_energy'], stats['avg_energy_sem'], '.3e')} eV")
        else:
            print("  No leaked weight — average escape energy undefined.")

    def get_results(self):
        return {
            "results": self.results,
            "absorbed_coordinates": self.absorbed_coordinates,
            "energy_spectrum": self.energy_spectrum,
            "leaked_weights": self.leaked_weights,
            "region_count": self.region_count,
        }
