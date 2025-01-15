# tally.py
from geometry import count_coordinates_in_boundary

class Tally:
    def __init__(self):
        # Initialize counts and data storage
        self.results = {"absorbed": 0, "escaped": 0}
        self.absorbed_coordinates = []
        self.energy_spectrum = []
        self.region_count = 0  # Tracks particles detected in a specific region

    def update(self, result, absorbed=None, final_energy=None, region_detected=False):
        # Update tally counts
        self.results[result] += 1

        # Track absorbed and fission coordinates
        if absorbed:
            self.absorbed_coordinates.extend(absorbed)

        # Track energy spectrum
        if final_energy is not None:
            self.energy_spectrum.append(final_energy)

        # Update region detection count
        if region_detected:
            self.region_count += 1

    def get_results(self):
        # Return a summary of the results
        return {
            "results": self.results,
            "absorbed_coordinates": self.absorbed_coordinates,
            "energy_spectrum": self.energy_spectrum,
            "region_count": self.region_count,
        }

    def print_summary(self, num_particles):
        # Print a summary of results
        print(f"Simulation Results:")
        print(f"  Absorbed: {self.results['absorbed']}")
        print(f"  Escaped: {self.results['escaped']}")
        print(f"  Detected within detection region (if specified): {self.region_count}")
        print(f"  Total particles simulated: {num_particles}")
        if self.energy_spectrum:
            avg_energy = sum(self.energy_spectrum) / len(self.energy_spectrum)
            print(f"  Average final energy: {avg_energy:.2f} eV")
        else:
            print(f"  No particles left to calculate average final energy.")
        #print(f"  Absorbed coordinates: {self.absorbed_coordinates[:10]}...")  # First 10 for brevity

    def merge_partial_results(self, partial_results):
        """
        Merge partial results from a worker process.
        """
        self.results[partial_results["result"]] += 1

        if partial_results["absorbed"]:
            self.absorbed_coordinates.extend(partial_results["absorbed"])
        if partial_results["final_energy"] is not None:
            self.energy_spectrum.append(partial_results["final_energy"])
        if partial_results["region_detected"]:
            self.region_count += 1

