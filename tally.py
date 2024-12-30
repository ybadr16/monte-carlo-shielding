# tally.py
class Tally:
    def __init__(self):
        # Initialize counts and data storage
        self.results = {"absorbed": 0, "fission": 0, "escaped": 0}
        self.absorbed_coordinates = []
        self.fission_coordinates = []
        self.energy_spectrum = []
        self.region_count = 0  # Tracks particles detected in a specific region

    def update(self, result, absorbed=None, fissioned=None, final_energy=None, region_detected=False):
        # Update tally counts
        self.results[result] += 1

        # Track absorbed and fission coordinates
        if absorbed:
            self.absorbed_coordinates.extend(absorbed)
        if fissioned:
            self.fission_coordinates.extend(fissioned)

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
            "fission_coordinates": self.fission_coordinates,
            "energy_spectrum": self.energy_spectrum,
            "region_count": self.region_count,
        }

    def print_summary(self, num_particles):
        # Print a summary of results
        print(f"Simulation Results:")
        print(f"  Absorbed: {self.results['absorbed']}")
        print(f"  Fission: {self.results['fission']}")
        print(f"  Escaped: {self.results['escaped']}")
        print(f"  Detected within detection region (if specified): {self.region_count}")
        print(f"  Total particles simulated: {num_particles}")
        if self.energy_spectrum:
            avg_energy = sum(self.energy_spectrum) / len(self.energy_spectrum)
            print(f"  Average final energy: {avg_energy:.2f} eV")
        else:
            print(f"  No particles left to calculate average final energy.")
        print(f"  Absorbed coordinates: {self.absorbed_coordinates[:10]}...")  # First 10 for brevity
        print(f"  Fission coordinates: {self.fission_coordinates[:10]}...")  # First 10 for brevity
        # Add this method to Tally
    def count_coordinates_in_boundary(self, coordinates, x_bounds, y_bounds, z_bounds):
        x_min, x_max = x_bounds
        y_min, y_max = y_bounds
        z_min, z_max = z_bounds
        return [
            coord for coord in coordinates
            if x_min <= coord[0] <= x_max and y_min <= coord[1] <= y_max and z_min <= coord[2] <= z_max
        ]
