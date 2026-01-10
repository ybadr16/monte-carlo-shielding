from src.cross_section_read import CrossSectionReader
import numpy as np

def check_levels():
    base_path = "./endfb"
    reader = CrossSectionReader(base_path)

    # Change this to "Pb208" or "Fe56" if you have them!
    # B10 might not have many discrete inelastic levels in some libraries.
    element_to_check = "Pb208"

    print(f"--- Scanning {element_to_check} for Inelastic Levels (MT 51-91) ---")

    # Trigger the scan
    # We use a high energy (14 MeV) to ensure we are above thresholds
    test_energy = 14.0e6
    N = 1.0 # Dummy density

    total_inelastic, components = reader.get_inelastic_components(element_to_check, test_energy, N)

    print(f"Total Inelastic Macroscopic XS at {test_energy/1e6} MeV: {total_inelastic:.4e}")

    if not components:
        print("RESULT: No inelastic levels found. (This might be normal for light isotopes like B10)")
    else:
        print(f"RESULT: Found {len(components)} discrete levels!")
        print(f"{'MT':<5} | {'XS (cm-1)':<12} | {'Q-Value (eV)':<15}")
        print("-" * 40)

        for mt, xs, q in components:
            print(f"{mt:<5} | {xs:.4e}   | {q:.4e}")

if __name__ == "__main__":
    check_levels()
