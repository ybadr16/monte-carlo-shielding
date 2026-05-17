import h5py
import os

# Adjust path if needed
base_path = "endfb/neutron"
files_to_check = [("Fe56.h5", 91), ("Be9.h5", 16)]

for fname, mt in files_to_check:
    path = os.path.join(base_path, fname)
    if not os.path.exists(path): continue

    print(f"\n--- Inspecting {fname} MT={mt} ---")
    try:
        with h5py.File(path, 'r') as f:
            # Find root
            if f"neutron/{fname.split('.')[0]}" in f: root = f[f"neutron/{fname.split('.')[0]}"]
            elif fname.split('.')[0] in f: root = f[fname.split('.')[0]]
            else: print("Root not found"); continue

            # Find Distribution
            dist_path = f"reactions/reaction_{mt:03}/product_0/distribution_0"
            if dist_path not in root: print("Distribution not found"); continue

            grp = root[dist_path]
            print(f"Group Keys: {list(grp.keys())}")

            # Check Energy Attributes
            if "energy" in grp:
                print(f"Energy Dataset Attributes: {list(grp['energy'].attrs.keys())}")

            # Check Group Attributes
            print(f"Group Attributes: {list(grp.attrs.keys())}")

            # Check if offsets exist as a dataset
            if "offsets" in grp:
                print("Found 'offsets' as a separate Dataset!")

    except Exception as e:
        print(f"Error: {e}")
