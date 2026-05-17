import h5py
import numpy as np
import os

BASE_PATH = "endfb"
TARGET = "Fe56"

def check_high_energy_structure():
    fpath = os.path.join(BASE_PATH, "neutron", f"{TARGET}.h5")
    if not os.path.exists(fpath): return

    with h5py.File(fpath, 'r') as f:
        # Navigate to distribution
        root = f"neutron/{TARGET}" if f"neutron/{TARGET}" in f else TARGET
        base = f"{root}/reactions/reaction_016/product_0/distribution_0"

        # Load Incident Energy to find the 14 MeV index
        E_in = f[base]["energy"][:]
        offsets = f[base]["offsets"] if "offsets" in f[base] else f[base]["energy_out"].attrs["offsets"]

        # Find 14 MeV bin
        idx = np.searchsorted(E_in, 14.0e6)
        # Check surrounding bins to be safe
        idx = min(idx, len(E_in) - 2)

        print(f"📊 Analyzing {TARGET} at E_in = {E_in[idx]/1e6:.1f} MeV (Index {idx})")

        start = offsets[idx]
        end = offsets[idx+1]

        # Load the Energy Out block for this high-energy neutron
        dset_out = f[base]["energy_out"]
        block = dset_out[:, start:end]

        print(f"   Block Size: {block.shape[1]} outgoing energy points")

        # Inspect Row 3 and Row 4 for the HIGHEST outgoing energy (most anisotropic)
        # Usually the last points in the block
        print("\n   Checking last 5 outgoing energy points (Highest E'):")
        print(f"   {'E_out':<12} | {'Row 3 (Val)':<12} | {'Row 4 (Val)':<12}")
        print("-" * 45)

        for i in range(max(0, block.shape[1]-5), block.shape[1]):
            e_val = block[0, i]
            r3    = int(block[3, i])
            r4    = int(block[4, i])
            print(f"   {e_val:<12.2e} | {r3:<12} | {r4:<12}")

        print("\n   [ANALYSIS]")
        r3_max = np.max(block[3, :])
        if r3_max > 5:
            print(f"   🚨 ROW 3 CHANGES! Max value is {r3_max}. It is the LENGTH.")
        else:
            print(f"   ✅ Row 3 stays small ({r3_max}). It is likely the FLAG.")

check_high_energy_structure()
