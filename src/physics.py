import numpy as np

# Constants
m_n = 1.674927471e-27  # Neutron mass in kg
eV_to_J = 1.60217663e-19 # Conversion factor

# ==========================================
#  VECTOR KINEMATICS HELPERS (The "Real" Physics)
# ==========================================

def get_vectors_and_vcm(initial_energy, A, sampler, rng):
    """
    Generates velocity vectors for the neutron and target,
    and calculates the Center of Mass (CM) velocity.
    """
    # 1. Neutron Scalar Speed (Lab Frame)
    v_n = np.sqrt(2 * initial_energy * eV_to_J / m_n)

    # 2. Target Scalar Speed (Sampled from Free Gas)
    # If energy is high (>10 eV), thermal motion is negligible (v_t -> 0)
    if initial_energy > 10.0:
        v_t = 0.0
    else:
        v_t = sampler.sample_velocity(vn=v_n)

    # 3. Create Vectors
    # We define the Neutron as moving along the Z-axis [0, 0, v_n]
    # The Target is sampled with an isotropic direction relative to the neutron.
    vec_vn = np.array([0.0, 0.0, v_n])

    if v_t > 0:
        # Sample random direction for target
        mu_t = 2 * rng.random() - 1
        phi_t = 2 * np.pi * rng.random()
        sin_t = np.sqrt(max(0.0, 1.0 - mu_t**2))

        vec_vt = np.array([
            v_t * sin_t * np.cos(phi_t),
            v_t * sin_t * np.sin(phi_t),
            v_t * mu_t
        ])
    else:
        vec_vt = np.zeros(3)

    # 4. Calculate Center of Mass Velocity (Vector)
    # v_cm = (v_n + A * v_t) / (A + 1)
    vec_vcm = (vec_vn + A * vec_vt) / (A + 1)

    # 5. Neutron Velocity in CM Frame (Vector)
    # V_n = v_n - v_cm
    vec_Vn = vec_vn - vec_vcm

    return vec_vn, vec_vcm, vec_Vn

def scatter_isotropic_cm(vec_Vn, magnitude, rng):
    """
    Scatters a vector isotropically in 3D, preserving a specific magnitude.
    """
    # Sample random direction in CM
    mu_cm = 2 * rng.random() - 1
    phi_cm = 2 * np.pi * rng.random()
    sin_cm = np.sqrt(max(0.0, 1.0 - mu_cm**2))

    # New direction components
    u_prime = sin_cm * np.cos(phi_cm)
    v_prime = sin_cm * np.sin(phi_cm)
    w_prime = mu_cm

    # New Vector in CM
    vec_Vn_prime = magnitude * np.array([u_prime, v_prime, w_prime])

    return vec_Vn_prime, mu_cm

# ==========================================
#  MAIN SCATTERING FUNCTIONS
# ==========================================

def elastic_scattering(initial_energy, A, sampler, rng):
    """
    Performs elastic scattering using full Vector Kinematics.
    Replaces the scalar approximation for higher fidelity.
    """
    # 1. Get Kinematics Vectors
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)

    # Magnitude of neutron velocity in CM (Conserved in Elastic)
    Vn_speed = np.linalg.norm(vec_Vn)

    # 2. Scatter in CM (Isotropic)
    vec_Vn_prime, mu_cm = scatter_isotropic_cm(vec_Vn, Vn_speed, rng)

    # 3. Transform back to Lab Frame
    # v'_n = V'_n + v_cm
    vec_vn_prime = vec_Vn_prime + vec_vcm

    # 4. Final Energy
    vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
    E_prime = 0.5 * m_n * vn_prime_speed_sq / eV_to_J

    # 5. Calculate Lab Scattering Angle (mu_lab)
    # Dot product of initial direction (0,0,1) and final direction
    # Since initial was Z-axis, mu_lab is just w component / magnitude
    vn_prime_speed = np.sqrt(vn_prime_speed_sq)

    if vn_prime_speed > 0:
        mu_lab = vec_vn_prime[2] / vn_prime_speed
    else:
        mu_lab = 1.0

    return max(1e-5, E_prime), mu_cm, max(-1.0, min(1.0, mu_lab))

def inelastic_scattering(initial_energy, A, Q_value, sampler, rng):
    """
    Performs inelastic scattering (Discrete Level).
    Q_value: The energy lost to excitation (positive float, e.g., 80000 eV).
    """
    # 1. Get Kinematics Vectors
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)

    # 2. Calculate Kinetic Energy available in CM
    Vn_speed = np.linalg.norm(vec_Vn)

    # Energy of the neutron in the CM frame
    E_n_cm_joules = 0.5 * m_n * Vn_speed**2
    E_n_cm = E_n_cm_joules / eV_to_J

    # Total Kinetic Energy in CM = E_n_cm * (A+1)/A
    Total_KE_cm = E_n_cm * (A + 1) / A

    # Check Threshold
    if Total_KE_cm <= Q_value:
        # Not enough energy to excite nucleus: Fallback to Elastic
        return elastic_scattering(initial_energy, A, sampler, rng)

    # Subtract Q (Energy Loss)
    Total_KE_cm_prime = Total_KE_cm - Q_value

    # Convert back to Neutron Velocity magnitude in CM
    # E_n_cm' = Total_KE_cm' * A / (A+1)
    E_n_cm_prime = Total_KE_cm_prime * A / (A + 1)

    Vn_prime_speed = np.sqrt(2 * E_n_cm_prime * eV_to_J / m_n)

    # 3. Scatter in CM (Isotropic) with REDUCED magnitude
    vec_Vn_prime, mu_cm = scatter_isotropic_cm(vec_Vn, Vn_prime_speed, rng)

    # 4. Transform back to Lab
    vec_vn_prime = vec_Vn_prime + vec_vcm

    # 5. Final Energy & Angle
    vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
    E_prime = 0.5 * m_n * vn_prime_speed_sq / eV_to_J

    vn_prime_speed = np.sqrt(vn_prime_speed_sq)
    if vn_prime_speed > 0:
        mu_lab = vec_vn_prime[2] / vn_prime_speed
    else:
        mu_lab = 1.0

    return max(1e-5, E_prime), mu_cm, max(-1.0, min(1.0, mu_lab))

# ==========================================
#  LEGACY / COMPATIBILITY FUNCTIONS
#  (Required for Unit Tests to pass)
# ==========================================

def calculate_mu_lab(mu_cm, E, E_prime, E_cm_prime, A):
    """
    Legacy scalar implementation.
    Kept because tests/test_physics.py imports it directly.
    """
    if E_prime <= 0: return 0.0
    term1 = mu_cm * np.sqrt(E_cm_prime / E_prime)
    term2 = (1 / (A + 1)) * np.sqrt(E / E_prime)
    mu_lab = term1 + term2
    return max(-1.0, min(1.0, mu_lab))

def calculate_E_cm_prime(initial_energy, A, sampler):
    """
    Used for Cross Section Lookup and Legacy Kinematics.
    """
    # High Energy Case (>10 eV)
    if initial_energy > 10:
        # At high energy, we must return the kinetic energy relative to the center of mass.
        # This is ALWAYS lower than the lab energy because we subtract the forward momentum.
        ratio = A / (A + 1)
        return initial_energy * (ratio ** 2)

    # Thermal Energy Case (<10 eV)
    # Estimate effective energy using simplified vector math
    v_n = np.sqrt(2 * initial_energy * eV_to_J / m_n)
    v_t = sampler.sample_velocity(vn=v_n)

    # 1. Scalar approximation of Center of Mass Velocity
    v_cm = (v_n + A * v_t) / (A + 1)

    # 2. Neutron velocity in CM frame
    v_l_cm = abs(v_n - v_cm)

    # 3. Energy in CM
    E_cm_prime_joules = 0.5 * m_n * v_l_cm**2
    return E_cm_prime_joules / eV_to_J

def calculate_E_prime(E_cm_prime, initial_energy, A, rng):
    """Legacy scalar implementation for tests."""
    E = initial_energy
    mu_cm = 2 * rng.random() - 1
    numerator = E + 2 * mu_cm * (A + 1) * np.sqrt(E * E_cm_prime)
    term2 = numerator / ((A + 1)**2)
    E_prime = E_cm_prime + term2
    return E_prime, mu_cm

def sample_new_direction_cosines(u, v, w, mu_lab, rng):
    """
    Standard rotation matrix. Used by simulation.py.
    """
    phi = 2 * np.pi * rng.random()
    sin_theta = np.sqrt(max(0.0, 1.0 - mu_lab**2))

    if abs(w) >= 0.999999:
        sign = 1.0 if w > 0 else -1.0
        u_new = sin_theta * np.cos(phi)
        v_new = sin_theta * np.sin(phi)
        w_new = sign * mu_lab
    else:
        denom = np.sqrt(max(1e-12, 1.0 - w**2))
        u_new = (mu_lab * u) + (sin_theta / denom) * (u * w * np.cos(phi) - v * np.sin(phi))
        v_new = (mu_lab * v) + (sin_theta / denom) * (v * w * np.cos(phi) + u * np.sin(phi))
        w_new = (mu_lab * w) - (sin_theta * denom * np.cos(phi))

    norm = np.sqrt(u_new**2 + v_new**2 + w_new**2)
    return u_new/norm, v_new/norm, w_new/norm, phi
