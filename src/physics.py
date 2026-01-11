import numpy as np

# Constants
m_n = 1.674927471e-27  # Neutron mass in kg
eV_to_J = 1.60217663e-19 # Conversion factor

# ==========================================
#  VECTOR KINEMATICS HELPERS
# ==========================================

def get_vectors_and_vcm(initial_energy, A, sampler, rng):
    """
    Generates velocity vectors for the neutron and target,
    and calculates the Center of Mass (CM) velocity.
    """
    # 1. Neutron Scalar Speed (Lab Frame)
    v_n = np.sqrt(2 * initial_energy * eV_to_J / m_n)

    # 2. Target Scalar Speed (Sampled from Free Gas)
    if initial_energy > 10.0:
        v_t = 0.0
    else:
        v_t = sampler.sample_velocity(vn=v_n)

    # 3. Create Vectors
    vec_vn = np.array([0.0, 0.0, v_n])

    if v_t > 0:
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

    # 4. Calculate Center of Mass Velocity
    vec_vcm = (vec_vn + A * vec_vt) / (A + 1)

    # 5. Neutron Velocity in CM Frame
    vec_Vn = vec_vn - vec_vcm

    return vec_vn, vec_vcm, vec_Vn

def scatter_isotropic_cm(vec_Vn, magnitude, rng):
    """
    Scatters a vector isotropically in 3D.
    """
    mu_cm = 2 * rng.random() - 1
    phi_cm = 2 * np.pi * rng.random()
    sin_cm = np.sqrt(max(0.0, 1.0 - mu_cm**2))

    u_prime = sin_cm * np.cos(phi_cm)
    v_prime = sin_cm * np.sin(phi_cm)
    w_prime = mu_cm

    vec_Vn_prime = magnitude * np.array([u_prime, v_prime, w_prime])
    return vec_Vn_prime, mu_cm

def scatter_anisotropic_cm(vec_Vn, magnitude, rng, mu_bar=0.0):
    """
    Scatters with a forward bias (Linear Anisotropy).
    mu_bar: Average scattering cosine (0=Isotropic, 1=Forward).
    """
    if abs(mu_bar) < 1e-3:
        return scatter_isotropic_cm(vec_Vn, magnitude, rng)

    # Sample mu from P(mu) = 1 + 3*mu_bar*mu (Linear Approximation)
    # Rejection sampling method:
    while True:
        mu = 2 * rng.random() - 1
        prob = 0.5 * (1 + 3 * mu_bar * mu)
        max_p = 0.5 * (1 + 3 * abs(mu_bar))
        if rng.random() * max_p < prob:
            mu_cm = mu
            break

    phi_cm = 2 * np.pi * rng.random()
    sin_cm = np.sqrt(max(0.0, 1.0 - mu_cm**2))

    u_prime = sin_cm * np.cos(phi_cm)
    v_prime = sin_cm * np.sin(phi_cm)
    w_prime = mu_cm

    vec_Vn_prime = magnitude * np.array([u_prime, v_prime, w_prime])
    return vec_Vn_prime, mu_cm

# ==========================================
#  MAIN SCATTERING FUNCTIONS
# ==========================================

def elastic_scattering(initial_energy, A, sampler, rng):
    """
    Performs elastic scattering using full Vector Kinematics + Anisotropy.
    """
    # 1. Get Kinematics Vectors
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)
    Vn_speed = np.linalg.norm(vec_Vn)

    # 2. ANISOTROPY CHECK
    # Lead (A=208) scatters forward at high energy.
    # We use a linear ramp approximation for mu_bar based on energy.
    mu_bar = 0.0
    if A > 50 and initial_energy > 50e3: # Only for heavy nuclei > 50 keV
        # Logarithmic ramp: 0.05 at 100keV -> 0.6 at 14MeV
        mu_bar = min(0.6, 0.05 * np.log10(initial_energy/1e4))
        mu_bar = max(0.0, mu_bar)

    # 3. Scatter in CM (Anisotropic)
    vec_Vn_prime, mu_cm = scatter_anisotropic_cm(vec_Vn, Vn_speed, rng, mu_bar)

    # 4. Transform back to Lab Frame
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

def inelastic_scattering(initial_energy, A, Q_value, sampler, rng):
    """
    Performs discrete inelastic scattering.
    """
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)
    Vn_speed = np.linalg.norm(vec_Vn)

    # Energy Analysis in CM
    E_n_cm_joules = 0.5 * m_n * Vn_speed**2
    E_n_cm = E_n_cm_joules / eV_to_J
    Total_KE_cm = E_n_cm * (A + 1) / A

    if Total_KE_cm <= Q_value:
        return elastic_scattering(initial_energy, A, sampler, rng)

    # Subtract Energy Cost
    Total_KE_cm_prime = Total_KE_cm - Q_value
    E_n_cm_prime = Total_KE_cm_prime * A / (A + 1)
    Vn_prime_speed = np.sqrt(2 * E_n_cm_prime * eV_to_J / m_n)

    # Inelastic is Isotropic in CM (Approximation)
    vec_Vn_prime, mu_cm = scatter_isotropic_cm(vec_Vn, Vn_prime_speed, rng)

    # Transform back
    vec_vn_prime = vec_Vn_prime + vec_vcm
    vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
    E_prime = 0.5 * m_n * vn_prime_speed_sq / eV_to_J

    vn_prime_speed = np.sqrt(vn_prime_speed_sq)
    mu_lab = vec_vn_prime[2] / vn_prime_speed if vn_prime_speed > 0 else 1.0

    return max(1e-5, E_prime), mu_cm, max(-1.0, min(1.0, mu_lab))

def sample_maxwellian(temperature_ev, rng):
    """
    Samples energy from Maxwellian distribution: f(E) ~ sqrt(E)*exp(-E/T)
    Used for Evaporation and Fission.
    """
    r1 = rng.random()
    r2 = rng.random()
    r3 = rng.random()

    # Rejectionless sampling (Rule C64)
    val = -np.log(r1) - np.log(r2) * (np.cos(np.pi * r3 / 2.0) ** 2)
    return temperature_ev * val

# ==========================================
#  LEGACY / COMPATIBILITY (For Tests)
# ==========================================

def calculate_mu_lab(mu_cm, E, E_prime, E_cm_prime, A):
    if E_prime <= 0: return 0.0
    term1 = mu_cm * np.sqrt(E_cm_prime / E_prime)
    term2 = (1 / (A + 1)) * np.sqrt(E / E_prime)
    mu_lab = term1 + term2
    return max(-1.0, min(1.0, mu_lab))

def calculate_E_cm_prime(initial_energy, A, sampler):
    if initial_energy > 10:
        ratio = A / (A + 1)
        return initial_energy * (ratio ** 2)
    v_n = np.sqrt(2 * initial_energy * eV_to_J / m_n)
    v_t = sampler.sample_velocity(vn=v_n)
    v_cm = (v_n + A * v_t) / (A + 1)
    v_l_cm = abs(v_n - v_cm)
    return (0.5 * m_n * v_l_cm**2) / eV_to_J

def calculate_E_prime(E_cm_prime, initial_energy, A, rng):
    E = initial_energy
    mu_cm = 2 * rng.random() - 1
    numerator = E + 2 * mu_cm * (A + 1) * np.sqrt(E * E_cm_prime)
    term2 = numerator / ((A + 1)**2)
    E_prime = E_cm_prime + term2
    return E_prime, mu_cm

def sample_new_direction_cosines(u, v, w, mu_lab, rng):
    phi = 2 * np.pi * rng.random()
    sin_theta = np.sqrt(max(0.0, 1.0 - mu_lab**2))
    if abs(w) >= 0.999999:
        sign = 1.0 if w > 0 else -1.0
        return sin_theta * np.cos(phi), sin_theta * np.sin(phi), sign * mu_lab, phi
    else:
        denom = np.sqrt(max(1e-12, 1.0 - w**2))
        u_new = (mu_lab * u) + (sin_theta / denom) * (u * w * np.cos(phi) - v * np.sin(phi))
        v_new = (mu_lab * v) + (sin_theta / denom) * (v * w * np.cos(phi) + u * np.sin(phi))
        w_new = (mu_lab * w) - (sin_theta * denom * np.cos(phi))
        norm = np.sqrt(u_new**2 + v_new**2 + w_new**2)
        return u_new/norm, v_new/norm, w_new/norm, phi

def get_nuclear_temperature(energy, A):
    """
    Calculates Nuclear Temperature (T) using the Fermi Gas Model.
    Approximation: a = A / 8.0 (Level density parameter)
    """
    if energy <= 0: return 1.0  # Safety floor

    # Level density parameter
    a = A / 8.0

    # T = sqrt(U/a) where U is excitation energy ~ incident energy
    return np.sqrt(energy / a)
