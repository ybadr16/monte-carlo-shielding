import numpy as np

m_n = 1.674927471e-27  # neutron mass, kg
eV_to_J = 1.60217663e-19


def get_vectors_and_vcm(initial_energy, A, sampler, rng):
    """Lab/CM velocity vectors for a neutron on a free-gas target."""
    v_n = np.sqrt(2 * initial_energy * eV_to_J / m_n)
    vec_vn = np.array([0.0, 0.0, v_n])

    # Target motion is negligible above ~10 eV. Below it the target speed and
    # its cosine to the neutron are sampled together, keeping the speed-angle
    # correlation that lets the neutron upscatter toward thermal equilibrium.
    if initial_energy > 10.0:
        vec_vt = np.zeros(3)
    else:
        v_t, mu_t = sampler.sample_velocity(vn=v_n, return_mu=True)
        if v_t > 0:
            phi_t = 2 * np.pi * rng.random()
            sin_t = np.sqrt(max(0.0, 1.0 - mu_t**2))
            vec_vt = np.array([
                v_t * sin_t * np.cos(phi_t),
                v_t * sin_t * np.sin(phi_t),
                v_t * mu_t,
            ])
        else:
            vec_vt = np.zeros(3)

    vec_vcm = (vec_vn + A * vec_vt) / (A + 1)
    vec_Vn = vec_vn - vec_vcm
    return vec_vn, vec_vcm, vec_Vn

def scatter_isotropic_cm(vec_Vn, magnitude, rng):
    mu_cm = 2 * rng.random() - 1
    phi_cm = 2 * np.pi * rng.random()
    sin_cm = np.sqrt(max(0.0, 1.0 - mu_cm**2))

    u_prime = sin_cm * np.cos(phi_cm)
    v_prime = sin_cm * np.sin(phi_cm)
    w_prime = mu_cm

    vec_Vn_prime = magnitude * np.array([u_prime, v_prime, w_prime])
    return vec_Vn_prime, mu_cm

def elastic_scattering(initial_energy, A, sampler, rng):
    """Free-gas elastic scattering, isotropic in CM.

    The transport loop draws the elastic CM cosine from the tabulated ENDF
    angular data instead; this self-contained kernel is used by the discrete
    fallbacks and the unit tests.
    """
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)
    Vn_speed = np.linalg.norm(vec_Vn)

    vec_Vn_prime, mu_cm = scatter_isotropic_cm(vec_Vn, Vn_speed, rng)
    vec_vn_prime = vec_Vn_prime + vec_vcm

    vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
    E_prime = 0.5 * m_n * vn_prime_speed_sq / eV_to_J

    vn_prime_speed = np.sqrt(vn_prime_speed_sq)
    if vn_prime_speed > 0:
        mu_lab = vec_vn_prime[2] / vn_prime_speed
    else:
        mu_lab = 1.0

    return max(1e-5, E_prime), mu_cm, max(-1.0, min(1.0, mu_lab))

def inelastic_scattering(initial_energy, A, Q_value, sampler, rng):
    """Discrete-level inelastic scattering, isotropic in CM."""
    vec_vn, vec_vcm, vec_Vn = get_vectors_and_vcm(initial_energy, A, sampler, rng)
    Vn_speed = np.linalg.norm(vec_Vn)

    E_n_cm_joules = 0.5 * m_n * Vn_speed**2
    E_n_cm = E_n_cm_joules / eV_to_J
    Total_KE_cm = E_n_cm * (A + 1) / A

    if Total_KE_cm <= Q_value:
        return elastic_scattering(initial_energy, A, sampler, rng)

    Total_KE_cm_prime = Total_KE_cm - Q_value
    E_n_cm_prime = Total_KE_cm_prime * A / (A + 1)
    Vn_prime_speed = np.sqrt(2 * E_n_cm_prime * eV_to_J / m_n)

    vec_Vn_prime, mu_cm = scatter_isotropic_cm(vec_Vn, Vn_prime_speed, rng)

    vec_vn_prime = vec_Vn_prime + vec_vcm
    vn_prime_speed_sq = np.dot(vec_vn_prime, vec_vn_prime)
    E_prime = 0.5 * m_n * vn_prime_speed_sq / eV_to_J

    vn_prime_speed = np.sqrt(vn_prime_speed_sq)
    mu_lab = vec_vn_prime[2] / vn_prime_speed if vn_prime_speed > 0 else 1.0

    return max(1e-5, E_prime), mu_cm, max(-1.0, min(1.0, mu_lab))

def sample_maxwellian(temperature_ev, rng):
    """Evaporation/fission energy from f(E) ~ sqrt(E) exp(-E/T) (algorithm C64)."""
    r1 = rng.random()
    r2 = rng.random()
    r3 = rng.random()

    val = -np.log(r1) - np.log(r2) * (np.cos(np.pi * r3 / 2.0) ** 2)
    return temperature_ev * val

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
    """Fermi-gas nuclear temperature, T = sqrt(U/a) with a = A/8 MeV^-1. energy, T in eV."""
    if energy <= 0:
        return 1e3

    a = A / 8.0
    energy_MeV = energy / 1e6
    T_MeV = np.sqrt(energy_MeV / a)
    return T_MeV * 1e6
