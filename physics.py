# physics.py
import numpy as np

def calculate_mu_lab(mu_cm, E, E_prime, E_cm_prime, A):
    term1 = mu_cm * np.sqrt(E_cm_prime / E_prime)
    term2 = 1 / (A + 1) * np.sqrt(E / E_prime)
    mu_lab = term1 + term2
    #print(mu_lab)
    return mu_lab


def calculate_E_cm_prime(initial_energy, A, sampler):
    m_n = 1.674927471 * 10**-27  # neutron mass in kg
    E = initial_energy
    v_l = np.sqrt(2 * E / m_n)  # neutron velocity in the lab frame

    # Sample target velocity if energy is below a threshold
    v_t = sampler.sample_velocity(vn=v_l) if E < 10 else 0

    # Calculate center-of-mass velocity and neutron velocity in the CM frame
    v_cm = (v_l + A * v_t) / (A + 1)
    v_l_cm = v_l - v_cm

    # Calculate E_cm_prime
    E_cm_prime = 0.5 * m_n * v_l_cm**2
    return E_cm_prime

def calculate_E_prime(E_cm_prime, initial_energy, A, rng):
    E = initial_energy
    mu_cm = 2 * rng.uniform(0, 1) - 1  # Isotropic distribution in CM frame

    # Calculate E_prime (neutron energy in the lab frame after scattering)
    E_prime = E_cm_prime + (E + 2 * mu_cm * (A + 1) * np.sqrt(E * E_cm_prime)) / ((A + 1)**2)
    return E_prime, mu_cm

# Example of how to use these functions
def elastic_scattering(initial_energy, A, sampler, rng):
    E_cm_prime = calculate_E_cm_prime(initial_energy, A, sampler)
    E_prime, mu_cm = calculate_E_prime(E_cm_prime, initial_energy, A, rng)

    mu_lab = calculate_mu_lab(mu_cm, initial_energy, E_prime, E_cm_prime, A)
    return E_prime, mu_cm, mu_lab


def sample_new_direction_cosines(u, v, w, mu_lab, rng):
    phi = 2 * np.pi * rng.random()
    root_term = np.sqrt(max(0, 1 - mu_lab**2))
    denom = np.sqrt(max(1e-8, 1 - w**2))
    u_new = mu_lab * u + root_term * (u * w * np.cos(phi) - v * np.sin(phi)) / denom
    v_new = mu_lab * v + root_term * (v * w * np.cos(phi) + u * np.sin(phi)) / denom
    w_new = mu_lab * w - root_term * np.cos(phi) * denom
    assert abs(u_new**2 + v_new**2 + w_new**2 - 1.0) < 1e-6
    return u_new, v_new, w_new, phi
