# physics.py
import math
import random
#import numpy as np

def calculate_mu_lab(mu_cm, E, E_prime, A):
    term1 = mu_cm * math.sqrt(E_prime / E)
    term2 = 1 / (A + 1) * math.sqrt(E / E_prime)
    mu_lab = term1 + term2
    return max(min(mu_lab, 1.0), -1.0)


def calculate_E_cm_prime(initial_energy, A, sampler):
    m_n = 1.674927471 * 10**-27  # neutron mass in kg
    E = initial_energy
    v_l = math.sqrt(2 * E / m_n)  # neutron velocity in the lab frame

    # Sample target velocity if energy is below a threshold
    v_t = sampler.sample_velocity(vn=v_l) if E < 10 else 0

    # Calculate center-of-mass velocity and neutron velocity in the CM frame
    v_cm = (v_l + A * v_t) / (A + 1)
    v_l_cm = v_l - v_cm

    # Calculate E_cm_prime
    E_cm_prime = 0.5 * m_n * v_l_cm**2
    return E_cm_prime

def calculate_E_prime(E_cm_prime, initial_energy, A):
    E = initial_energy
    mu_cm = 2 * random.uniform(0, 1) - 1  # Isotropic distribution in CM frame

    # Calculate E_prime (neutron energy in the lab frame after scattering)
    E_prime = E_cm_prime + (E + 2 * mu_cm * (A + 1) * math.sqrt(E * E_cm_prime)) / ((A + 1)**2)
    return E_prime, mu_cm

# Example of how to use these functions
def elastic_scattering(initial_energy, A, sampler):
    E_cm_prime = calculate_E_cm_prime(initial_energy, A, sampler)
    E_prime, mu_cm = calculate_E_prime(E_cm_prime, initial_energy, A)

    # Assuming calculate_mu_lab is a defined function
    mu_lab = calculate_mu_lab(mu_cm, initial_energy, E_prime, A)
    return E_prime, mu_cm, mu_lab


def sample_new_direction_cosines(u, v, w, mu_lab):
    mu_lab = max(min(mu_lab, 1.0), -1.0)
    phi = 2 * math.pi * random.random()
    root_term = math.sqrt(max(0, 1 - mu_lab**2))
    denom = math.sqrt(max(1e-8, 1 - w**2))
    u_new = mu_lab * u + root_term * (u * w * math.cos(phi) - v * math.sin(phi)) / denom
    v_new = mu_lab * v + root_term * (v * w * math.cos(phi) + u * math.sin(phi)) / denom
    w_new = mu_lab * w - root_term * math.cos(phi) / denom
    return u_new, v_new, w_new, phi, mu_lab
