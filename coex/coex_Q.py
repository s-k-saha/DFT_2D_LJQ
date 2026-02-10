import numpy as np
from scipy.optimize import fsolve
import math

pi = np.pi
R = 0.5
R2 = R * R
R3 = R * R2

eps_si = 1.503#1e-2   # solvent-ion interaction
eps_ss = 1.0     # solvent-solvent interaction
alpha_si = 1.171861897
alpha_ss = 1.171861897

rhoig = 5.5619e-05    # fixed ion density in gas

def mu_s(rhos, rhoi):
    eta = 4./3 * pi * R3 * (rhos + 2 * rhoi)
    return (
        np.log(rhos)
        + (14*eta - 13*eta**2 + 5*eta**3) / (2*(1-eta)**3)
        - np.log(1-eta)
        + 4*pi*rhos*(-alpha_ss)*eps_ss
        + 8*pi*rhoi*(-alpha_si)*eps_si
    )

def mu_i(rhos, rhoi):
    eta = 4./3 * pi * R3 * (rhos + 2 * rhoi)
    return (
        np.log(rhoi)
        + (14*eta - 13*eta**2 + 5*eta**3) / (2*(1-eta)**3)
        - np.log(1-eta)
        + 4*pi*rhos*(-alpha_si)*eps_si
    )

def P(rhos, rhoi):
    eta = 4./3 * pi * R3 * (rhos + 2 * rhoi)
    return (
        (rhos + 2*rhoi) * (1 + eta + eta**2) / (1-eta)**3
        + 8*pi*rhoi*rhos*(-alpha_si)*eps_si
        + 2*pi*rhos*rhos*(-alpha_ss)*eps_ss
    )

# ---- Modified system: use log-density for ions in liquid ----

def system(vars):
    rho_s_l, rho_s_g, x_i_l = vars
    rho_i_l = np.exp(x_i_l)   # <-- key fix
    rho_i_g = rhoig

    return [
        mu_s(rho_s_l, rho_i_l) - mu_s(rho_s_g, rho_i_g),
        mu_i(rho_s_l, rho_i_l) - mu_i(rho_s_g, rho_i_g),
        P(rho_s_l, rho_i_l)    - P(rho_s_g, rho_i_g)
    ]

if __name__ == "__main__":

    # Initial guesses
    rho_s_l0 = 0.5331
    rho_s_g0 = 0.0113
    rho_i_l0 = 7.8771e-02

    init_guess = [
        rho_s_l0,
        rho_s_g0,
        np.log(rho_i_l0)   # <-- log variable
    ]

    sol = fsolve(system, init_guess)

    rho_s_l, rho_s_g, x_i_l = sol
    rho_i_l = np.exp(x_i_l)
    print(
        eps_ss, eps_si,
        rho_s_l, rho_s_g, rho_i_l,rhoig,
        mu_s(rho_s_l, rho_i_l),
        mu_s(rho_s_g, rhoig),
        mu_i(rho_s_l, rho_i_l),
        mu_i(rho_s_g, rhoig),
        P(rho_s_l, rho_i_l),
        P(rho_s_g, rhoig)
    )

