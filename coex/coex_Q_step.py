import numpy as np
from scipy.optimize import fsolve
import math

pi = np.pi
R = 0.5
R3 = R**3

# Free energy parameters
eps_ss = 1.0          # solvent-solvent attraction
alpha_ss = 1.171861897
alpha_si = 1.171861897

# Target chemical potential of ions (semi-grand)
mu_i_target = -12.0   # set low initially for dilute ions, can adjust

def mu_s(rhos, rhoi, eps_si):
    """Solvent chemical potential"""
    eta = 4./3 * pi * R3 * (rhos + 2*rhoi)
    return (
        np.log(rhos)
        + (14*eta - 13*eta**2 + 5*eta**3) / (2*(1-eta)**3)
        - np.log(1-eta)
        + 4*pi*rhos*(-alpha_ss)*eps_ss
        + 8*pi*rhoi*(-alpha_si)*eps_si
    )

def mu_i_ex(rhos, eps_si):
    """Excess part of ion chemical potential (independent of rho_i)"""
    eta = 4./3 * pi * R3 * (0 + 2*0)  # rho_i negligible here
    # Actually, we can include HS contribution with rho_i small
    eta = 4./3 * pi * R3 * rhos
    return (
        (14*eta - 13*eta**2 + 5*eta**3) / (2*(1-eta)**3)
        - np.log(1-eta)
        + 4*pi*rhos*(-alpha_si)*eps_si
    )

def P(rhos, rhoi, eps_si):
    eta = 4./3 * pi * R3 * (rhos + 2*rhoi)
    return (
        (rhos + 2*rhoi) * (1 + eta + eta**2) / (1-eta)**3
        + 8*pi*rhoi*rhos*(-alpha_si)*eps_si
        + 2*pi*rhos**2*(-alpha_ss)*eps_ss
    )

# ------------------------
# System to solve: equality of mu_s and P between liquid and gas
# Unknowns: rho_s_l, rho_s_g
# rho_i_l and rho_i_g are determined from mu_i_target after solution
# ------------------------
def system(vars, eps_si):
    rho_s_l, rho_s_g = vars
    # Compute rho_i in liquid and gas from mu_i_target
    rho_i_l = np.exp(mu_i_target - mu_i_ex(rho_s_l, eps_si))
    rho_i_g = np.exp(mu_i_target - mu_i_ex(rho_s_g, eps_si))
    
    return [
        mu_s(rho_s_l, rho_i_l, eps_si) - mu_s(rho_s_g, rho_i_g, eps_si),
        P(rho_s_l, rho_i_l, eps_si)    - P(rho_s_g, rho_i_g, eps_si)
    ]

# ------------------------
# Continuation in eps_si
# ------------------------
if __name__ == "__main__":
    # initial guess for single-component limit
    rho_s_l0, rho_s_g0 = 0.7, 0.01
    init_guess = [rho_s_l0, rho_s_g0]

    eps_steps = np.linspace(1e-10, 1.75, 100)  # continue to eps_si = 0.5
    results = []
    fl=open(f"../data/coex_mui_{mu_i_target}.dat","w")
    for eps_si in eps_steps:
        # Solve for solvent densities
        solution = fsolve(system, init_guess, args=(eps_si,))
        rho_s_l, rho_s_g = solution
        
        # Compute ion densities from semi-grand relation
        rho_i_l = np.exp(mu_i_target - mu_i_ex(rho_s_l, eps_si))
        rho_i_g = np.exp(mu_i_target - mu_i_ex(rho_s_g, eps_si))
        
        results.append([eps_si, rho_s_l, rho_s_g, rho_i_l, rho_i_g])
        
        # Use this solution as initial guess for next step
        init_guess = [rho_s_l, rho_s_g]
        
        print(f"eps_si={eps_si:.3f} | rho_s_l={rho_s_l:.4f}, rho_s_g={rho_s_g:.4f} | rho_i_l={rho_i_l:.4e}, rho_i_g={rho_i_g:.4e}")
        fl.write(f"{eps_si} {rho_s_l} {rho_s_g} {rho_i_l} {rho_i_g}\n")
    fl.close()

# After running, results contains the full eps_si sweep

