import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import least_squares

#kB * T=1

eps_si =1. #solvent-ion interaction parameter
eps_ss =10.#solvent-solvent interaction parameter
alpha_si = 1. #integrated potential of solvent-ion
alpha_ss = 1. #integrated potential of solvent-solvent

Lambda3 = 1.# m = 2*pi*h_cut^2
R = .5
R3 = R*R*R#
pi=np.pi
rhoig=0.001

def f(x): #Rosenfeld FMT
	eta=4./3*pi*R3*x
	return  1./(pi*R3)*(-3./4 * eta * np.log(1-eta) + 9./4 * eta*eta/(1-eta) + 9./8 * eta*eta*eta/(1-eta)/(1-eta))
	  
def g(x): #delf/delx
	eta = 4./3*pi*R3*x
	return 4./3*(-3./4*np.log(1-eta) + 21./4*eta/(1-eta) + 45./8*eta*eta/(1-eta)/(1-eta)+9./4*eta*eta*eta/(1-eta)/(1-eta)/(1-eta))
	
def get_mu_i(rhos,rhoi):
	return np.log(Lambda3*rhoi)+g(rhos+2*rhoi)-eps_si*rhos*alpha_si
def get_mu_s(rhos,rhoi):
	return np.log(Lambda3*rhos)+g(rhos+2*rhoi)-eps_ss*rhos*alpha_ss-2*eps_si*rhoi*alpha_si
def get_P(rhos,rhoi):
	return ((2*rhoi+rhos) + (2*rhoi+rhos)*g(2*rhoi+rhos)-f(2*rhoi+rhos)- 2*eps_si*rhos*rhoi*alpha_si - .5* eps_ss * rhos*rhos*alpha_ss)

def solve_coexistence(get_mu_s, get_mu_i, get_P, rhos_guess, rhoi_guess):
    """
    Solve for liquid-gas coexistence densities for symmetric electrolyte.
    
    Parameters
    ----------
    get_mu_s : function
        Chemical potential of solvent: get_mu_s(rhos, rhoi)
    get_mu_i : function
        Chemical potential of ions: get_mu_i(rhos, rhoi)
    get_P : function
        Pressure: get_P(rhos, rhoi)
    rhos_guess : tuple
        Initial guesses for (rhos_liq, rhos_gas)
    rhoi_guess : tuple
        Initial guesses for (rhoi_liq, rhoi_gas)
    
    Returns
    -------
    rhos_liq, rhos_gas, rhoi_liq, rhoi_gas : floats
        Bulk coexistence densities
    """
    
    # Unknowns: [rhos_liq, rhos_gas, rhoi_liq, rhoi_gas]
    x0 = np.array([rhos_guess[0], rhos_guess[1], rhoi_guess[0], rhoig])
    
    def equations(x):
        rhos_l, rhos_g, rhoi_l, rhoi_g = x
        
        # Chemical potential differences
        mu_s_diff = get_mu_s(rhos_l, rhoi_l) - get_mu_s(rhos_g, rhoi_g)
        mu_i_diff = get_mu_i(rhos_l, rhoi_l) - get_mu_i(rhos_g, rhoi_g)
        
        # Pressure difference
        P_diff = get_P(rhos_l, rhoi_l) - get_P(rhos_g, rhoi_g)
        
        return [mu_s_diff, mu_i_diff, P_diff, rhos_l - rhos_g + 1e-3]  # last eq: small constraint to stabilize
        
    # Solve system
    sol, info, ier, mesg = fsolve(equations, x0, full_output=True)
    
    if ier != 1:
        raise RuntimeError(f"Coexistence solver did not converge: {mesg}")
    
    rhos_liq, rhos_gas, rhoi_liq, rhoi_gas = sol
    return rhos_liq, rhos_gas, rhoi_liq, rhoi_gas
    
    
    
if __name__=="__main__":
	rhos_guess=[0.8,.01]
	rhoi_guess=[.1,rhoig]
	rhos_liq, rhos_gas, rhoi_liq, rhoi_gas=solve_coexistence(get_mu_s, get_mu_i, get_P, rhos_guess, rhoi_guess)
	print(f"{rhos_liq} {rhos_gas} {rhoi_liq} {rhoi_gas}")
