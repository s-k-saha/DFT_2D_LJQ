import numpy as np
from scipy.optimize import fsolve
import math

pi=np.pi
R=0.5
R2=R*R
R3=R*R2

eps_si =0.0 #solvent-ion interaction parameter
eps_ss =1.#solvent-solvent interaction parameter
alpha_si = 1.171861897 #integrated potential of solvent-ion
alpha_ss = 1.171861897 #integrated potential of solvent-solvent

rhoig=0.0

def mu_s(rhos,rhoi):
	eta=4./3*pi*R3*(rhos+2*rhoi)
	return np.log(rhos)+(14*eta-13*eta**2+5*eta**3)/(2*(1-eta)**3)-np.log(1-eta) + 4*pi*rhos*(-alpha_ss)*eps_ss +  8*pi*rhoi*(-alpha_si)*eps_si

def mu_i(rhos,rhoi):
	eta=4./3*pi*R3*(rhos+2*rhoi)
	return np.log(rhoi)+(14*eta-13*eta**2+5*eta**3)/(2*(1-eta)**3)-np.log(1-eta)  +  4*pi*rhos*(-alpha_si)*eps_si
	
def P(rhos,rhoi):
	eta=4./3*pi*R3*(rhos+2*rhoi)
	return (rhos+2*rhoi)*(1+eta+eta**2)/(1-eta)**3 + 8*pi*rhoi*rhos*(-alpha_si)*eps_si + 2*pi*rhos*rhos*(-alpha_ss)*eps_ss


def system(vars):
	rho_s_l,rho_s_g=vars
	rho_i_l=0.0
	rho_i_g=rhoig
	return [mu_s(rho_s_l,rho_i_l)-mu_s(rho_s_g,rho_i_g),P(rho_s_l,rho_i_l)-P(rho_s_g,rho_i_g)]
	
if __name__=="__main__":
	#eps_r=np.linspace(1.5,0.76,1000)
	rho_s_l,rho_s_g=0.6,0.01
	init_guess=[rho_s_l,rho_s_g]
	solution=fsolve(system,init_guess)
	print(f"{eps_ss} {eps_si} {solution[0]} {solution[1]} {mu_s(solution[0],0)}")
