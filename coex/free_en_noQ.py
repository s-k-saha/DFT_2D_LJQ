import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
import sys

pi=np.pi
R=0.5
R2=R*R
R3=R*R2

eps_si =1e-10 #solvent-ion interaction parameter
eps_ss =float(sys.argv[1])#1.#solvent-solvent interaction parameter
alpha_si = 1.171861897 #integrated potential of solvent-ion
alpha_ss = 1.171861897 #integrated potential of solvent-solvent
tol=1.e-5


def f_FMT(x):
	n3=4./3*pi*R3*x
	n2=4.*pi*R2*x
	n1=R*x
	n0=x
	
	return -n0*np.log(1-n3)+n1*n2/(1-n3)+n2*n2*n2/(24*pi*(1-n3)*(1-n3))
	

def f(rhos):
	return rhos*(np.log(rhos)-1) + f_FMT(rhos) + eps_ss*rhos*rhos*(-alpha_ss)*2*pi
	

rho_r=np.linspace(0.00001,1.2,100000)
plt.plot(rho_r,f(rho_r))
plt.show()


