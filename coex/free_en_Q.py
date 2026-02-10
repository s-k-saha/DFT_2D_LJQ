import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
import sys

pi=np.pi
R=0.5
R2=R*R
R3=R*R2

eps_si =float(sys.argv[1]) #solvent-ion interaction parameter
eps_ss =1.#solvent-solvent interaction parameter
alpha_si = 1.171861897 #integrated potential of solvent-ion
alpha_ss = 1.171861897 #integrated potential of solvent-solvent
tol=1.e-5


def f_FMT(x):
	n3=4./3*pi*R3*x
	n2=4.*pi*R2*x
	n1=R*x
	n0=x
	
	return -n0*np.log(1-n3)+n1*n2/(1-n3)+n2*n2*n2/(24*pi*(1-n3)*(1-n3))
	

def f(rhos,rhoi):
	return rhos*(np.log(rhos)-1)+ 2*rhoi*(np.log(rhoi)-1)+ f_FMT(rhos+2*rhoi) + eps_ss*rhos*rhos*(-alpha_ss)*2*pi + eps_si*rhos*rhoi*(-alpha_si)*8*pi
	

rhos_r=np.linspace(0.00001,1.2,500)
rhoi_r=np.linspace(0.00001,0.05,500)

with open(f"../data/coex_feng_eps_ss_{eps_ss}_eps_si_{eps_si}.dat","w") as fl:
	for rhos in rhos_r:
		for rhoi in rhoi_r:
			fl.write(f"{rhos} {rhoi} {f(rhos,rhoi)}\n")
		fl.write("\n");



