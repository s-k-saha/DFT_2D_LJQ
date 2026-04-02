import numpy as np
import matplotlib.pyplot as plt

fname="../data/newrho1Deps1.000000eps_i1.403153ew1.000000dx0.025000rhob_solvent0.603936rhob_solute0.019845lambdaB4.500000Vq1.000000.dat"	

data=np.loadtxt(fname,max_rows=1999)
x=data[:,0]
rho0=data[:,1]
rho1=data[:,2]
rho2=data[:,3]


rho0b=0.603936
rho1b=0.019845



plt.plot(x,rho0/rho0b,label="solvent",color='cornflowerblue')
plt.plot(x,rho1/rho1b,label="cation",color='red')
plt.plot(x,rho2/rho1b,label="anion",color='green')
plt.xlim(0.0,25)
plt.ylim(0.0,)
plt.ylabel(r"$\rho(z)/\rho_b$",size=20)
plt.xlabel(r"$z/\sigma$",size=20)
plt.legend(fontsize=15)
plt.tight_layout()
plt.show()
