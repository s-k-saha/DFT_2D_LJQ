import numpy as np
import matplotlib.pyplot as plt

mu=-8.0
eps_ll=1.0
eps_nl=0.75
ew=1.0
lambdaB=0.5
Vq=0.1
L=50.

rhob_solvent=0.611184
rhob_solute=0.000968
# --- Load data ---
fname=f"../data/newrho2DNeps{eps_ll:.6f}eps_i{eps_nl:.6f}ew{ew:.6f}dx0.025000rhob_solvent_l{rhob_solvent:.6f}rhob_solute_l{rhob_solute:.6f}lambdaB{lambdaB:.6f}Vq{Vq:.6f}.dat"

print(fname)
data = np.loadtxt(
    fname,max_rows=1999*1999
)

# --- Extract ---
x = data[:, 0]
y = data[:, 1]
rhol = data[:, 2]
rhon1 = data[:, 3]
rhon2 = data[:, 4]


# --- Get grid size ---
Nx = len(np.unique(x))
Ny = len(np.unique(y))

# --- Reshape ---
X = x.reshape(Nx, Ny)
Y = y.reshape(Nx, Ny)
rhol_grid = rhol.reshape(Nx, Ny)
rhon1_grid = rhon1.reshape(Nx, Ny)
rhon2_grid = rhon2.reshape(Nx, Ny)

dx = 0.025
step = 200

xticks = np.arange(x.min(), x.max(), step)
yticks = np.arange(y.min(), y.max(), step)

# --- Plot rhol ---
fig1, ax1 = plt.subplots()
im1 = ax1.imshow(
    rhol_grid.T,
    extent=(x.min(), x.max(), y.min(), y.max()),
    origin='lower',
    cmap='Blues',
    aspect='auto'
)
fig1.colorbar(im1, ax=ax1, label=r'$\rho_l$')
ax1.set_title(rf'$\rho_l(x,z); \epsilon_w ={ew}; \lambda_B = {lambdaB}; V_q = {Vq}$')
ax1.set_xlabel('x')
ax1.set_ylabel('z')

ax1.set_xticks(xticks)
ax1.set_yticks(yticks)
ax1.set_xticklabels(xticks * dx - L/2)
ax1.set_yticklabels(yticks * dx)

plt.tight_layout()
plt.savefig(f'../figs/2DN_rhol_ew{ew:.6f}eps{eps_ll:.6f}eps_i{eps_nl:.6f}_rhob_solvent_l{rhob_solvent:.6f}rhob_solute_l{rhob_solute:.6f}lambdaB{lambdaB:.6f}Vq{Vq:.6f}_dx0.025000Lx49.950000.png')


vmin1=min(rhon1)
vmin2=min(rhon2)
vmin=min([vmin1,vmin2])

vmax1=max(rhon1)
vmax2=max(rhon2)
vmax=max([vmax1,vmax2])


# --- Plot rhon1 ---
fig2, ax2 = plt.subplots()
im2 = ax2.imshow(
    rhon1_grid.T,
    extent=(x.min(), x.max(), y.min(), y.max()),
    origin='lower',
    cmap='Greys',
    vmin=vmin,
    vmax=vmax,
    aspect='auto'
)
fig2.colorbar(im2, ax=ax2, label=r'$\rho_n1$')
ax2.set_title(rf'$\rho_{{n1}}(x,z); \epsilon_w ={ew}; \lambda_B = {lambdaB}; V_q = {Vq}$')
ax2.set_xlabel('x')
ax2.set_ylabel('z')

ax2.set_xticks(xticks)
ax2.set_yticks(yticks)
ax2.set_xticklabels(xticks * dx - L/2)
ax2.set_yticklabels(yticks * dx)

plt.tight_layout()

plt.savefig(f'../figs/2DN_rhon1_ew{ew:.6f}eps{eps_ll:.6f}eps_i{eps_nl:.6f}_rhob_solvent_l{rhob_solvent:.6f}rhob_solute_l{rhob_solute:.6f}lambdaB{lambdaB:.6f}Vq{Vq:.6f}_dx0.025000Lx49.950000.png')

# --- Plot rhon2 ---
fig3, ax3 = plt.subplots()
im3 = ax3.imshow(
    rhon2_grid.T,
    extent=(x.min(), x.max(), y.min(), y.max()),
    origin='lower',
    cmap='Greys',
    vmin=vmin,
    vmax=vmax,
    aspect='auto'
)
fig2.colorbar(im3, ax=ax3, label=r'$\rho_n2$')
ax3.set_title(rf'$\rho_{{n2}}(x,z); \epsilon_w ={ew}; \lambda_B = {lambdaB}; V_q = {Vq}$')
ax3.set_xlabel('x')
ax3.set_ylabel('z')

ax3.set_xticks(xticks)
ax3.set_yticks(yticks)
ax3.set_xticklabels(xticks * dx - L/2)
ax3.set_yticklabels(yticks * dx)

plt.tight_layout()


plt.savefig(f'../figs/2DN_rhon2_ew{ew:.6f}eps{eps_ll:.6f}eps_i{eps_nl:.6f}_rhob_solvent_l{rhob_solvent:.6f}rhob_solute_l{rhob_solute:.6f}lambdaB{lambdaB:.6f}Vq{Vq:.6f}_dx0.025000Lx49.950000.png')


plt.show()

