import sys
# Add the directory containing the module to sys.path
sys.path.append('/workspaces/Physnum')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import math as m
import utils_v2 as u

from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation
from matplotlib import animation

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

ext = "pdf"


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice6_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


#for input file
tfin = 0.1
xL = -1
xR = 1
n_v = 2
x0 = -0.5
n = 16
sigma_norm = 0.04

dt = 1e-4

V0 = 0
Nintervals = 256

output = "V0_0"
    
# # run the simulations
# cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')

x = np.linspace(xL, xR, Nintervals+1)

#extract data
data_obs = np.loadtxt(f"{output}_obs")
data_pot = np.loadtxt(f"{output}_pot")
data_psi2 = np.loadtxt(f"{output}_psi2")

t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T

psi_module2 = data_psi2[:,::3]
psi_module = np.sqrt(psi_module2)

psi_real = data_psi2[:,1::3]
psi_imag = data_psi2[:,2::3]

#-----------------plot xt-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r"$\langle x \rangle (t)$")

ax.plot(t,xmoy, label=r"$\langle x \rangle (t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"xmoyV0", ext = ext)

#-----------------plot pt-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r"$\langle p \rangle (t)$")

ax.plot(t,pmoy, label=r"$\langle p \rangle (t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"pmoyV0", ext = ext)

#-----------------plot dxt-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r"$\langle \Delta x \rangle (t)$")

ax.plot(t,xincertitude, label=r"$\langle \Delta x \rangle (t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"dxtV0", ext = ext)

#-----------------plot dpt-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r"$\langle \Delta p \rangle (t)$")

ax.plot(t,pincertitude, label=r"$\langle \Delta p \rangle (t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"dptV0", ext = ext)

#-----------------plot psi2-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
# levels1 = np.array([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
# levels2 = np.linspace(1,np.max(psi_module2),20)
# levels = np.concatenate((levels1,levels2))
c = ax.contourf(x, t, psi_module2, levels=20, cmap = "magma")
fig.colorbar(c, ax=ax,label=r'$|\psi(x,t)|^2$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, f"psi2V0", ext = ext)

#-----------------plot psi-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
# levels1 = np.array([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
# levels2 = np.linspace(1,np.max(psi_module),20)
# levels = np.concatenate((levels1,levels2))
c = ax.contourf(x, t, psi_module, levels=20, cmap = "magma")
fig.colorbar(c, ax=ax,label=r'$|\psi(x,t)|$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, f"psiV0", ext = ext)

#-----------------plot psi_real-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')

c = ax.contourf(x, t, psi_real, levels=20, cmap = "magma")
fig.colorbar(c, ax=ax,label=r'$Re(\psi(x,t))$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, f"psi_realV0", ext = ext)

#-----------------plot psi_imag-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
c = ax.contourf(x, t, psi_imag, levels=20, cmap = "magma")
fig.colorbar(c, ax=ax,label=r'$Im(\psi(x,t))$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, f"psi_imagV0", ext = ext)

#-----------------plot prob-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'Probabilité')


ax.plot(t,prob_left+prob_right, label=r"Pr$_{tot}(t)$", linestyle="--", color = "black")
plt.ylim(0, 1.1*np.max(prob_left+prob_right))
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"probV0", ext = ext)

#-----------------plot energy-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'$E$')
    
ax.plot(t,E, label=r"$E(t)$")
    
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"EV0", ext = ext)

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'$E$')
ax.plot(t,E, label=r"$E(t)$")
plt.ylim(0, 1.1*np.max(E))  # Set y-axis limit to show a larger scale
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"E_large_scaleV0", ext = ext)

#-----------------Eisenberg uncertainty principle-----------------

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'$\langle \Delta x \rangle (t) \langle \Delta p\rangle(t)$')

ax.plot(t,xincertitude*pincertitude, label=r"$\langle \Delta x \rangle (t) \langle \Delta p\rangle(t)$")
ax.plot(t, 0.5*np.ones_like(t), linestyle='--', color='black', label=r"$\hbar/2$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"HeinsebergV0", ext = ext)

v_0 = np.sqrt(2*E[0])

x_t = np.zeros_like(t)
x_t[0] = x0
for i in range(len(x_t)-1):
    x_t[i+1] = x_t[i] + v_0*dt
    if  x_t[i] >= xR and v_0 > 0:
        v_0 = -v_0
    elif x_t[i] <= xL and v_0 < 0:
        v_0 = -v_0

v_01 = pmoy[0]

x_t1 = np.zeros_like(t)
x_t1[0] = x0
for i in range(len(x_t1)-1):
    x_t1[i+1] = x_t1[i] + v_01*dt
    if  x_t1[i] >= xR and v_01 > 0:
        v_01 = -v_01
    elif x_t1[i] <= xL and v_01 < 0:
        v_01 = -v_01    

#-----------------plot x(t)-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'$x(t)$')

ax.plot(t,xmoy, label=r"$\langle x \rangle (t)$ quantique")
ax.plot(t, x_t, label=r'$x(t)$ classique même E(0)', linestyle='--')
ax.plot(t, x_t1, label=r'$x(t)$ classique même p(0)', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"x_t", ext = ext)

#-----------------plot psi_real and psi_module at different times-----------------
ax, fig = u.create_figure_and_apply_format(figsize=(8, 6), xlabel=r"$x$ [m]", ylabel=r'$Re(\psi)$ / $|\psi|$')
selected_times = [0, 0.015, 0.025]  # Choose the times you want to plot

for t_index in selected_times:
    t_index = int(t_index / dt)  # Convert time to index
    ax.plot(x, psi_real[t_index], label=r'$Re(\psi(x, t={})$'.format(np.round(t_index*dt,3)))
    ax.plot(x, psi_module[t_index], label=r'$|\psi(x, t={})$'.format(np.round(t_index*dt,3)))
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$Re(\psi)$ / $|\psi|$')

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18, ncol=2)
u.savefig(fig, f"psi_real_module_t", ext=ext)
