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
Nintervals = 256

x = np.linspace(xL, xR, Nintervals+1)

#-----------------V0 >E-----------------
V0 = 100
output = "pointa_V0+"

# run the simulation

# cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')


#extract data
data_obs = np.loadtxt(f"./outputs/{output}_obs.out")
data_pot = np.loadtxt(f"./outputs/{output}_pot.out")
data_psi2 = np.loadtxt(f"./outputs/{output}_psi2.out")

t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T

#extract avery 3 colomns
psi_module = data_psi2[:,::3]

print("mean E, V0>E : ",np.mean(E))

#-----------------plot potential-----------------
x, V = data_pot.T

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$V$ [J]')

ax.plot(x,V, label=r"$V(x)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "POINT_A_V0>E_V(x)", ext = ext)


#-----------------plot P(t) left and right-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'Probabilité')

ax.plot(t,prob_left, label=r"Pr$_{x<x_c}(t)$")
ax.plot(t,prob_right, label=r"Pr$_{x>x_c}(t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "POINT_A_P(t)_V0>E", ext = ext)



#---------------plot contours of psi2(x,t)-----------------

#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
c = ax.contourf(x, t, psi_module, levels=10, cmap = "plasma")
fig.colorbar(c, ax=ax,label=r'$|\psi(x,t)|^2$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, "POINT_A_contour_V0>E", ext = ext)







#-----------------V0<E-----------------
V0 = 1e-5
output = "pointa_V0-"

# run the simulation

# cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')


#extract data
data_obs = np.loadtxt(f"./outputs/{output}_obs.out")
data_pot = np.loadtxt(f"./outputs/{output}_pot.out")
data_psi2 = np.loadtxt(f"./outputs/{output}_psi2.out")


t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T

#extract avery 3 colomns
psi_module = data_psi2[:,::3]

print("mean E, V0<E : ",np.mean(E))

#-----------------plot P(t) left and right-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'Probabilité')

ax.plot(t,prob_left, label=r"Pr$_{x<x_c}(t)$")
ax.plot(t,prob_right, label=r"Pr$_{x>x_c}(t)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "POINT_A_P(t)_V0<E", ext = ext)

#---------------plot contours of psi2(x,t)-----------------

#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
c = ax.contourf(x, t, psi_module, levels=10, cmap = "plasma")
fig.colorbar(c, ax=ax,label=r'$|\psi(x,t)|^2$')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$t$ [s]')
plt.tight_layout()
u.savefig(fig, "POINT_A_contour_V0<E", ext = ext)



