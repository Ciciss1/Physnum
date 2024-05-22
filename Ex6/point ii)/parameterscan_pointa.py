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

V0s = [1e4, 1000]
outputs = ["pointa_V0-", "pointa_V0+"]
cases = ["V0<E", "V0>E"]

Ncases = len(V0s)

for i in range(Ncases):
    V0 = V0s[i]
    output = outputs[i]
    case = cases[i]

    # run the simulation

    cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


    #extract data
    data_obs = np.loadtxt(f"./outputs/{output}_obs.out")
    data_pot = np.loadtxt(f"./outputs/{output}_pot.out")
    data_psi2 = np.loadtxt(f"./outputs/{output}_psi2.out")

    t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T

    #extract avery 3 colomns
    psi_module = data_psi2[:,::3]

    print("mean E, ",case," : ",np.mean(E))

    #-----------------plot potential-----------------
    x, V = data_pot.T

    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$V$ [J]')

    ax.plot(x,V, label=r"$V(x)$")

    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18)
    u.savefig(fig, f"POINT_A_{case}_V(x)", ext = ext)


    #-----------------plot P(t) left and right-----------------
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'Probabilité')

    ax.plot(t,prob_left, label=r"Pr$_{x<x_c}(t)$")
    ax.plot(t,prob_right, label=r"Pr$_{x>x_c}(t)$")
    ax.plot(t,prob_left+prob_right, label=r"Pr$_{x<x_c}(t)$ + Pr$_{x>x_c}(t)$", linestyle="--", color = "black")

    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18)
    u.savefig(fig, f"POINT_A_P(t)_{case}", ext = ext)



    #---------------plot contours of psi2(x,t)-----------------

    #plot
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'$t$ [s]')
    c = ax.contourf(x, t, psi_module, levels=20, cmap = "plasma")
    fig.colorbar(c, ax=ax,label=r'$|\psi(x,t)|^2$')
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$t$ [s]')
    plt.tight_layout()
    u.savefig(fig, f"POINT_A_contour_{case}", ext = ext)
    
    #plot energy
    
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$t$ [s]", ylabel=r'$E$ [J]')
    
    ax.plot(t,E, label=r"$E(t)$")
    
    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18)
    u.savefig(fig, f"POINT_A_E(t)_{case}", ext = ext)
    
    # #---------------Animate-----------------
    
    # # Initialisation de la figure
    # fig, ax = plt.subplots()
    # line, = ax.plot([], [], lw=2)
    # line2, = ax.plot([], [], color="black", linestyle="--")

    # # Définition de la fonction d'initialisation de l'animation
    # def init():
    #     ax.set_xlim(-1.1, 1.1)
    #     ax.set_ylim(-1, 7)
    #     # ax.set_ylim(-1,1)
    #     return line,

    # # Fonction d'animation
    # def animate(i):
    #     line.set_data(x, psi_module[i])
    #     ax.set_title(f"Case : {cases[i]}")
    #     return line,

    # #------animate-------
    # # Création de l'animation
    # ani = FuncAnimation(fig, animate, frames=len(psi_module), init_func=init, blit=True, interval=40)

    # # Enregistrer l'animation au format GIF
    # # ani.save("simulation_vagues.gif", writer='pillow')

    # FFwriter=animation.FFMpegWriter(fps=2000/40, extra_args=['-vcodec', 'libx264'])
    # ani.save(f"./animations/simulation_{cases[i]}.mp4", writer=FFwriter)



