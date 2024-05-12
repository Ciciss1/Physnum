import sys
# Add the directory containing the module to sys.path
sys.path.append('/workspaces/Physnum')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import math as m
import utils_v2 as u

from matplotlib.animation import FuncAnimation

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

ext = "pdf"


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice5_students'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'input_example'  # Name of the input file



#for input file
cb_gauche = "fixe"
cb_droite = "libre"
v_uniform = "true"


def k_n(n,xL,xR) : 
    return np.pi*(n+0.5)/(xR-xL)

def modeanalytique(x,t,n,g,h00,xL,xR,A) : 
    k = k_n(n,xL,xR)
    omega = k*np.sqrt(g*h00)
    return A*np.sin(k*(x-xL))*np.cos(-omega*t)

A = 1
x1 = 2
x2 = 6

equation_type="Eq1"
nx=20
n_init=1
initialization="mode"
initial_state="static"

CFL=1.0
nsteps=30
impose_nsteps="true" 

output="./outputs/erreurmodepropre.out"
n_stride=1
ecrire_f=1

hL = 7000.0
hR = 200.0
hC = 35.0
h00= 3.0
xa= 3e5
xb= 7e5
xc= 7.2e5
xd= 8.5e5
xL= 0.0 
xR= 10

g = 9.81

# frq = k_n(n_init,xL,xR)*np.sqrt(g*h00)/(2*np.pi)

# tfin = 1/frq
# print("tfin : ", tfin)

# dt = tfin/nsteps
# dx = (xR-xL)/nx

# u_ = np.sqrt(g*h00)
# beta_CFL = u_ * dt/dx

# print("betaCFL : ", beta_CFL)





# #run the simulation

# cmd = f"{repertoire}{executable} {input_filename} cb_gauche={cb_gauche} cb_droite={cb_droite} v_uniform={v_uniform} tfin={tfin} A={A} x1={x1} x2={x2} equation_type={equation_type} nx={nx} n_init={n_init} initialization={initialization} initial_state={initial_state} CFL={CFL} nsteps={nsteps} impose_nsteps={impose_nsteps} output={output} n_stride={n_stride} ecrire_f={ecrire_f} hL={hL} hR={hR} hC={hC} h00={h00} xa={xa} xb={xb} xc={xc} xd={xd} xL={xL} xR={xR} g={g}"

# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')



# #extract
# data = np.loadtxt(output + "_f")
# t = data[:, 0]

# N=nx+1

# #plot
# ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')

# x = np.linspace(xL, xR, N)

# mode = modeanalytique(x, tfin, n_init, g, h00, xL, xR, A)

# erreur = np.trapz(np.abs(data[-1, 1:] - mode), x)
# print("erreur =", erreur)

# ax.plot(x, mode, label=f"Mode analytique t = {tfin:.2f} s", linestyle='--')

# for i in [100]:
#     ax.plot(x,data[i,1:], label=f"t = {t[i]:.2f} s")
    

# Add arrows to show wave propagation direction

# ax.annotate("", xy=(3.2, 0.2), xytext=(4.2, 0.2),
#                 arrowprops=dict(arrowstyle="->", color="black"))
# ax.annotate("", xy=(8.5, -0.2), xytext=(7.5, -0.2),
#                 arrowprops=dict(arrowstyle="->", color="black"))
# ax.annotate("", xy=(8, -0.8), xytext=(9, -0.8),
#                 arrowprops=dict(arrowstyle="->", color="black"))
# ax.annotate("", xy=(4, -0.4), xytext=(3, -0.4),
#                 arrowprops=dict(arrowstyle="->", color="black"))

# plt.tight_layout()
# u.set_legend_properties(ax, fontsize=18,ncol=2)
# u.savefig(fig, "erreurmodepropre", ext = ext)



# #animate
# # Initialisation de la figure
# fig, ax = plt.subplots()
# line, = ax.plot([], [], lw=2)


# # Définition de la fonction d'initialisation de l'animation
# def init():
#     ax.set_xlim(xL, xR)
#     ax.set_ylim(-4, 3)
#     return line,

# # Fonction d'animation
# def animate(i):
#     line.set_data(x, data[i, 1:])
#     ax.set_title(f"t = {t[i]:.2f} s")

#     ax.hlines(-3, xL, xR, color='brown', linestyle='--')
#     return line,

# # Création de l'animation
# ani = FuncAnimation(fig, animate, frames=len(data), init_func=init, blit=True, interval=40)

# # Enregistrer l'animation au format GIF
# ani.save("simulation_vagues.gif", writer='pillow')

# # Affichage de l'animation
# plt.show()


# Initialize lists to store the nsteps values and corresponding errors
dt_list = []
dx_list = []
error_list = []

# Loop over the nsteps values
for i in range(10):

    frq = k_n(n_init,xL,xR)*np.sqrt(g*h00)/(2*np.pi)

    tfin = 1/frq
    print("tfin : ", tfin)

    dt = tfin/nsteps
    dx = (xR-xL)/nx

    u_ = np.sqrt(g*h00)
    beta_CFL = u_ * dt/dx

    print("betaCFL : ", beta_CFL)

    cmd = f"{repertoire}{executable} {input_filename} cb_gauche={cb_gauche} cb_droite={cb_droite} v_uniform={v_uniform} tfin={tfin} A={A} x1={x1} x2={x2} equation_type={equation_type} nx={nx} n_init={n_init} initialization={initialization} initial_state={initial_state} CFL={CFL} nsteps={nsteps} impose_nsteps={impose_nsteps} output={output} n_stride={n_stride} ecrire_f={ecrire_f} hL={hL} hR={hR} hC={hC} h00={h00} xa={xa} xb={xb} xc={xc} xd={xd} xL={xL} xR={xR} g={g}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')
       
    # Load the data
    data = np.loadtxt(output + "_f")
       
    # Calculate the mode analytique
    N=nx+1
    x = np.linspace(xL, xR, N)
    mode = modeanalytique(x, tfin, n_init, g, h00, xL, xR, A)
      
    # Calculate the error
    erreur = np.trapz(np.abs(data[-1, 1:] - mode), x)
        
    # Append the nsteps value and error to the lists
    dt_list.append(dt)
    dx_list.append(dx)
    error_list.append(erreur)

    nx = 2*nx
    nsteps = 2*nsteps

# # Plot the error vs nsteps
# plt.plot(1/(np.array(nsteps_list)**3), error_list)
# plt.xlabel("1/nsteps^2")
# plt.ylabel("Error")
# plt.title("Error vs 1/nsteps^2")
# plt.show()

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\Delta t$", ylabel=r'Erreur position finale')

ax.loglog(np.array(dt_list), error_list,'-x', label=f"Erreur position finale")

x_line = np.linspace(0, 0.1, 100)
y_line = 5*x_line
plt.loglog(x_line, y_line, linestyle='--', color='red', label=r"$\sim \Delta t$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "erreurmodepropredt", ext = ext)

ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\Delta x$", ylabel=r'Erreur position finale')

ax.loglog(np.array(dx_list), error_list,'-x', label=f"Erreur position finale")

x_line = np.linspace(0, 1, 100)
y_line = 0.5*x_line
plt.loglog(x_line, y_line, linestyle='--', color='red', label=r"$\sim \Delta x$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "erreurmodepropredx", ext = ext)