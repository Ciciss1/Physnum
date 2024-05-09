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

ext = "png"


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice5_students'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'input_example'  # Name of the input file


#for input file
cb_gauche = "sortie"
cb_droite = "sortie"
v_uniform = "false"

tfin = 14000

A = 1
x1 = 50e3
x2 = 250e3

equation_type="Eq2"
nx=1000
n_init=3
initialization="pas mode"
initial_state="right"

CFL=1.0
nsteps=40
impose_nsteps="false" 

output="./outputs/10000.out"
n_stride=20
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
xR= 1000e3

g = 9.81

#run the simulation

cmd = f"{repertoire}{executable} {input_filename} cb_gauche={cb_gauche} cb_droite={cb_droite} v_uniform={v_uniform} tfin={tfin} A={A} x1={x1} x2={x2} equation_type={equation_type} nx={nx} n_init={n_init} initialization={initialization} initial_state={initial_state} CFL={CFL} nsteps={nsteps} impose_nsteps={impose_nsteps} output={output} n_stride={n_stride} ecrire_f={ecrire_f} hL={hL} hR={hR} hC={hC} h00={h00} xa={xa} xb={xb} xc={xc} xd={xd} xL={xL} xR={xR} g={g}"

print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')



#extract
data = np.loadtxt(output + "_f")
t = data[:, 0]

N=nx+1

#-----------------plot la vague à différents t-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')

x = np.linspace(xL, xR, N)

for i in [0,8,20,30,60,110,150]:
    ax.plot(x,data[i,1:], label=f"t = {t[i]:.0f} s")
        
# u.add_ticks("x",[xa,xb,xc,xc],[r'$x_a$',r'$x_b$',r'$x_c$',r'$x_d$'],ax=ax,divide =5)
# ax.vlines((xc+xd)/2,-1,4,linestyles="--",color="black",label=r"$Pic du corail$")
# ax.vlines(xb,-1,4,linestyles="-.",color="red",label=r"$x_b$")
ax.vlines(xc,-1,4,linestyles="--",color="black",label=r"$x_c$")
ax.vlines(xd,-1,4,linestyles="-.",color="black",label=r"$x_d$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "Vague_difft", ext = ext)



#-----------------trouver amplitude et vitesse en fonction de x-----------------

t = data[:, 0]

t_crete = np.zeros(N)
A = np.zeros(N)
v = np.zeros(N)

#on itère sur les x
for i in range(N):
    f_i = data[:,i+1]
    t_crete_i_index = np.argmax(f_i)
    t_crete_i = t[t_crete_i_index]
        
    
    if t_crete_i_index == 0 or t_crete_i_index == 1:
        continue
    
    # print("i : ",i)
    # print(x[i])
    # print(t_crete_i_index)
    
    
    #on interpole quadratiquement autours du max de la crête
    a,b,c = np.polyfit(t[t_crete_i_index-2:t_crete_i_index+2],data[t_crete_i_index-2:t_crete_i_index+2,i+1],2)
    
    #l'amplitude de la crête est le max de la parabole
    t_crete_i = -b/(2*a)
    t_crete[i] = t_crete_i
    
    A_i = a*t_crete_i**2 + b*t_crete_i + c
    A[i] = A_i
    
    
    
#la vitesse est def comme v= (x_i+k - x_i-k)/(t_crete_i+k - t_crete_i-k) avec k un nombre entier choisi pour réduire les oscillations

# print(np.sum(t_crete==0))

ks = [5,10,15,20]

vs = np.zeros((len(ks),N))

j= 0
for k in ks:
    for i in range(k,N-k):
        v_i = (x[i+k] - x[i-k])/(t_crete[i+k] - t_crete[i-k])
        vs[j][i] = v_i
    j+=1

# print(vs)

#plot amplitude
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')

ax.plot(x,A, label=f"Amplitude", color="black")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "Amplitude", ext = ext)

 
#plot vitesse
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Vitesse [m/s]')

for i in range(len(ks)):
    ax.plot(x,vs[i], label=f"Vitesse k={ks[i]}")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "Vitesse", ext = ext)


#-----------------plot le profil h_0-----------------
dx = (xR - xL) / (N-1)
 

h0 = np.zeros(N)


for i in range(N):
    if x[i] >= 0 and x[i] <= xa:
        h0[i] = hL
    if x[i] > xa and x[i] < xb:
        h0[i] = 0.5*(hL + hC) + 0.5*(hL - hC)*m.cos(m.pi*((x[i]-xa)/(xb -xa)))
    if x[i] >= xb and x[i] <= xc:
        h0[i] = hC
    if x[i] > xc and x[i] < xd:
        h0[i] = 0.5*(hR + hC) - 0.5*(hR - hC)*m.cos(m.pi*((x[i]-xc)/(xd - xc)))
    if x[i] >= xd:
        h0[i] = hR
    
h0=-h0

#plot le profil de hauteur
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'h [m]')

ax.plot(x,h0, label=f"t = {t[0]:.2f} s", color="black")

# ax.vlines((xc+xd)/2,0,-7000,linestyles="--",color="black",label=r"$Pic du corail$")

ax.vlines(xc,0,-7000,linestyles="--",color="black",label=r"$x_c$")
ax.vlines(xd,0,-7000,linestyles="--",color="red",label=r"$x_d$")
ax.vlines(xa,0,-7000,linestyles="--",color="blue",label=r"$x_a$")
ax.vlines(xb,0,-7000,linestyles="--",color="green",label=r"$x_b$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "profil_h0", ext = ext)








# #animate
# # Initialisation de la figure
# fig, ax = plt.subplots()
# line, = ax.plot([], [], lw=2)
# line2, = ax.plot([], [], color="black", linestyle="--")

# # Définition de la fonction d'initialisation de l'animation
# def init():
#     ax.set_xlim(xL, xR)
#     ax.set_ylim(-4, 4)
#     # ax.set_ylim(-1,1)
#     return line,

# # Fonction d'animation
# def animate(i):
#     line.set_data(x, data[i, 1:])
#     line2.set_data(x, h0)
#     ax.set_title(f"t = {t[i]:.2f} s")
#     return line,

# #------animate-------
# # Création de l'animation
# ani = FuncAnimation(fig, animate, frames=len(data), init_func=init, blit=True, interval=40)

# # Enregistrer l'animation au format GIF
# ani.save("simulation_vagues.gif", writer='pillow')

# FFwriter=animation.FFMpegWriter(fps=1000/40, extra_args=['-vcodec', 'libx264'])
# ani.save("simulation_vagues.mp4", writer=FFwriter)




    

