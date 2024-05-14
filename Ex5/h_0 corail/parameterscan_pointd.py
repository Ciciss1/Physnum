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
executable = 'Exercice5_students'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'input_example'  # Name of the input file


#for input file
cb_gauche = "sortie"
cb_droite = "sortie"
v_uniform = "false"

xas = [6e5,6.9e5]

tfin = 14000

A = 1
x1 = 50e3
x2 = 250e3

equation_type="Eq2"
nx=5000
n_init=3
initialization="pas mode"
initial_state="right"

CFL=1.0
nsteps=40
impose_nsteps="false" 

output="./outputs/5000_300000.0.out"
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


output_files = []
#create the output files 
for xa in xas:
    output_files.append(f"./outputs/5000_{xa}.out")


# #run the simulation
# for i in range(len(xas)):
#     cmd = f"{repertoire}{executable} {input_filename} cb_gauche={cb_gauche} cb_droite={cb_droite} v_uniform={v_uniform} tfin={tfin} A={A} x1={x1} x2={x2} equation_type={equation_type} nx={nx} n_init={n_init} initialization={initialization} initial_state={initial_state} CFL={CFL} nsteps={nsteps} impose_nsteps={impose_nsteps} output={output_files[i]} n_stride={n_stride} ecrire_f={ecrire_f} hL={hL} hR={hR} hC={hC} h00={h00} xa={xas[i]} xb={xb} xc={xc} xd={xd} xL={xL} xR={xR} g={g}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')

#extract
DATA = np.zeros((len(xas),919,nx+2))
for i in range(len(xas)):   
    data = np.loadtxt(output_files[i] + "_f")
    DATA[i]=data


t = data[:, 0]

N=nx+1

x = np.linspace(xL, xR, N)

for i in range(len(xas)):
    data = DATA[i]
    
    #-----------------plot la vague à différents t-----------------
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')

    for j in [75,150,275,600]:
        ax.plot(x,data[j,1:], label=f"t = {t[j]:.0f} s")
            
    # u.add_ticks("x",[xa,xb,xc,xc],[r'$x_a$',r'$x_b$',r'$x_c$',r'$x_d$'],ax=ax,divide =5)
    # ax.vlines((xc+xd)/2,-1,4,linestyles="--",color="black",label=r"$Pic du corail$")
    # ax.vlines(xb,-1,4,linestyles="-.",color="red",label=r"$x_b$")
    ax.vlines(xas[i],-1,4,linestyles="-.",color="red",label=r"$x_a$")
    ax.vlines(xb,-1,4,linestyles="--",color="blue",label=r"$x_b$")
    ax.vlines(xc,-1,4,linestyles="-.",color="purple",label=r"$x_c$")
    ax.vlines(xd,-1,4,linestyles="--",color="black",label=r"$x_d$")
    
    #plot a little black arrow from right to left
    ax.annotate("", xy=(0.08e6, 0.4), xytext=(0, 0.4), arrowprops=dict(arrowstyle="<-",lw=2))
    
    
    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18,ncol=2, loc = "upper left")
    u.savefig(fig, f"pointd_Vague_difft_{xas[i]}", ext = ext)



#-----------------trouver amplitude et vitesse en fonction de x-----------------

t = DATA[0,:, 0]

#-----------------profil h_0-----------------
dx = (xR - xL) / (N-1)

As = np.zeros((len(xas),N))
Vs = np.zeros((len(xas),N))

A0 =np.zeros(N)
v0 = np.zeros(N)






for j in range(len(xas)):
    print(xas[j])
    h0 = np.zeros(N)

    for k in range(N):
        if x[k] >= 0 and x[k] <= xas[j]:
            h0[k] = hL
        if x[k] > xas[j] and x[k] < xb:
            h0[k] = 0.5*(hL + hC) + 0.5*(hL - hC)*m.cos(m.pi*((x[k]-xas[j])/(xb -xas[j])))
        if x[k] >= xb and x[k] <= xc:
            h0[k] = hC
        if x[k] > xc and x[k] < xd:
            h0[k] = 0.5*(hR + hC) - 0.5*(hR - hC)*m.cos(m.pi*((x[k]-xc)/(xd - xc)))
        if x[k] >= xd:
            h0[k] = hR  
    h0=-h0

    u_WKB = np.sqrt(g*np.abs(h0))

    A0_1 =1*pow(g*hL,1/4) 
    A0_2 =1*pow(g*hL,-1/4)

    A_WKB_1 = A0_1*pow(g*np.abs(h0),-1/4)
    A_WKB_2 = A0_2*pow(g*np.abs(h0),1/4)
    
    
    
    u_corail = np.sqrt(g*np.abs(hC))
    u_apres_corail = np.sqrt(g*np.abs(hR))

    A_corail = A0_1*pow(g*np.abs(hC),-1/4)
    A_apres = A0_1*pow(g*np.abs(hR),-1/4)
    
    
    data = DATA[j]
    
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
        
        #on interpole quadratiquement autours du max de la crête
        a,b,c = np.polyfit(t[t_crete_i_index-2:t_crete_i_index+2],data[t_crete_i_index-2:t_crete_i_index+2,i+1],2)
        
        #l'amplitude de la crête est le max de la parabole
        t_crete_i = -b/(2*a)
        t_crete[i] = t_crete_i
        
        A_i = a*t_crete_i**2 + b*t_crete_i + c
        A[i] = A_i
        if j==0 : 
            # print("j : ",j, "i : ",i)
            A0[i] = A_i
    
    #la vitesse est def comme v= (x_i+k - x_i-k)/(t_crete_i+k - t_crete_i-k) avec k un nombre entier choisi pour réduire les oscillations
    k = 10
    for i in range(k,N-k):
        v_i = (x[i+k] - x[i-k])/(t_crete[i+k] - t_crete[i-k])
        v[i] = v_i    
        if j==0 : 
            v0[i] = v_i   
        
    
        
    #plot amplitude
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')
    
    ax.plot(x,A,label=fr"Amplitude simulée", color="black")

    ax.plot(x,A_WKB_1, label=f"Amplitude WKB éq.1", color="red",linestyle="--")
    # ax.plot(x,A_WKB_2, label=f"Amplitude WKB éq.2", color="blue",linestyle="--")
    # ax.plot(x,A0, label=fr"Amplitude $x_a$={xas[0]:.1e}", color="blue",linestyle="--")
    
    ax.vlines(xas[j],-1,4.2,linestyles="-.",color="orange",label=rf"$x_a = {xas[j]:.1e}$")
    ax.vlines(xb,-1,4.2,linestyles="-.",color="blue",label=r"$x_b$")
    ax.vlines(xc,-1,4.2,linestyles="-.",color="purple",label=r"$x_c$")
    ax.vlines(xd,-1,4.2,linestyles="-.",color="black",label=r"$x_d$")

    ax.set_xlim(0.185e6,1e6)
    ax.set_ylim(0.8,4.2)
    
    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18,ncol=1)
    u.savefig(fig, f"pointd_Amplitude_{xas[j]}", ext = ext)


    
    #plot vitesse
    ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$x$ [m]", ylabel=r'Vitesse [m/s]')

    ax.plot(x,v, label=fr"Vitesse simulée", color="black")
        
    ax.plot(x,u_WKB, label=f"Vitesse théorique", color="red",linestyle="--")

    ax.vlines(xas[j],-200,3500,linestyles="-.",color="orange",label=rf"$x_a = {xas[j]:.1e}$")
    ax.vlines(xb,-200,3500,linestyles="-.",color="blue",label=r"$x_b$")
    ax.vlines(xc,-200,3500,linestyles="-.",color="purple",label=r"$x_c$")
    ax.vlines(xd,-200,3500,linestyles="-.",color="black",label=r"$x_d$")

    ax.set_xlim(0.185e6,1e6)
    ax.set_ylim(-200,3500)

    plt.tight_layout()
    u.set_legend_properties(ax, fontsize=18,ncol=1,loc="upper left")
    u.savefig(fig, f"pointd_Vitesse_{xas[j]}", ext = ext)


    A_corail_simul = []
    u_corail_simul = []
    A_apres_corail_simul = []
    u_apres_corail_simul = []
    for i in range(N) : 
        if x[i] >= xb and x[i] <= xc:
            A_corail_simul.append(A[i])
            A_corail_simul.append(A[i])
            u_corail_simul.append(v[i])

        if x[i]>=xd : 
            A_apres_corail_simul.append(A[i])
            if v[i] != 0.0:
                u_apres_corail_simul.append(v[i])
            
            
    A_corail_simul = np.array(A_corail_simul)
    A_apres_corail_simul = np.array(A_apres_corail_simul)
    u_corail_simul = np.array(u_corail_simul)
    u_apres_corail_simul = np.array(u_apres_corail_simul)

        
        
    A_corail_simul_mean = np.mean(A_corail_simul)
    A_corail_simul_std = np.std(A_corail_simul)
    A_apres_corail_simul_mean = np.mean(A_apres_corail_simul)
    A_apres_corail_simul_std = np.std(A_apres_corail_simul)
    u_corail_simul_mean = np.mean(u_corail_simul)
    u_corail_simul_std = np.std(u_corail_simul)
    u_apres_corail_simul_mean = np.mean(u_apres_corail_simul)
    u_apres_corail_simul_std = np.std(u_apres_corail_simul)
        
    print("xa : ",xas[j])
    print("Amplitude sur le corail théorique : ",A_corail)
    print("Amplitude sur le corail simulée : ",A_corail_simul_mean,"+-",A_corail_simul_std)
    print("Erreur relative amplitude sur le corail simulée : ",(A_corail_simul_mean-A_corail)/A_corail*100,"%")

    print("Amplitude après le corail théorique : ",A_apres)
    print("Amplitude après le corail simulée : ",A_apres_corail_simul_mean,"+-",A_apres_corail_simul_std)
    print("Erreur relative amplitude après le corail simulée : ",(A_apres_corail_simul_mean-A_apres)/A_apres*100,"%")

    print("Vitesse sur le corail théorique : ",u_corail)
    print("Vitesse sur le corail simulée : ",u_corail_simul_mean,"+-",u_corail_simul_std)
    print("Errer relative vitesse sur le corail simulée : ",(u_corail_simul_mean-u_corail)/u_corail*100,"%")

    print("Vitesse après le corail théorique : ",u_apres_corail)
    print("Vitesse après le corail simulée : ",u_apres_corail_simul_mean,"+-",u_apres_corail_simul_std)
    print("Erreur relative vitesse après le corail simulée : ",(u_apres_corail_simul_mean-u_apres_corail)/u_apres_corail*100,"%")
        



# #--------------------------------------ANIMATE--------------------------------------

# for j in range(len(xas)):
#     # Initialisation de la figure
#     fig, ax = plt.subplots()
#     line, = ax.plot([], [], lw=2)
#     line2, = ax.plot([], [], color="black", linestyle="--")

#     # Définition de la fonction d'initialisation de l'animation
#     def init():
#         ax.set_xlim(xL, xR)
#         ax.set_ylim(-4, 4)
#         # ax.set_ylim(-1,1)
#         return line,

#     # Fonction d'animation
#     def animate(i):
#         line.set_data(x, DATA[j][i, 1:])
#         line2.set_data(x, h0)
#         ax.set_title(f"t = {t[i]:.2f} s")
#         return line,

#     #------animate-------
#     # Création de l'animation
#     ani = FuncAnimation(fig, animate, frames=len(DATA[j]), init_func=init, blit=True, interval=40)

#     # Enregistrer l'animation au format GIF
#     # ani.save("simulation_vagues.gif", writer='pillow')

#     FFwriter=animation.FFMpegWriter(fps=2000/40, extra_args=['-vcodec', 'libx264'])
#     ani.save(f"./animations/simulation_vagues_{xas[j]}.mp4", writer=FFwriter)
    




    

