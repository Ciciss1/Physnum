import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import utils_v2 as u

from scipy.signal import find_peaks


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Ex2_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file
ext = "pdf"

nsteps = np.array([600,650,700,750,800,860,900,950,1000,1100,1250,1400,1450,1500,1750,2000,2250,2500,2750,3000,4000,5000,10000,15000,20000]) # TODO change as you need
nsimul = len(nsteps)  # Number of simulations to perform

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need


theta_0 = 0.5
# Analysis
# TODO: Insert here the expressions for the exact final solution
final_theta_th = 0.5
final_thetadot_th = 0

L = 0.2
g = 9.81
alpha = -0.02
T = 0.96*L/abs(alpha)


paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  

for i in range(nsimul):
    output_file = f"./outputs/{paramstr}={param[i]}.out"
    outputs.append(output_file)
    
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


final_theta_list = []
final_thetadot_list = []

thetasimul_list = []
thetadotsimul_list = []
Esimul_list = []

for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t  = data[:, 0]
    final_theta = data[-1, 1] 
    final_thetadot = data[-1, 2]  
    final_theta_list.append(final_theta)
    final_thetadot_list.append(final_thetadot)

    
    if i == nsimul-1 : 
        thetasimul = data[:, 1]
        thetadotsimul = data[:, 2]
        Esimul = data[:, 3]
    
    thetasimul_list.append(data[:, 1])
    thetadotsimul_list.append(data[:, 2])
    Esimul_list.append(data[:, 3])
    
    
    
    
#-----------plotting the error in theta------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\Delta t^2$ [s]", ylabel=r'$\theta$ [rad]', grid_bool=False)
ax.plot(dt**2, final_theta_list,linestyle = "-", marker = "+", color = "blue", markersize = 15) # this is an example, change it if you need

#-----------plotting the expected order of convergence------------
x = np.linspace(np.min(dt),np.max(dt),1000)
# ax.plot(x, 0.00001*x**4, color = "red", label = r"~ $\Delta_t^4$", linestyle = "--")

#in case we need linear fit
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_theta), color = "red", ax = ax, precisions = [4,4], plot = False)
# ax.plot(x, (x**a)*np.exp(b), color = "red") 
# print(f"a = {a:.3f} et b = {np.exp(b):.7f}")

#-----------formatting------------
ax.legend()    
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "RectracteFil_convergence",ext)


#-----------plot positions------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]', grid_bool=False)
#theta_simul!!!!!!
t = np.linspace(0,T,len(thetasimul))
l = alpha*t + L
y = -l*np.cos(thetasimul)
x = l*np.sin(thetasimul)
ax.axis('equal')
ax.plot(x,y,linestyle = "-", color = "navy") # this is an example, change it if you need
# plt.savefig(f"png/positions.png")
plt.tight_layout()
ax.grid()
u.savefig(fig, "RectracteFil_xy",ext)


#-----------create a function to compute dE/dt------------
def compute_derivative(E, dt):
    n = len(E)
    
    P = np.zeros(n)
    
    for i in range(1, n - 1):
        P[i] = (E[i+1] - E[i-1]) / (2 * dt)
    
    P[0] = (E[1] - E[0]) / dt
    P[n - 1] = (E[n - 1] - E[n - 2]) / dt
    
    return P


m = 0.1
t = np.linspace(0,T,len(thetasimul))

def length(t):
    return L + alpha*t

def lendot(t):
    return alpha

def lendotdot(t):
    return 0

#-----------compute theorical Pnonc with parameters from the simulation------------
thetadotdot_ = (-g*np.sin(thetasimul) - 2 * lendot(t) * thetadotsimul)/length(t)
P_th = m*(lendot(t) * lendotdot(t) + length(t) * lendot(t) * pow(thetadotsimul,2) + thetadotsimul * thetadotdot_ * pow(length(t),2) - g * lendot(t) * np.cos(thetasimul) + g * length(t) * thetadotsimul * np.sin(thetasimul))

#-----------compute dE/dt------------
DT = t[1]-t[0]
P_simul = compute_derivative(Esimul, DT)

#-----------plot dE/dt------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\frac{dE}{dt}$ [W]', grid_bool=False)

perc = 0.98
for i in [2,20,23] : 
    # print("i = ",i)
    t = np.linspace(T*perc,T,int(len(Esimul_list[i])*(1-perc))+1)
    DT = t[1]-t[0]
    # print(t)
    # print(np.shape(t),np.shape(Esimul_list[i]),np.shape(compute_derivative(Esimul_list[i],DT)))
    dE_dt = compute_derivative(Esimul_list[i][int(len(Esimul_list[i])*perc):], DT)
    if i ==23 :
        ax.plot(t,dE_dt,linestyle = "-",label = r"$N_{steps}$" + f" = {nsteps[i]}",color = "blue") 
    else : 
        ax.plot(t,dE_dt,linestyle = "-",label = r"$N_{steps}$" + f" = {nsteps[i]}") 

t = np.linspace(T*perc,T,int(len(thetasimul)*(1-perc))+1)
ax.plot(t,P_th[int(len(P_th)*perc):],linestyle = "--", color = "black", label = r"$P_{nonc}$ théorique") 

plt.tight_layout()
ax.grid()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
u.savefig(fig, "RectracteFil_Pnonc",ext)


#-----------plot the difference between the theorical and the simulated Pnonc------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 8),ylabel=r"$\delta_P$ [W]", xlabel=r'$t$ [s]', grid_bool=False)

perc = 0.9
for i in [20,21,22,23] : 
    # print("i = ",i)
    t_ = np.linspace(T*perc,T,int(len(Esimul_list[i])*(1-perc))+1)
    DT = t_[1]-t_[0]
    # print(t)
    # print(np.shape(t),np.shape(Esimul_list[i]),np.shape(compute_derivative(Esimul_list[i],DT)))
    dE_dt = compute_derivative(Esimul_list[i][int(len(Esimul_list[i])*perc):], DT)
    
    t = np.linspace(0,T,len(thetasimul_list[i]))
    thetadotdot_ = (-g*np.sin(thetasimul_list[i]) - 2 * lendot(t) * thetadotsimul_list[i])/length(t)
    
    P_th = m*(lendot(t) * lendotdot(t) + length(t) * lendot(t) * pow(thetadotsimul_list[i],2) + thetadotsimul_list[i] * thetadotdot_ * pow(length(t),2) - g * lendot(t) * np.cos(thetasimul_list[i]) + g * length(t) * thetadotsimul_list[i] * np.sin(thetasimul_list[i]))
    
    P = P_th[int(len(P_th)*perc):] 
    ax.plot(t_,dE_dt-P,linestyle = "-",label = r"$N_{steps}$" + f" = {nsteps[i]}") 

# t = np.linspace(0,T,len(thetasimul))
# DT = t[1]-t[0]
# dE_dt = compute_derivative(Esimul, DT)

# thetadotdot_ = (-g*np.sin(thetasimul) - 2 * lendot(t) * thetadotsimul)/length(t)
# P_th = m*(lendot(t) * lendotdot(t) + length(t) * lendot(t) * pow(thetadotsimul,2) + thetadotsimul * thetadotdot_ * pow(length(t),2) - g * lendot(t) * np.cos(thetasimul) + g * length(t) * thetadotsimul * np.sin(thetasimul))
# P = P_th
# ax.plot(t,dE_dt-P,linestyle = "-",color = "navy",label = r"$\delta_P$")


plt.tight_layout()
ax.grid()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
u.savefig(fig, "RectracteFil_deltaP",ext)



#-----------plot theta(t)----------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\theta$ [rad]')

t = np.linspace(0,T,len(thetasimul))
ax.plot(t,thetasimul,linestyle = "-", color = "navy", label = "Numerical solution")

#plot h line at 2pi and -2pi
ax.axhline(y=pi, color='r', linestyle='--', label = r"$y = \pi$")
ax.axhline(y=-pi, color='r', linestyle='--', label = r"$y = -\pi$")

u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)

plt.tight_layout()
u.savefig(fig, "RectracteFil_theta(t)",ext)


#-----------plot thetadot(t)----------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\dot{\theta}$ [rad/s]')

t = np.linspace(0,T,len(thetadotsimul))
ax.plot(t,thetadotsimul,linestyle = "-", color = "navy")

plt.tight_layout()
u.savefig(fig, "RectracteFil_thetadot(t)",ext)

#-----------plot E(t)----------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$E$ [J]')

t = np.linspace(0,T,len(Esimul))
ax.plot(t,Esimul,linestyle = "-", color = "navy")

plt.tight_layout()
u.savefig(fig, "RectracteFil_E(t)",ext)



