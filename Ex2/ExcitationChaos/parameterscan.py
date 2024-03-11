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

nsteps = np.array([500,550,600,650,700,750,800,850,900,950,1000,1100,1250,1500,1750,2000,2250,2500,2750,3000,4000,5000,10000,15000,20000]) # TODO change as you need
nsimul = len(nsteps)  # Number of simulations to perform

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need


theta_0 = 0.01
thetadot_0 = 0.0

L = 0.2
g = 9.81
alpha = 0
d = 0.01

omega_0 = np.sqrt(g/L)
omega = 2*omega_0



T_excitation = 2*pi/omega

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan


#------------------------------------------------
#-------------N=10-------------------------------
#------------------------------------------------

# Simulations
outputs = []  

#periode d excitation
N = 10
tfin = N*T_excitation

for i in range(nsimul):
    output_file = f"./outputs/{paramstr}={param[i]}_N={N}.out"
    outputs.append(output_file)
    
# for i in range(nsimul):
#     output_file = outputs[i]
#     cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file} N={N}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')


final_theta_list = []

for i in range(nsimul) : 
    data = np.loadtxt(outputs[i])
    
    final_theta = data[-1,1]
    final_theta_list.append(final_theta)

data = np.loadtxt(outputs[-1])
thetasimul = data[:, 1]
thetadotsimul = data[:, 2]
Esimul = data[:, 3]

    
# -----------plotting the error in theta------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta t^2$ [s$^2$]", ylabel=r'$\theta_{fin}$ [rad]', grid_bool=False)
ax.plot(dt**2, final_theta_list,linestyle = "-", marker = "+", color = "blue", markersize = 15) # this is an example, change it if you need

# # -----------plotting the expected order of convergence------------
# x = np.linspace(np.min(dt),np.max(dt),1000)
# ax.plot(x, 0.00001*x**4, color = "red", label = r"~ $\Delta_t^4$", linestyle = "--")

# in case we need linear fit
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_theta), color = "red", ax = ax, precisions = [4,4], plot = False)
# ax.plot(x, (x**a)*np.exp(b), color = "red") 
# print(f"a = {a:.3f} et b = {np.exp(b):.7f}")

# -----------formatting------------
plt.tight_layout()
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "ExcitationChaos_convergence_N10",ext)


#-----------plot positions------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]', grid_bool=False)
t = np.linspace(0,tfin,len(thetasimul))
l = alpha*t + L + d*np.sin(omega*t)
y = -l*np.cos(thetasimul)
x = l*np.sin(thetasimul)
ax.plot(x,y,linestyle = "-", color = "navy") # this is an example, change it if you need
# plt.savefig(f"png/positions.png")
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_xy_N10",ext)


#----------plot phase diagramm------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta$ [rad]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)

ax.plot(thetasimul,thetadotsimul,linestyle = "-", color = "navy") # this is an example, change it if you need
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_PortraitPhase_N10",ext)

#-----------plot theta(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\theta$ [rad]', grid_bool=False)

t = np.linspace(0,tfin,len(thetasimul))
ax.plot(t,thetasimul,linestyle = "-", color = "navy") 
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_theta(t)_N10",ext)


#-----------plot thetadot(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)

t = np.linspace(0,tfin,len(thetadotsimul))
ax.plot(t,thetadotsimul,linestyle = "-", color = "navy") 
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_thetadot(t)_N10",ext)

#-----------plot E(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$E$ [J]', grid_bool=False)

t = np.linspace(0,tfin,len(Esimul))
ax.plot(t,Esimul,linestyle = "-", color = "navy")
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_E(t)_N10",ext)







#------------------------------------------------
#-------------N=20-------------------------------
#------------------------------------------------

# Simulations
outputs = []  

#periode d excitation
N = 20
tfin = N*T_excitation

for i in range(nsimul):
    output_file = f"./outputs/{paramstr}={param[i]}_N={N}.out"
    outputs.append(output_file)
    
# for i in range(nsimul):
#     output_file = outputs[i]
#     cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file} N={N}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')


final_theta_list = []

for i in range(nsimul) : 
    data = np.loadtxt(outputs[i])
    
    final_theta = data[-1,1]
    final_theta_list.append(final_theta)

data = np.loadtxt(outputs[-1])
thetasimul = data[:, 1]
thetadotsimul = data[:, 2]
Esimul = data[:, 3]

    
# -----------plotting the error in theta------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta t$ [s]", ylabel=r'$\theta_{fin}$ [rad]')
ax.plot(dt, final_theta_list,linestyle = "-", marker = "+", color = "blue", markersize = 15) # this is an example, change it if you need

# # -----------plotting the expected order of convergence------------
# x = np.linspace(np.min(dt),np.max(dt),1000)
# ax.plot(x, 0.00001*x**4, color = "red", label = r"~ $\Delta_t^4$", linestyle = "--")

# in case we need linear fit
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_theta), color = "red", ax = ax, precisions = [4,4], plot = False)
# ax.plot(x, (x**a)*np.exp(b), color = "red") 
# print(f"a = {a:.3f} et b = {np.exp(b):.7f}")

# -----------formatting------------
xticks = [0,2.5e-4,5e-4,7.5e-4,1e-3,1.25e-3,1.5e-3,1.75e-3,2e-3]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{i*1e3}" for i in xticks])
#display a tiny 1e3 on the x axis
ax.text(1.0, -0.1, r"$\times 10^{-3}$", fontsize=12, ha='center', va='center', transform=ax.transAxes)
plt.tight_layout()
u.savefig(fig, "ExcitationChaos_convergence_N20",ext)


#-----------plot positions------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]', grid_bool=False)
t = np.linspace(0,tfin,len(thetasimul))
l = alpha*t + L + d*np.sin(omega*t)
y = -l*np.cos(thetasimul)
x = l*np.sin(thetasimul)
ax.plot(x,y,linestyle = "-", color = "navy") # this is an example, change it if you need
# plt.savefig(f"png/positions.png")
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_xy_N20",ext)


#----------plot phase diagramm------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta$ [rad]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)

ax.plot(thetasimul,thetadotsimul,linestyle = "-", color = "navy") # this is an example, change it if you need
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_PortraitPhase_N20",ext)

#-----------plot theta(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\theta$ [rad]', grid_bool=False)

t = np.linspace(0,tfin,len(thetasimul))
ax.plot(t,thetasimul,linestyle = "-", color = "navy") 
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_theta(t)_N20",ext)


#-----------plot thetadot(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)

t = np.linspace(0,tfin,len(thetadotsimul))
ax.plot(t,thetadotsimul,linestyle = "-", color = "navy") 
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_thetadot(t)_N20",ext)

#-----------plot E(t)------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$E$ [J]', grid_bool=False)

t = np.linspace(0,tfin,len(Esimul))
ax.plot(t,Esimul,linestyle = "-", color = "navy")
plt.tight_layout()
ax.grid()
u.savefig(fig, "ExcitationChaos_E(t)_N20",ext)







