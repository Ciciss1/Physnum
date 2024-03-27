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

nsteps = np.array([20,50,100,200,500,1000,2000,5000]) # TODO change as you need
nsimul = len(nsteps)  # Number of simulations to perform

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need


theta_0 = 1e-10
# Analysis
# TODO: Insert here the expressions for the exact final solution
final_theta_th = 1e-10
final_thetadot_th = 0

l = 0.2
g = 9.81
omega_0 = np.sqrt(g/l)
T = 2*np.pi/omega_0


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

error_theta = np.zeros(nsimul)
error_thetadot = np.zeros(nsimul)

final_theta_list = []
final_thetadot_list = []

for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t  = data[:, 0]
    final_theta = data[-1, 1] 
    final_thetadot = data[-1, 2]  
    final_theta_list.append(final_theta)
    final_thetadot_list.append(final_thetadot)
    # compute the error if you can ..
    error_theta[i] = abs(final_theta_th - final_theta)
    error_thetadot[i] = abs(final_thetadot_th - final_thetadot)
    
    
    
    
    
#-----------plotting the error in theta------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_\theta$ [rad]', grid_bool=False)
ax.loglog(dt, error_theta,linestyle = "-", marker = "+", label = r"Erreur sur $\theta$", color = "blue", markersize = 15) # this is an example, change it if you need

#-----------plotting the expected order of convergence------------
x = np.linspace(np.min(dt),np.max(dt),1000)
ax.plot(x, 0.00001*x**4, color = "red", label = r"~ $\Delta_t^4$", linestyle = "--")

#in case we need linear fit
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_theta), color = "red", ax = ax, precisions = [4,4], plot = False)
# ax.plot(x, (x**a)*np.exp(b), color = "red") 
# print(f"a = {a:.3f} et b = {np.exp(b):.7f}")

# #-----------formatting------------
# ax.legend()    
# u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
# plt.grid(True, which="both", linestyle='--')
# u.save_png(fig, "PasExcitation_erreur_theta" + "." + ext)


#-----------plotting the error in thetadot------------
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_\dot{\theta}$ [m]', grid_bool=False)
ax.loglog(dt, error_thetadot,linestyle = "-", marker = "+", label = r"Erreur sur $\dot{\theta}$", markersize = 15, color = "black") # this is an example, change it if you need

#-----------plotting the expected order of convergence------------
x = np.linspace(np.min(dt),np.max(dt),1000)
ax.plot(x, 0.0000006*x**2, color = "green", label = r"~ $\Delta_t^2$", linestyle = "--")

#in case we need linear fit
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_thetadot), color = "red", ax = ax, precisions = [4,4], plot = False)
# ax.plot(x, (x**a)*np.exp(b), color = "red") 
# print(f"a = {a:.3f} et b = {np.exp(b):.9f}")

#-----------formatting------------
ax.legend()    
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "PasExcitation_convergence", ext = ext)





# #-----------portrait de phase------------
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta$ [rad]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)


# #-----------trajectoires numériques------------
# data = np.loadtxt(outputs[-1])  
# theta = data[:, 1]
# thetadot = data[:, 2]
# ax.plot(theta,thetadot, label = r"Trajectoire numérique pour N_{steps} = 5000$", color = "blue", linestyle = "-")

# #-----------trajectoire théorique------------
# t = np.linspace(0,2*T,10000)
# theta_th = theta_0*np.cos(omega_0*t)
# thetadot_th = -theta_0*omega_0*np.sin(omega_0*t)
# ax.plot(theta_th,thetadot_th, label = "Trajectoire théorique", color = "black", linestyle = "--")

# #-----------formatting------------
# ax.legend()    
# u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
# plt.grid(True, which="both", linestyle='--')
# u.savefig(fig, "PasExcitation_PortraitPhase",ext = ext)


#-------------plot theta(t) for different Nsteps-------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\theta$ [rad]', grid_bool=False)

#trajectoires numérique
for i in [0,1,3] : 
    data = np.loadtxt(outputs[i])  
    t = data[:, 0]
    theta = data[:, 1]
    ax.plot(t,theta, label = r"$N_{steps}$" + f" = {nsteps[i]}", linestyle = "-")

#trajectoire théorique
t = np.linspace(0,3*T,10000)
theta_th = theta_0*np.cos(omega_0*t)
ax.plot(t,theta_th, label = "Trajectoire théorique", color = "black", linestyle = "--")

ax.legend()    
ax.set_ylim(-1.1e-10,1.75e-10)
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "PasExcitation_theta(t)",ext = ext)


#-------------plot thetadot(t) for different Nsteps-------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\dot{\theta}$ [rad/s]', grid_bool=False)

#trajectoires numérique
for i in [0,1,3] : 
    data = np.loadtxt(outputs[i])  
    t = data[:, 0]
    thetadot = data[:, 2]
    ax.plot(t,thetadot, label = r"$N_{steps}$" + f" = {nsteps[i]}", linestyle = "-")

#trajectoire théorique
t = np.linspace(0,3*T,10000)
thetadot_th = -theta_0*omega_0*np.sin(omega_0*t)
ax.plot(t,thetadot_th, label = "Trajectoire théorique", color = "black", linestyle = "--")

ax.legend()    
ax.set_ylim(-8e-10,12e-10)
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "PasExcitation_thetadot(t)",ext = ext)



#-----------simulate for different initial conditions and fixed dt------------
Nsteps = 5000
T_ = 4*2*np.pi/omega_0
thetas_0 = np.linspace(pi/100,pi,100,endpoint=False)
thetadot_0 = 0


outputs = []  
#-----------------simulation-----------------
paramstr = 'theta0'  # Parameter name to scan
param = thetas_0  # Parameter values to scan
for i in range(len(thetas_0)):
    output_file = f"./outputs/{paramstr}={param[i]}.out"
    outputs.append(output_file)
    
#à décommenter pour simuler
# for i in range(len(thetas_0)):   
#     cmd = f"{repertoire}{executable} {input_filename} nsteps={5000:.15g} output={output_file} theta0={thetas_0[i]} thetadot0={thetadot_0} tFin={T_}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')


#-----------------finding periods and thetas, thetasdot----------------
Ts = np.zeros(len(thetas_0))
thetas_list = []
thetadots_list = []
for i in range(len(thetas_0)):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t  = data[:, 0]
    theta = data[:, 1]
    thetadot = data[:, 2]
    
    thetas_list.append(theta)
    thetadots_list.append(thetadot)
    
    # Calcul des différences temporelles
    dt = t[1]-t[0]
    # Trouver les pics dans les positions
    peaks, _ = find_peaks(theta)

    if peaks.size != 0:
        T_ = t[peaks[0]]
        Ts[i] = T_


#-----------------plotting----------------
Omegas = 2*np.pi/Ts

ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta_0$ [rad]", ylabel=r'$\Omega$ [rad/s]')
ax.plot(thetas_0,Omegas, label = r"Fréquence propre $\Omega$ en fonction de $\theta_0$", color = "blue", linestyle = "-")


#formule approchée de Borda
T_theta0 = T*(1 + (1/16)*thetas_0**2)
ax.plot(thetas_0,2*np.pi/T_theta0, label = r"$\Omega(\theta_0) = 1/\left(\sqrt{\frac{l}{g}} \left( 1 + \frac{\theta_0^2}{16} \right)\right)$", color = "red", linestyle = "--")
#plot at y = omega_0
ax.plot(thetas_0,omega_0*np.ones(len(thetas_0)), color = "black", linestyle = "--", label = r"Fréquence propre théorique pour de petits angles")
#formule MAG-2
T_theta0 = T*(1 + (1/16)*thetas_0**2 + (11/3072)*thetas_0**4)
ax.plot(thetas_0,2*np.pi/T_theta0, label = r"$\Omega(\theta_0) = 1/\left(\sqrt{\frac{l}{g}} \left( 1 + \frac{\theta_0^2}{16} + \frac{11 \theta_0^4}{3072} \right)\right)$", color = "green", linestyle = "--")

#-----------formatting------------
ax.legend()    
u.set_legend_properties(ax,colors=["black","navy","black"], fontsize = 20, loc = "best",markers=["+","-","."],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "PasExcitation_FreqPropre",ext = ext)



#-----------portrait de phase avec différents theta0------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta$ [rad]", ylabel=r'$\dot{\theta}$ [rad/s]')

for i in range(int(len(thetas_list)/10)) : 
    ax.plot(thetas_list[10*i],thetadots_list[10*i], linestyle = "-", color="navy", linewidth = 1.2)

u.set_legend_properties(ax,colors=["black","navy","black"], fontsize = 20, loc = "best",markers=["+","-","."],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.savefig(fig, "PasExcitation_PortraitPhaseTheta0",ext = ext)