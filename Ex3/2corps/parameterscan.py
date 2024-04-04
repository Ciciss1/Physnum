import sys
# Add the directory containing the module to sys.path
sys.path.append('/workspaces/Physnum')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import utils_v2 as u

from scipy.signal import find_peaks
from scipy.optimize import curve_fit



# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Ex3_2024'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([5000,7500,10000,12500,15000,17500,20000,25000,30000,35000,50000,60000,75000,90000,110000])
nsimul_nsteps = len(nsteps)

tols = np.array([0.01,0.02,0.03,0.05,0.08,0.1,0.2,0.3,0.5,0.8,1,2,3,5,8,10,25,35,50,75,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])
nsimul_ntols = len(tols)

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need


M_earth = 5.972e24
r0 = 6800e3
r1 = 1.5e9
G = 6.67430e-11

x_final_th = r1
y_final_th = 0
Emec_th = - (G*M_earth)/(r1+r0)

print(Emec_th)


# Simulations
outputs_nsteps = []  
outputs_tols = []  

for i in range(nsimul_ntols):
    output_file = f"./outputs/tols={tols[i]}.out"
    outputs_tols.append(output_file)

for i in range(nsimul_nsteps):
    output_file = f"./outputs/Nsteps={nsteps[i]}.out"
    outputs_nsteps.append(output_file)
        
# for i in range(nsimul_nsteps):
#     output_file = outputs_nsteps[i]
#     cmd = f"{repertoire}{executable} {input_filename} adapt=false nsteps={nsteps[i]:.15g} output={output_file}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')
    
# for i in range(nsimul_ntols):
#     output_file = outputs_tols[i]
#     cmd = f"{repertoire}{executable} {input_filename} adapt=true tol={tols[i]:.15g} output={output_file}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')


    

N = 5000
#simulation with varying nsteps
x_nsteps = np.zeros(N)
y_nsteps = np.zeros(N)
t_nsteps = np.zeros(N)
Emec_nsteps = np.zeros(N)

Xs_nsteps=[]
Ys_nsteps=[]
Vs_nsteps=[]
Es_nsteps=[]
ts = []

Emec_errors_nsteps = np.zeros(nsimul_nsteps)
pos_errors_nsteps = np.zeros(nsimul_nsteps)

for i in range(nsimul_nsteps):
    data = np.loadtxt(outputs_nsteps[i])
    
    if nsteps[i] == N:
        t_nsteps = data[:,0]*3.8052e-7
        x_nsteps = data[:,1]  
        y_nsteps = data[:,2]
        Emec_nsteps = data[:,5]
    
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    vx = data[:,3]
    vy = data[:,4]
    E = data[:,5]
    
    Xs_nsteps.append(x)
    Ys_nsteps.append(y)
    Vs_nsteps.append(vx**2 + vy**2)
    ts.append(t*3.8052e-7)
    Es_nsteps.append(E)
    
    Emec = data[:,5]
    Epot = -G*M_earth/np.sqrt(x**2 + y**2)
    Ekin = 1/2 * (vx**2 + vy**2)   
    
    Emec_errors_nsteps[i] = np.abs(np.max(Emec) - np.min(Emec))
    
    lastx = data[-1,1]
    lasty = data[-1,2]
    
    pos_errors_nsteps[i] = np.sqrt((lastx - x_final_th)**2 + (lasty - y_final_th)**2)
    
    
    
    
    
#simulation with varying tol
x_tols = 0
y_tols = 0
t_tols = 0
Emec_tols = 0

tols_nsteps = np.zeros(nsimul_ntols)

Emec_errors_tols = np.zeros(nsimul_ntols)
pos_errors_tols = np.zeros(nsimul_ntols)

for i in range(nsimul_ntols):
    data = np.loadtxt(outputs_tols[i])
    
    if tols[i] == 100:
        t_tols = 3.8052e-7*data[:,0]
        x_tols = data[:,1]  
        y_tols = data[:,2]
        v2_tols = data[:,3]**2 + data[:,4]**2
        Emec_tols = data[:,5]
    
    Emec = data[:,5]
    Emec_errors_tols[i] = np.abs(np.max(Emec) - np.min(Emec))

    lastx = data[-1,1]
    lasty = data[-1,2]
    
    pos_errors_tols[i] = np.sqrt((lastx - x_final_th)**2 + (lasty - y_final_th)**2)
    
    tols_nsteps[i] = data[-1,-1]
    
    
    

ext = "png"    
    
    
rect_xlim = [-0.02e9, 0.02e9]
rect_ylim = [-0.05e8, 0.25e8]
#---------------------------------NSTEPS---------------------------------#
#plotting the position for varying nsteps
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')

ax.plot(x_nsteps, y_nsteps, linestyle='-', marker = "+" , color='navy')

# #plotting a rectangle on the zoom region
# ax.add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))


# ax.axis('equal')
plt.tight_layout()
u.savefig(fig, "2corps_xy_nsteps", ext = ext)



# #plotting the zoomed region
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')
# ax.plot(x_nsteps[:500], y_nsteps[:500], linestyle='-', marker = "+" , color='navy')
# ax.scatter(0, 0, s=250, color='black', marker='o', label='Earth',zorder =15)
# ax.text(0, 0.2e7, 'Earth', fontsize=18, ha='right', va='bottom', zorder = 15)
# ax.set_xlim(rect_xlim)
# ax.set_ylim(rect_ylim)
# plt.tight_layout()
# u.savefig(fig, "2corps_xy_nsteps_zoom", ext = ext)


#plotting the energy for varying nsteps
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [mois]", ylabel=r'$E$ [J]')

ax.plot(t_nsteps, Emec_nsteps, linestyle='-', marker = "+" , color='navy', label = "Valeur numérique")
#plotting the theoretical value
ax.hlines(Emec_th, t_nsteps[0], t_nsteps[-1], linestyle='--', color='black', label = "Valeur théorique")
u.set_legend_properties(ax, fontsize=18)
plt.tight_layout()

u.savefig(fig, "2corps_E(t)_nsteps", ext = ext)



# #---------------------------------TOLS---------------------------------#
#plotting the position for varying tols
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')

# #plotting a rectangle on the zoom region
# ax.add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))


ax.plot(x_tols, y_tols, linestyle='-', marker = "+" , color='navy')

# ax.axis('equal')
plt.tight_layout()
u.savefig(fig, "2corps_xy_ntols", ext = ext)

# #plotting the zoomed region
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')
# ax.plot(x_tols, y_tols, linestyle='-', marker = "+" , color='navy')
# ax.scatter(0, 0, s=250, color='black', marker='o', label='Earth',zorder =15)
# ax.text(0, 0.2e7, 'Earth', fontsize=18, ha='right', va='bottom', zorder = 15)
# ax.set_xlim(rect_xlim)
# ax.set_ylim(rect_ylim)
# plt.tight_layout()
# u.savefig(fig, "2corps_xy_ntols_zoom", ext = ext)


#plotting the energy for varying tols
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [mois]", ylabel=r'$E$ [J]')

ax.plot(t_tols, Emec_tols, linestyle='-', marker = "+" , color='navy', label = "Valeur numérique")
ax.hlines(Emec_th, t_nsteps[0], t_nsteps[-1], linestyle='--', color='black', label = "Valeur théorique")
plt.tight_layout()

u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "2corps_E(t)_tols", ext = ext)

print("Nsteps : ", tols_nsteps[4])


rect_xlim2 = [-0.01e9, 0.065e9]
rect_ylim2 = [-0.05e8, 0.5e8]

#---------------------------------Plot xy zoomed for both methods---------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')
ax.plot(x_tols, y_tols, linestyle='-', marker = "d" ,markersize = 10, color='blue', label = r"RK4 $\Delta t$ adaptatif")
ax.plot(x_nsteps[:500], y_nsteps[:500], linestyle='-', marker = "^" ,markersize = 10, color='red', label = r"RK4 $\Delta t$ fixe")
ax.scatter(0, 0, s=250, color='black', marker='o', zorder =15)
ax.text(0.7e7, 0.2e7, 'Earth', fontsize=18, ha='right', va='bottom', zorder = 15)
ax.set_xlim(rect_xlim2)
ax.set_ylim(rect_ylim2)
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "2corps_xy_zoom", ext = ext)




#---------------------------------Plot the trajectory for various nsteps---------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')
for i in [0,1,2,13]:
    ax.plot(Xs_nsteps[i], Ys_nsteps[i], linestyle='-', label = f"Nsteps = {nsteps[i]}")

ax.plot(x_tols, y_tols, linestyle='--', color = "black",label = r"RK4 $\Delta t$ adaptatif")

# ax.axis('equal')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "2corps_xy_var", ext = ext)

#---------------------------------Plot r for various nsteps--------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$t$ [mois]", ylabel=r'$r$ [m]')
for i in [0,1,2,13]:
    ax.plot(ts[i], np.sqrt(Xs_nsteps[i]**2 + Ys_nsteps[i]**2), linestyle='-', label = f"Nsteps = {nsteps[i]}")
    
ax.plot(t_tols, np.sqrt(x_tols**2 + y_tols**2), linestyle='--', color = "black",label = r"RK4 $\Delta t$ adaptatif")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "2corps_r_var", ext = ext)

#---------------------------------Plot the energy for various nsteps---------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [mois]", ylabel=r'$E$ [J]')
for i in [0,1,2,13]:
    ax.plot(ts[i], Es_nsteps[i], linestyle='-', label = f"Nsteps = {nsteps[i]}")
    
ax.plot(t_tols, Emec_tols, linestyle='--', color = "black",label = r"RK4 $\Delta t$ adaptatif")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "2corps_E_var", ext = ext)


# #---------------------------------Plot the speed---------------------------------#
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$v^2$ [m/s]')

# ax.plot(t_tols, v2_tols, linestyle='-', color='blue', label = r"RK4 $\Delta t$ adaptatif")
# for i in [0,1]:
#     ax.plot(ts[i], Vs_nsteps[i], linestyle='-', label = r"N_{steps} = " + f"{nsteps[i]}")

# plt.tight_layout()
# u.set_legend_properties(ax, fontsize=18)
# u.savefig(fig, "2corps_v", ext = ext)


#---------------------------------ERRORS---------------------------------#
#pos
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r'Erreur sur la position finale [m]')

ax.loglog(nsteps, pos_errors_nsteps, linestyle='-' ,marker = "+",  label = r"RK4 $\Delta t$ fixe", color='blue')
ax.loglog(tols_nsteps, pos_errors_tols, linestyle='-', marker = "x", label = r"RK4 $\Delta t$ adaptatif" , color='red')

xlim = [30, 150000]
ylim = [180,3e9]


x = np.linspace(10, 1000, 1000)
y = 5e12/x**4
ax.loglog(x,y, linestyle='--', color='black', label=r'$ \sim N^{-4}$')

x = np.linspace(10, 130000, 1000)
y = 1e23/x**4
ax.loglog(x,y, linestyle='--', color='black')

ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.grid(which='both', linestyle='--', linewidth=0.5)
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "Erreur_position_finale", ext = ext)


#Emec
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r"Erreur sur l'énergie mécanique [J]")

ax.loglog(nsteps, Emec_errors_nsteps, linestyle='-' ,marker = "+",  label = r"RK4 $\Delta t$ fixe", color='blue')
ax.loglog(tols_nsteps, Emec_errors_tols, linestyle='-', marker = "x", label = r"RK4 $\Delta t$ adaptatif" , color='red')


xlim = [30, 150000]
ylim = [0.1, 2e7]


x = np.linspace(10, 1000, 1000)
y = 5e8/x**4
ax.loglog(x,y, linestyle='--', color='black', label=r'$ \sim N^{-4}$')

x = np.linspace(10, 130000, 1000)
y = 4e19/x**4
ax.loglog(x,y, linestyle='--', color='black')

ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.grid(which='both', linestyle='--', linewidth=0.5)
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "Erreur_Energie", ext = ext)

