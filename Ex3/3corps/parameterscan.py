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


nsteps = np.array([5000,5125,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,8750,10000,12500,15000,17500,20000,25000,30000,35000,50000,60000])
nsimul_nsteps = len(nsteps)

tols = np.array([0.1,0.2,0.3,0.5,0.8,1,2,3,5,8,10,20,30,55,85,115,170,200,255,300,355,405,450,500,593,645,700,750,820])
nsimul_ntols = len(tols)

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need



M_S = 1.98892e30
M_T = 5.9736e24
d = 149.598023e9
G = 6.674e-11
L2x = 1.51099098533039886475e+11
L2y = 0
alpha = M_T / (M_S + M_T)
beta = M_S / (M_S + M_T)
X_S = -alpha * d
X_T = beta * d

x0 = L2x
y0 = L2y
vx0 = 0
vy0 = 0.1

Omega = 1.99119425021e-07


Ekin0 = 0.5*((vx0**2 + vy0**2) + Omega**2*(x0**2 + y0**2) + 2*Omega*(x0*vy0 - y0*vx0))
Ekin0 = 0.5*((vx0**2 + vy0**2) - Omega**2*(x0**2 + y0**2))
# Ekin0 = 0
Epot0 = -G*M_T/np.sqrt((x0-X_T)**2 + y0**2) - G*M_S/np.sqrt((x0-X_S)**2 + y0**2)
Emec_th = Ekin0 + Epot0


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
ts = []

Emec_errors_nsteps = np.zeros(nsimul_nsteps)
pos_errors_nsteps = np.zeros(nsimul_nsteps)

lastxs_nsteps = np.zeros(nsimul_nsteps)
lastys_nsteps = np.zeros(nsimul_nsteps)

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
    
    Xs_nsteps.append(x)
    Ys_nsteps.append(y)
    Vs_nsteps.append(vx**2 + vy**2)
    ts.append(t*3.8052e-7)
    
    Emec = data[:,5]
    Epot = -G*M_T/np.sqrt(x**2 + y**2)
    Ekin = 1/2 * (vx**2 + vy**2)   
    
    Emec_errors_nsteps[i] = np.abs(np.max(Emec) - np.min(Emec))
    
    lastx = data[-1,1]
    lasty = data[-1,2]
    
    lastxs_nsteps[i] = lastx
    lastys_nsteps[i] = lasty
    
    # pos_errors_nsteps[i] = np.sqrt((lastx - x_final_th)**2 + (lasty - y_final_th)**2)
    
    
    
    
#simulation with varying tol
x_tols = 0
y_tols = 0
t_tols = 0
Emec_tols = 0

lastxs_tols = np.zeros(nsimul_ntols)
lastys_tols = np.zeros(nsimul_ntols)

tols_nsteps = np.zeros(nsimul_ntols)

Emec_errors_tols = np.zeros(nsimul_ntols)

max_speed_list = np.zeros(nsimul_ntols)

Xs_ntols=[]
Ys_ntols=[]
Vs_ntols=[]
ts_ntols = []

for i in range(nsimul_ntols):
    data = np.loadtxt(outputs_tols[i])
    
    if tols[i] == 115:
        t_tols = 3.8052e-7*data[:,0]
        x_tols = data[:,1]  
        y_tols = data[:,2]
        v2_tols = data[:,3]**2 + data[:,4]**2
        Emec_tols = data[:,5]
    
    Emec = data[:,5]
    Emec_errors_tols[i] = np.abs(np.max(Emec) - np.min(Emec))

    lastx = data[-1,1]
    lasty = data[-1,2]
    
    lastxs_tols[i] = lastx
    lastys_tols[i] = lasty
    
    max_speed_list[i] = np.max(data[:,3]**2 + data[:,4]**2)
    
    # pos_errors_tols[i] = np.sqrt((lastx - x_final_th)**2 + (lasty - y_final_th)**2)
    
    tols_nsteps[i] = data[-1,6]
    
    Xs_ntols.append(data[:,1])
    Ys_ntols.append(data[:,2])
    ts_ntols.append(data[:,0]*3.8052e-7)
    
    
    
    
ext = "pdf"    
    
    
rect_xlim = [-0.02e9, 0.02e9]
rect_ylim = [-0.05e8, 0.25e8]
#---------------------------------NSTEPS---------------------------------#
#plotting the position for varying nsteps
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')

ax.plot(x_nsteps, y_nsteps, linestyle='-', marker = "+" , color='navy')
ax.plot(X_T,0, marker = "o", color='black', markersize = 10, label = "Terre", zorder = 15)
ax.plot(X_S,0, marker = "o", color='black', markersize = 10, label = "Soleil", zorder = 15)


ax.axis('equal')
plt.tight_layout()
u.savefig(fig, "3corps_xy_nsteps", ext = ext)


#plotting the energy for varying nsteps
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [mois]", ylabel=r'$E$ [J]')

ax.plot(t_nsteps, Emec_nsteps, linestyle='-', marker = "+" , color='navy', label = "Valeur numérique")
ax.hlines(Emec_th, t_nsteps[0], t_nsteps[-1], linestyle='--', color='black', label = "Energie théorique")
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "3corps_E(t)_nsteps", ext = ext)



# #---------------------------------TOLS---------------------------------#
#plotting the position for varying tols
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')


ax.plot(x_tols, y_tols, linestyle='-', marker = "+" , color='navy')
ax.plot(X_T,0, marker = "o", color='black', markersize = 10, label = "Terre", zorder = 15)
ax.text(X_T + 0.2e11,0, "Terre", fontsize=18, verticalalignment='bottom', horizontalalignment='right', zorder = 15)
ax.plot(X_S,0, marker = "o", color='black', markersize = 10, label = "Soleil", zorder = 15)
ax.text(X_S + 0.2e11,0, "Soleil", fontsize=18, verticalalignment='bottom', horizontalalignment='right', zorder = 15)

# ax.set_xlim(1.45e11, 1.55e11)
# ax.set_ylim(-0.05e11, 0.05e11)

ax.axis('equal')
plt.tight_layout()
u.savefig(fig, "3corps_xy_ntols", ext = ext)

#plotting the energy for varying tols
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [mois]", ylabel=r'$E$ [J]')

ax.plot(t_tols, Emec_tols, linestyle='-', marker = "+" , color='navy', label = "Valeur numérique")
ax.hlines(Emec_th, t_nsteps[0], t_nsteps[-1], linestyle='--', color='black', label = "Energie théorique")
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "3corps_E(t)_tols", ext = ext)

print("Nsteps : ", tols_nsteps[4])


rect_xlim2 = [-0.01e9, 0.065e9]
rect_ylim2 = [-0.05e8, 0.5e8]




#---------------------------------Plot the trajectory for various nsteps---------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$x$ [m]", ylabel=r'$y$ [m]')
for i in [0,3,13]:
    ax.plot(Xs_nsteps[i], Ys_nsteps[i], linestyle='-', label = f"Nsteps = {nsteps[i]}")

ax.plot(x_tols, y_tols, linestyle='--', color = "black",label = r"RK4 $\Delta t$ adaptatif, $N_{steps}$ = " + f"{tols_nsteps[4]}")

# ax.set_xlim(1.45e11, 1.55e11)


ax.scatter(X_T,0, marker = "o", color='black', s = 100, label = "Terre", zorder = 15)

# ax.axis('equal')
ax.set_ylim(-10e8, 15e8)
ax.set_xlim(1.48e11,1.515e11)
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "3corps_xy_var", ext = ext)



#---------------------------------Plot r (log))--------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$t$ [mois]", ylabel=r'Distance à L2 [m]', grid_bool=False)

# for i in [8,13,14,15,16,17]:
#     ax.plot(ts[i], np.sqrt((Xs_nsteps[i]-L2x)**2 + Ys_nsteps[i]**2), linestyle='-', label = f"Nsteps = {nsteps[i]}")

# ax.plot(ts[8], np.sqrt((Xs_nsteps[8]-L2x)**2 + Ys_nsteps[8]**2), linestyle='-', label = f"Nsteps = {nsteps[8]}")

dist = np.sqrt((Xs_ntols[0]-L2x)**2 + Ys_ntols[0]**2)
ax.plot(ts_ntols[0][10:], dist[10:], linestyle='-', color = "navy",label = rf"RK4 $\Delta t$ adaptatif, Nsteps = {tols_nsteps[0]}")

N1 = 10
N2 = 450

print(ts_ntols[0][N1]*30,ts_ntols[0][N2]*30)

ax.scatter([ts_ntols[0][N1],ts_ntols[0][N2]],[dist[N1],dist[N2]], marker = "^", color='black', s = 100, label = "Extrémités du fit", zorder = 15)
# ax.scatter(ts_ntols[0][N2],dist[N2], marker = "^", color='black', s = 100, label = "Fin du fit", zorder = 15)

X = ts_ntols[0][N1:N2]
Y = dist[N1:N2]
a,_,b,_ = u.linear_fit(X,np.log(Y),"red",[3,3],ax=ax,plot=False)

t =np.linspace(0,ts_ntols[0][N2+500],1000)

ax.plot(t,np.exp(a*t+b), linestyle='--', color='red', label = rf"$y =e^{{ {a:.2f}  t + {b:.0f}}}$")

ax.set_yscale('log')

ax.grid(which='both', linestyle='--', linewidth=1)
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1,loc = "lower right")
u.savefig(fig, "3corps_dist_L2_log", ext = ext)



#---------------------------------Plot r (non log)--------------------------------#
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$t$ [mois]", ylabel=r'Distance à L2 [m]')

# for i in [8,13,14,15,16,17]:
#     ax.plot(ts[i], np.sqrt((Xs_nsteps[i]-L2x)**2 + Ys_nsteps[i]**2), linestyle='-', label = f"Nsteps = {nsteps[i]}")

# ax.plot(ts[8], np.sqrt((Xs_nsteps[8]-L2x)**2 + Ys_nsteps[8]**2), linestyle='-', label = f"Nsteps = {nsteps[8]}")

dist = np.sqrt((Xs_ntols[0]-L2x)**2 + Ys_ntols[0]**2)
ax.plot(ts_ntols[0][10:], dist[10:], linestyle='-', color = "navy",label = rf"RK4 $\Delta t$ adaptatif, Nsteps = {tols_nsteps[0]}")

N1 = 10
N2 = 450

ax.scatter([ts_ntols[0][N1],ts_ntols[0][N2]],[dist[N1],dist[N2]], marker = "^", color='black', s = 100, label = "Extrémités du fit", zorder = 15)
# ax.scatter(ts_ntols[0][N2],dist[N2], marker = "^", color='black', s = 100, label = "Fin du fit", zorder = 15)

X = ts_ntols[0][N1:N2]
Y = dist[N1:N2]
a,_,b,_ = u.linear_fit(X,np.log(Y),"red",[3,3],ax=ax,plot=False)

t =np.linspace(0,ts_ntols[0][N2+150],1000)

ax.plot(t,np.exp(a*t+b), linestyle='--', color='red', label = rf"$y =e^{{ {a:.2f}  t + {b:.0f}}}$")

# ax.set_yscale('log')


plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "3corps_dist_L2_nlog", ext = ext)






#---------------------------------ERRORS---------------------------------#
#pos
final_pos_nsteps = np.sqrt(lastxs_nsteps**2 + lastys_nsteps**2)
final_pos_tols = np.sqrt(lastxs_tols**2 + lastys_tols**2)

#nsteps
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\frac{1}{N_{steps}^4}$", ylabel=r'Position finale [m]')
ax.plot(1/nsteps**4, final_pos_nsteps, linestyle='-' ,marker = "+",  label = r"RK4 $\Delta t$ fixe", color='blue')

ax.grid(which='both', linestyle='--', linewidth=0.5)
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "3corps_Erreur_position_finale_Nsteps", ext = ext)

#ntols
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\frac{1}{N_{steps}^5}$", ylabel=r'Position finale [m]')

ax.plot(1/tols_nsteps**5, final_pos_tols, linestyle='-', marker = "x", label = r"RK4 $\Delta t$ adaptatif" , color='red')

ax.grid(which='both', linestyle='--', linewidth=0.5)
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "3corps_Erreur_position_finale_Ntols", ext = ext)



# #---------------------------------ERRORS 2---------------------------------#
# #pos
# final_pos_nsteps = np.sqrt((lastxs_nsteps-L2x)**2 + (lastys_nsteps-L2y)**2)
# final_pos_tols = np.sqrt((lastxs_tols-L2x)**2 + (lastys_tols-L2y)**2)

# #nsteps
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r'Distance finale par rapport à L2 [m]')
# ax.plot(nsteps, final_pos_nsteps, linestyle='-' ,marker = "+",  label = r"RK4 $\Delta t$ fixe", color='blue')

# ax.grid(which='both', linestyle='--', linewidth=0.5)
# u.set_legend_properties(ax, fontsize=18)
# u.savefig(fig, "Erreur_position_finale_L2_Nsteps", ext = ext)

# #ntols
# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r'Distance finale par rapport à L2 [m]')
# ax.plot(tols_nsteps, final_pos_tols, linestyle='-', marker = "x", label = r"RK4 $\Delta t$ adaptatif" , color='red')

# ax.grid(which='both', linestyle='--', linewidth=0.5)
# u.set_legend_properties(ax, fontsize=18)
# u.savefig(fig, "Erreur_position_finale_L2_Ntols", ext = ext)






#Emec
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r"Erreur sur l'énergie par unité de masse [J/kg]")

ax.loglog(nsteps, Emec_errors_nsteps, linestyle='-' ,marker = "+",  label = r"RK4 $\Delta t$ fixe", color='blue')
# ax.grid(which='both', linestyle='--', linewidth=0.5)
# u.set_legend_properties(ax, fontsize=18)
# u.savefig(fig, "3corps_Erreur_Energie_Nsteps", ext = ext)


# ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$N_{steps}$", ylabel=r"Erreur sur l'énergie mécanique [J]")

ax.loglog(tols_nsteps, Emec_errors_tols, linestyle='-', marker = "x", label = r"RK4 $\Delta t$ adaptatif" , color='red')
ax.grid(which='both', linestyle='--', linewidth=0.5)

xlim = [30, 150000]
ylim = [5e-4,3e3]


x = np.linspace(10, 10000, 1000)
y = 2e12/x**4
ax.loglog(x,y, linestyle='--', color='black', label=r'$ \sim N^{-4}$')

x = np.linspace(10, 13000000, 1000)
y = 1e18/x**4
ax.loglog(x,y, linestyle='--', color='black')

ax.set_xlim(xlim)
ax.set_ylim(ylim)

u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "3corps_Erreur_Energie", ext = ext)

# xlim = [30, 150000]
# ylim = [0.1, 2e7]


# # x = np.linspace(10, 1000, 1000)
# # y = 5e8/x**4
# # ax.loglog(x,y, linestyle='--', color='black', label=r'$ \sim N^{-4}$')

# # x = np.linspace(10, 130000, 1000)
# # y = 4e19/x**4
# # ax.loglog(x,y, linestyle='--', color='black')

# ax.set_xlim(xlim)
# ax.set_ylim(ylim)


print(max_speed_list)