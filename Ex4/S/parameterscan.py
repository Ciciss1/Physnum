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

ext = "pdf"


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice4_2024'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


S0 = 1000
kappa0 = 0.001
R = 0.05
T_R = 293

# def Temperature(r) : 
#     C1 = (S0*R**2)/(2*kappa0)
#     C2 = T_R +((S0*R**2)/(2*kappa0))*(0.5 -np.log(R))
#     return -S0*r**2/(4*kappa0) + C1*np.log(r) + C2

def Temperature(r) : 
    C2 = T_R +((S0*R**2)/(4*kappa0))
    return -S0*r**2/(4*kappa0) + C2


def heat_flux(r) : 
    return S0*r/(2) 
    
    


Ns = np.logspace(np.log10(10), np.log10(10000), 100, dtype=int)

nsimul_Ns = len(Ns)
# Simulations
outputs_Ns_uniform = []   
outputs_Ns_nonuniform = [] 

#create output files for uniform case
for i in range(nsimul_Ns):
    output_file = f"./outputs/N={Ns[i]}_uniform"
    outputs_Ns_uniform.append(output_file)

#create output files for nonuniform case
for i in range(nsimul_Ns):
    output_file = f"./outputs/N={Ns[i]}_nonuniform"
    outputs_Ns_nonuniform.append(output_file)
        

# #simulate for uniform cases
# for i in range(nsimul_Ns):
#     output_file = outputs_Ns_uniform[i]
#     cmd = f"{repertoire}{executable} {input_filename} N={ Ns[i]:.15g} output={output_file} IsUniform=true"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')
    

# #simulate for nonuniform cases
# for i in range(nsimul_Ns):
#     output_file = outputs_Ns_nonuniform[i]
#     cmd = f"{repertoire}{executable} {input_filename} N={ Ns[i]:.15g} output={output_file} IsUniform=false"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')
    
    
#-------------------(i) uniform case T(r)-------------------

#read data
outputs_Ns_uniform_T = [f"./outputs/N={Ns[i]}_uniform_T.out" for i in range(nsimul_Ns)]
outputs_Ns_uniform_heat = [f"./outputs/N={Ns[i]}_uniform_heat.out" for i in range(nsimul_Ns)]

data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[0])
data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[0])

T_r_num_uniform = data_uniform_T[:,1]
r_T_uniform = data_uniform_T[:,0]

heat_r_num = data_uniform_heat[:,1]
r_heat = data_uniform_heat[:,0]



#read data
outputs_Ns_nonuniform_T = [f"./outputs/N={Ns[i]}_nonuniform_T.out" for i in range(nsimul_Ns)]
outputs_Ns_nonuniform_heat = [f"./outputs/N={Ns[i]}_nonuniform_heat.out" for i in range(nsimul_Ns)]

data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[0])
data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[0])

T_r_num_nonuniform = data_nonuniform_T[:,1]
r_T_nonuniform = data_nonuniform_T[:,0]

heat_r_num = data_nonuniform_heat[:,1]
r_heat = data_nonuniform_heat[:,0]

#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'$T$ [K]')


ax.plot(r_T_uniform,T_r_num_uniform, marker="+", label = "Équidistant", color = "green")

r = np.linspace(0,R,1000)
ax.plot(r_T_nonuniform,T_r_num_nonuniform, marker="+", label = "Quadratique", color = "navy")
ax.plot(r,Temperature(r),  label = "Analytique", color = "red", linestyle = "--")

ax.set_ylim(200,1050)

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "N10_T(r)", ext = ext)




#-------------------(i) error on numerical and analytical T(r)-------------------

error_uniform = np.abs(T_r_num_uniform - Temperature(r_T_uniform))
error_nonuniform = np.abs(T_r_num_nonuniform - Temperature(r_T_nonuniform))

#plot error
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'Erreur sur $T$ [K]')

ax.plot(r_T_uniform,error_uniform, marker="+", label = "Équidistant", color = "navy")
ax.plot(r_T_nonuniform,error_nonuniform, marker="x", label = "Quadratique", color = "red")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "N10_error_T", ext = ext)


#-------------------(i) heat(r)-------------------
#read data
outputs_Ns_uniform_heat = [f"./outputs/N={Ns[i]}_uniform_heat.out" for i in range(nsimul_Ns)]
data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[0])
heat_r_num_uniform = data_uniform_heat[:,1]
r_heat_uniform = data_uniform_heat[:,0]

#read data
outputs_Ns_nonuniform_heat = [f"./outputs/N={Ns[i]}_nonuniform_heat.out" for i in range(nsimul_Ns)]
data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[0])
heat_r_num_nonuniform = data_nonuniform_heat[:,1]
r_heat_nonuniform = data_nonuniform_heat[:,0]

#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'$j_Q$ [W/m$^2$]')

ax.plot(r_heat_uniform, heat_r_num_uniform, marker="+", label = "Équidistant", color = "navy")
ax.plot(r_heat_nonuniform, heat_r_num_nonuniform, marker="+", label = "Quadratique", color = "green")
ax.plot(r_heat, heat_flux(r_heat),  label = "Analytique", color = "red", linestyle = "--")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "N10_heat(r)", ext = ext)

#-------------------(i) error on numerical and analytical T(r)-------------------
error_uniform = np.abs(heat_r_num_uniform - heat_flux(r_heat_uniform))
error_nonuniform = np.abs(heat_r_num_nonuniform - heat_flux(r_heat_nonuniform))
#plot error
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'Erreur sur $j_Q$ [W/m$^2$]')
ax.plot(r_heat_uniform,error_uniform, marker="+", label = "Équidistant", color = "navy")
ax.plot(r_heat_nonuniform,error_nonuniform, marker="x", label = "Quadratique", color = "red")
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "N10_error_heat", ext = ext)


#-------------------(i) power balance-------------------
P_tot = S0*np.pi*R**2

boundary_heat_flux_uniform = 2*np.pi*R*heat_r_num_uniform[-1]
boundary_heat_flux_nonuniform = 2*np.pi*R*heat_r_num_nonuniform[-1]

print("P_tot = ", P_tot)
print("boundary_heat_flux_uniform = ", boundary_heat_flux_uniform)
print("boundary_heat_flux_nonuniform = ", boundary_heat_flux_nonuniform)
print("relative error uniform = ",100*abs(P_tot - boundary_heat_flux_uniform)/P_tot)
print("relative error nonuniform = ",100*abs(P_tot - boundary_heat_flux_nonuniform)/P_tot)


#-------------------(ii) convergence study-------------------

errors_T0_uniform = np.zeros(nsimul_Ns)
errors_heat_R_uniform = np.zeros(nsimul_Ns)
errors_boundary_heat_R_uniform = np.zeros(nsimul_Ns)

errors_T0_nonuniform = np.zeros(nsimul_Ns)
errors_heat_R_nonuniform = np.zeros(nsimul_Ns)
errors_boundary_heat_R_nonuniform = np.zeros(nsimul_Ns)

first_hs = np.zeros(nsimul_Ns)

#Iterating throught
for i in range(nsimul_Ns):
    data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[i])
    data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[i])
    data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[i])
    data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[i])
    
    
    T0_th = Temperature(0)
    heat_R_th = heat_flux(R)
    boundary_heat_R_th = P_tot
    
    #uniform case
    T0 = data_uniform_T[0,1]
    heat_R = data_uniform_heat[-1,1]
    boundary_heat_R = 2*np.pi*R*heat_R
    
        #errors
    error_T0 = abs(T0 - T0_th)
    error_heat_R = abs(heat_R - heat_R_th)
    error_boundary_heat_R = abs(boundary_heat_R - boundary_heat_R_th)

    errors_T0_uniform[i] = error_T0
    errors_heat_R_uniform[i] = error_heat_R
    errors_boundary_heat_R_uniform[i] = error_boundary_heat_R
    
    #nonuniform case
    T0 = data_nonuniform_T[0,1]
    heat_R = data_nonuniform_heat[-1,1]
    boundary_heat_R = 2*np.pi*R*heat_R
    
    first_hs[i] = data_nonuniform_T[1,0]- data_nonuniform_T[0,0]
        
        #errors 
    error_T0 = abs(T0 - T0_th)
    error_heat_R = abs(heat_R - heat_R_th)
    error_boundary_heat_R = abs(boundary_heat_R - boundary_heat_R_th)
    
    errors_T0_nonuniform[i] = error_T0
    errors_heat_R_nonuniform[i] = error_heat_R
    errors_boundary_heat_R_nonuniform[i] = error_boundary_heat_R
    
#plot error T0
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'Erreur sur $T(0)$ [K]', grid_bool = False)

ax.loglog(1/Ns,errors_T0_uniform, marker="+", label = "Équidistant", color = "navy")
ax.loglog(1/Ns,errors_T0_nonuniform, marker="x", label = "Quadratique", color = "red")

x = np.linspace(1e-4,1e-1,100)
ax.loglog(x,5e2*x**2,linestyle = "--", color = "black")
ax.text(8e-3,1e-2,r"~ $\frac{1}{N^2}$", fontsize = 18)

x = np.linspace(1e-4,1e-1,100)
ax.loglog(x,2e2*x,linestyle = "--", color = "black")
ax.text(8e-4,5e-2,r"~ $\frac{1}{N}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "convergence_T0", ext = ext)



#error on heat_R
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'Erreur sur $j_Q(R)$ [W/m$^2$]', grid_bool = False)

ax.loglog(1/Ns,errors_heat_R_uniform, marker="+", label = "Équidistant", color = "navy")
ax.loglog(1/Ns,errors_heat_R_nonuniform, marker="x", label = "Quadratique", color = "red")

x = np.linspace(1e-4,1e-1,100)
ax.loglog(x,3e0*x,linestyle = "--", color = "black")
ax.text(6e-3,1e-2,r"~ $\frac{1}{N}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "convergence_heat_R", ext = ext)


#error on boundary_heat_R
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'Erreur sur $P_{\mathrm{tot}}$ [W]', grid_bool = False)

ax.loglog(1/Ns,errors_boundary_heat_R_uniform, marker="+", label = "Équidistant", color = "navy")  
ax.loglog(1/Ns,errors_boundary_heat_R_nonuniform, marker="x", label = "Quadratique", color = "red")

x = np.linspace(1e-4,1e-1,100)
ax.loglog(x,9e-1*x,linestyle = "--", color = "black")
ax.text(2e-2,1e-2,r"~ $\frac{1}{N}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "convergence_boundary_heat_R", ext = ext)



#error on T0 with first_Hs
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$h_1$ [m]", ylabel=r'Erreur sur $T(0)$ [K]', grid_bool = False)

ax.loglog(first_hs,errors_T0_nonuniform, marker="x", label = "Quadratique", color = "red")

x = np.linspace(np.min(first_hs),np.max(first_hs),100)
ax.loglog(x,9e4*x**2,linestyle = "--", color = "black")
ax.text(4.5e-3,1e0,r"~ $\frac{1}{N^2}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "convergence_T0_h1", ext = ext)





