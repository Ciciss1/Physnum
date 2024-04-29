import sys
# Add the directory containing the module to sys.path
sys.path.append('/workspaces/Physnum')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import math as m
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

r0 = 0.03
sigma = 0.005


Ns = np.logspace(np.log10(100), np.log10(10000), 100, dtype=int)

# Ns=np.array([10,15,20,30,40,50,75,100,200,300,400,500,750,1000,1500,2000,3250,5000,10000])
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

#read data for uniform case
outputs_Ns_uniform_T = [f"./outputs/N={Ns[i]}_uniform_T.out" for i in range(nsimul_Ns)]
outputs_Ns_uniform_heat = [f"./outputs/N={Ns[i]}_uniform_heat.out" for i in range(nsimul_Ns)]


#N10
N10_data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[0])
N10_data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[0])

N10_T_r_num_uniform = N10_data_uniform_T[:,1]
N10_r_T_uniform = N10_data_uniform_T[:,0]

N10_heat_r_num = N10_data_uniform_heat[:,1]
N10_r_heat = N10_data_uniform_heat[:,0]


#N1000
N1000_data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[-1])
N1000_data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[-1])

N1000_T_r_num_uniform = N1000_data_uniform_T[:,1]
N1000_r_T_uniform = N1000_data_uniform_T[:,0]

N1000_heat_r_num = N1000_data_uniform_heat[:,1]
N1000_r_heat = N1000_data_uniform_heat[:,0]



#read data for non uniform case
outputs_Ns_nonuniform_T = [f"./outputs/N={Ns[i]}_nonuniform_T.out" for i in range(nsimul_Ns)]
outputs_Ns_nonuniform_heat = [f"./outputs/N={Ns[i]}_nonuniform_heat.out" for i in range(nsimul_Ns)]


#N10
N10_data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[0])
N10_data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[0])

N10_T_r_num_nonuniform = N10_data_nonuniform_T[:,1]
N10_r_T_nonuniform = N10_data_nonuniform_T[:,0]

N10_heat_r_num = N10_data_nonuniform_heat[:,1]
N10_r_heat = N10_data_nonuniform_heat[:,0]


#N1000
N1000_data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[-1])
N1000_data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[-1])

N1000_T_r_num_nonuniform = N1000_data_nonuniform_T[:,1]
N1000_r_T_nonuniform = N1000_data_nonuniform_T[:,0]

N1000_heat_r_num = N1000_data_nonuniform_heat[:,1]
N1000_r_heat = N1000_data_nonuniform_heat[:,0]



#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'$T$ [K]')


ax.plot(N10_r_T_uniform,N10_T_r_num_uniform, marker="+", label = "Équidistant, N=10", color = "navy")
ax.plot(N10_r_T_nonuniform,N10_T_r_num_nonuniform, marker="+", label = "Quadratique, N=10", color = "red")

ax.plot(N1000_r_T_nonuniform,N1000_T_r_num_nonuniform, linestyle = "--" ,label = "Quadratique, N=1000", color = "black")

ax.set_ylim(200,1050)

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "pointb_N10_N1000_T(r)", ext = ext)






#-------------------(i) heat(r)-------------------
#read data uniform case
outputs_Ns_uniform_heat = [f"./outputs/N={Ns[i]}_uniform_heat.out" for i in range(nsimul_Ns)]

#N10
N10_data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[0])
N10_heat_r_num_uniform = N10_data_uniform_heat[:,1]
N10_r_heat_uniform = N10_data_uniform_heat[:,0]

#N1000
N1000_data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[-1])
N1000_heat_r_num_uniform = N1000_data_uniform_heat[:,1]
N1000_r_heat_uniform = N1000_data_uniform_heat[:,0]


#read data non uniform case
outputs_Ns_nonuniform_heat = [f"./outputs/N={Ns[i]}_nonuniform_heat.out" for i in range(nsimul_Ns)]

#N10
N10_data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[0])
N10_heat_r_num_nonuniform = N10_data_nonuniform_heat[:,1]
N10_r_heat_nonuniform = N10_data_nonuniform_heat[:,0]

#N1000
N1000_data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[-1])
N1000_heat_r_num_nonuniform = N1000_data_nonuniform_heat[:,1]
N1000_r_heat_nonuniform = N1000_data_nonuniform_heat[:,0]



#plot
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$r$ [m]", ylabel=r'$j_Q$ [W/m$^2$]')

ax.plot(N10_r_heat_uniform, N10_heat_r_num_uniform, marker="+", label = "Équidistant, N=10", color = "navy")
ax.plot(N10_r_heat_nonuniform, N10_heat_r_num_nonuniform, marker="+", label = "Quadratique, N=10", color = "red")
ax.plot(N1000_r_heat_nonuniform, N1000_heat_r_num_nonuniform, linestyle = "--", label = "Quadratique, N=1000", color = "black")


plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=1)
u.savefig(fig, "pointb_N10_N1000_heat(r)", ext = ext)



#-------------------(i) power balance-------------------
exp_part = np.exp(-r0**2/(sigma**2)) - np.exp(-(R-r0)**2/(sigma**2))
erf_part = np.sqrt(np.pi)*r0/(sigma) * (m.erf((R-r0)/sigma) - m.erf(-r0/sigma))
P_tot = np.pi*S0*sigma**2 * (exp_part + erf_part)

#N10
boundary_heat_flux_uniform = 2*np.pi*R*N10_heat_r_num_uniform[-1]
boundary_heat_flux_nonuniform = 2*np.pi*R*N10_heat_r_num_nonuniform[-1]

print("N10 : ")
print("P_tot = ", P_tot)
print("boundary_heat_flux_uniform = ", boundary_heat_flux_uniform)
print("boundary_heat_flux_nonuniform = ", boundary_heat_flux_nonuniform)
print("relative error uniform = ",100*abs(P_tot - boundary_heat_flux_uniform)/P_tot)
print("relative error nonuniform = ",100*abs(P_tot - boundary_heat_flux_nonuniform)/P_tot)


#N1000
boundary_heat_flux_uniform = 2*np.pi*R*N1000_heat_r_num_uniform[-1]
boundary_heat_flux_nonuniform = 2*np.pi*R*N1000_heat_r_num_nonuniform[-1]

print("N1000 : ")
print("P_tot = ", P_tot)
print("boundary_heat_flux_uniform = ", boundary_heat_flux_uniform)
print("boundary_heat_flux_nonuniform = ", boundary_heat_flux_nonuniform)
print("relative error uniform = ",100*abs(P_tot - boundary_heat_flux_uniform)/P_tot)
print("relative error nonuniform = ",100*abs(P_tot - boundary_heat_flux_nonuniform)/P_tot)





#-------------------(ii) convergence study-------------------

lasts_T0_uniform = np.zeros(nsimul_Ns)
lasts_heat_R_uniform = np.zeros(nsimul_Ns)
errors_boundary_heat_R_uniform = np.zeros(nsimul_Ns)

lasts_T0_nonuniform = np.zeros(nsimul_Ns)
lasts_heat_R_nonuniform = np.zeros(nsimul_Ns)
errors_boundary_heat_R_nonuniform = np.zeros(nsimul_Ns)

first_hs = np.zeros(nsimul_Ns)

#Iterating throught
for i in range(nsimul_Ns):
    data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[i])
    data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[i])
    data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[i])
    data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[i])
    
    
    boundary_heat_R_th = P_tot
    
    #uniform case
    T0 = data_uniform_T[0,1]
    heat_R = data_uniform_heat[-1,1]
    boundary_heat_R = 2*np.pi*R*heat_R
    
        #errors
    last_T0 = T0 
    last_heat_R = heat_R
    error_boundary_heat_R = abs(boundary_heat_R - boundary_heat_R_th)

    lasts_T0_uniform[i] = last_T0
    lasts_heat_R_uniform[i] = last_heat_R
    errors_boundary_heat_R_uniform[i] = error_boundary_heat_R
    
    #nonuniform case
    T0 = data_nonuniform_T[0,1]
    heat_R = data_nonuniform_heat[-1,1]
    boundary_heat_R = 2*np.pi*R*heat_R
    
    first_hs[i] = data_nonuniform_T[1,0]- data_nonuniform_T[0,0]
        
        #errors 
    last_T0 = T0
    last_heat_R = heat_R
    error_boundary_heat_R = abs(boundary_heat_R - boundary_heat_R_th)
    
    lasts_T0_nonuniform[i] = last_T0
    lasts_heat_R_nonuniform[i] = last_heat_R
    errors_boundary_heat_R_nonuniform[i] = error_boundary_heat_R
    
#plot error T0,uniforme
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N^2}$", ylabel=r'$T(0)$ [K]', grid_bool = False)

ax.plot(1/Ns[:]**2,lasts_T0_uniform[:], marker="+", label = "Équidistant", color = "navy")

# #polyfit of order 4 on the data
# coeffs = np.polyfit(1/Ns[:],lasts_T0_uniform[:],8)
# x_fit= np.linspace(0,0.1,100)
# ax.plot(x_fit,np.polyval(coeffs,x_fit),linestyle = "--", color = "green",label = "Fit sur équidistant")

# ax.grid(which='both', linestyle='--')
# plt.tight_layout()
# u.set_legend_properties(ax, fontsize=18,ncol=2)
# u.savefig(fig, "pointb_convergence_T0_uniforme", ext = ext)


# #plot error T0,non uniforme
# ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'$T(0)$ [K]', grid_bool = False)

ax.plot(1/Ns[:]**2,lasts_T0_nonuniform[:], marker="x", label = "Quadratique", color = "red")

# #polyfit of order 4 on the data
# coeffs = np.polyfit(1/Ns[:],lasts_T0_nonuniform[:],8)
# x_fit= np.linspace(0,0.1,100)
# ax.plot(x_fit,np.polyval(coeffs,x_fit),linestyle = "--", color = "orange",label = "Fit sur quadratique")


ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_T0", ext = ext)



#error on heat_R
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'$j_Q(R)$ [W/m$^2$]', grid_bool = False)

ax.plot(1/Ns,lasts_heat_R_uniform, marker="+", label = "Équidistant", color = "navy")

# #fit
# coeffs = np.polyfit(1/Ns,lasts_heat_R_uniform,8)
# x_fit= np.linspace(0,0.1,100)
# ax.plot(x_fit,np.polyval(coeffs,x_fit),linestyle = "--", color = "green",label = "Fit sur équidistant")


ax.plot(1/Ns,lasts_heat_R_nonuniform, marker="x", label = "Quadratique", color = "red")

# #fit
# coeffs = np.polyfit(1/Ns,lasts_heat_R_nonuniform,8)
# x_fit= np.linspace(0,0.1,100)
# ax.plot(x_fit,np.polyval(coeffs,x_fit),linestyle = "--", color = "orange",label = "Fit sur quadratique")


ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_heat_R", ext = ext)


#error on boundary_heat_R
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'Erreur sur $P_{\mathrm{tot}}$ [W]', grid_bool = False)

ax.loglog(1/Ns,errors_boundary_heat_R_uniform, marker="+", label = "Équidistant", color = "navy")  
ax.loglog(1/Ns,errors_boundary_heat_R_nonuniform, marker="x", label = "Quadratique", color = "red")

x = np.linspace(1e-4,1e-2,100)
ax.loglog(x,2e-1*x,linestyle = "--", color = "black")
ax.text(9e-4,1e-4,r"~ $\frac{1}{N}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_boundary_heat_R", ext = ext)



#error on T0 with first_Hs
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$h_1^4$ [m]", ylabel=r'$T(0)$ [K]', grid_bool = False)

ax.plot(first_hs[:]**4,lasts_T0_nonuniform[:], marker="x", label = "Quadratique", color = "red")

# #polyfit of order 4 on the data
# coeffs = np.polyfit(first_hs[:],lasts_T0_nonuniform[:],8)
# x_fit= np.linspace(0,0.0175,100)
# ax.plot(x_fit,np.polyval(coeffs,x_fit),linestyle = "--", color = "black",label = "Fit")


ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_T0_h1", ext = ext)




#------------------------------------(ii) valeurs convergées------------------------------------
coeffs = np.polyfit(1/Ns[:],lasts_T0_uniform[:],6)
T0_uniform_fit = np.polyval(coeffs,0)
print("Valeur convergée de T(0) pour équidistant : ",np.polyval(coeffs,0))

coeffs = np.polyfit(1/Ns[:],lasts_T0_nonuniform[:],6)
T0_nonuniform_fit = np.polyval(coeffs,0)
print("Valeur convergée de T(0) pour quadratique : ",np.polyval(coeffs,0))

coeffs = np.polyfit(1/Ns,lasts_heat_R_nonuniform,6)
heat_nonuniform_fit = np.polyval(coeffs,0)
print("Valeur convergée de j_Q(R) pour quadratique : ",np.polyval(coeffs,0))

coeffs = np.polyfit(1/Ns,lasts_heat_R_uniform,6)
heat_uniform_fit = np.polyval(coeffs,0)
print("Valeur convergée de j_Q(R) pour équidistant : ",np.polyval(coeffs,0))


#------------------------------------(ii) Erreur de convergence par rapport au fit------------------------------------

errors_T0_uniform = np.zeros(nsimul_Ns)
errors_heat_R_uniform = np.zeros(nsimul_Ns)

errors_T0_nonuniform = np.zeros(nsimul_Ns)
errors_heat_R_nonuniform = np.zeros(nsimul_Ns)

first_hs = np.zeros(nsimul_Ns)

#Iterating throught
for i in range(nsimul_Ns):
    data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[i])
    data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[i])
    data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[i])
    data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[i])
    
    
    #uniform case
    T0_th = T0_uniform_fit
    heat_R_th = heat_uniform_fit   
    
    T0 = data_uniform_T[0,1]
    heat_R = data_uniform_heat[-1,1]
    
        #errors
    error_T0 = abs(T0 - T0_th)
    error_heat_R = abs(heat_R - heat_R_th)

    errors_T0_uniform[i] = error_T0
    errors_heat_R_uniform[i] = error_heat_R
    
    #nonuniform case
    T0_th = T0_nonuniform_fit
    heat_R_th = heat_nonuniform_fit
    
    T0 = data_nonuniform_T[0,1]
    heat_R = data_nonuniform_heat[-1,1]
    
    first_hs[i] = data_nonuniform_T[1,0]- data_nonuniform_T[0,0]
        
        #errors 
    error_T0 = abs(T0 - T0_th)
    error_heat_R = abs(heat_R - heat_R_th)
    
    errors_T0_nonuniform[i] = error_T0
    errors_heat_R_nonuniform[i] = error_heat_R
    
    
#plot error T0
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{1}{N}$", ylabel=r'$|T(0)-T(0)_{fit}|$ [K]', grid_bool = False)

ax.loglog(1/Ns[:],errors_T0_uniform[:], marker="+", label = "Équidistant", color = "navy")

ax.loglog(1/Ns[:],errors_T0_nonuniform[:], marker="+", label = "Quadratique", color = "red")

x = np.linspace(1e-4,1e-2,100)
ax.loglog(x,2e4*x**2,linestyle = "--", color = "black")
ax.text(1.4e-3,2e-2,r"~$\frac{1}{N^2}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_error_T0", ext = ext)


#plot error T0 with first_hs
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$h_1$ [m]", ylabel=r'$|T(0)-T(0)_{fit}|$ [K]', grid_bool = False)

ax.loglog(first_hs[:],errors_T0_nonuniform[:], marker="+", label = "Quadratique", color = "red")

x = np.linspace(5e-4,5e-3,100)
ax.loglog(x,2e9*x**4,linestyle = "--", color = "black")
ax.text(2e-3,1e-2,r"~$\frac{1}{N^4}$", fontsize = 18)

ax.grid(which='both', linestyle='--')
plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "pointb_convergence_error_T0_h1", ext = ext)





# xlim=[0,0.005]
# ylim=[970,980]
# x_zoom=0.03
# y_zoom=0.03
# X=[1/Ns[:-1],1/Ns[:-1]]
# Y=[lasts_T0_uniform[:-1],lasts_T0_nonuniform[:-1]]
# colors=["navy","red"]
# markers=["+","+"]
# corner1 = ["top","top"]
# corner2 = ["low","low"]

# u.zoom_in_plot(ax,X,Y,xlim,ylim,x_zoom,y_zoom,colors,markers,["-", "-"],corner1,corner2,7.8)



# xlim=[0,0.00002]
# ylim=[973,978]
# x_zoom=0.03
# y_zoom=0.03
# X=[first_hs[:-1]**2]
# Y=[lasts_T0_nonuniform[:-1]]
# colors=["red"]
# markers=["+"]
# corner1 = ["top","top"]
# corner2 = ["low","low"]

# u.zoom_in_plot(ax,X,Y,xlim,ylim,x_zoom,y_zoom,colors,markers,["-", "-"],corner1,corner2,7.8)

