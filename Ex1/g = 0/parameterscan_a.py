import numpy as np
import subprocess
import matplotlib.pyplot as plt
import utils_v2 as u

# Parameters
# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice1_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file

prefix = "g0_Semi_"
ext = "pdf"

nsteps =np.array([5,10,20,50,100,200,500]) # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 60 # TODO: Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps

pi = np.pi
g = 9.81  # m/s^2

# Analysis
# TODO insert the values
m = 0.056  

v0 = 5
omega = 10
mu = 6
R = 0.033
rho = 1.2

E_th = 0.5 * m * v0**2 + 0.2 * m * (R**2) * (omega**2)

gamma = mu * (R**3) * rho * omega/m
# add the other variables
# TODO: Insert here the expressions for the exact final solution
x_th  = (v0/gamma)*np.sin(gamma*tfin)
y_th  = (v0/gamma)*(1- np.cos(gamma*tfin))
vx_th = v0*np.cos(gamma*tfin)
vy_th = v0*np.sin(gamma*tfin)
"""
... and other parameters
"""

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
vy_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')
    
    


error = np.zeros(nsimul)
error_mec = np.zeros(nsimul)
datas = []

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]
    
    datas.append(data)
    

    xx = data[-1, 1]  # final position, velocity, energy
    yy = data[-1, 2]
    vx = data[-1, 3]
    vy = data[-1, 4]
    En = data[-1, 5]
    
    En_min = np.min(data[:, 5])
    En_max = np.max(data[:, 5])
    
    vy_list.append(vy)
    # TODO compute the error for each simulation
    error[i] = np.sqrt((x_th-xx)**2 + (y_th-yy)**2)
    error_mec[i] = En_max-En_min

lw = 1.5 # line width. TODO: adjust if needed
fs = 16  # font size. TODO: adjust if needed
ms= 10 # marker size. TODO: adjust if needed

 
fig, ax = plt.subplots(constrained_layout=True)
ax.set_xlabel('x [m]', fontsize=fs)
ax.set_ylabel('y [m]', fontsize=fs)
for i in range(len(datas)) : 
    ax.plot(datas[i][:, 1], datas[i][:, 2], label = r"$N_{steps}$" + f" = {nsteps[i]}")
    
t = np.linspace(0,tfin,1000)
x = (v0/gamma)*np.sin(gamma*t)
y = (v0/gamma)*(1- np.cos(gamma*t))

ax.plot(x, y, label = "Solution analytique", color = "black", linestyle = "--")
   
plt.tight_layout()
plt.axis('equal')
plt.legend()
plt.grid()
plt.savefig("./" + ext + "/" + prefix+'Positions' + "." + ext)

fig, ax = plt.subplots(constrained_layout=True)
ax.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
ax.set_xlabel('t [s]', fontsize=fs)
ax.set_yscale('log')
ax.grid(which='both')
for i in range(len(datas)) : 
    ax.plot(datas[i][:, 0], datas[i][:, 5], label = r"$N_{steps}$" + f" = {nsteps[i]}")
ax.axhline(y=E_th, color='black', linestyle='--', label = r"$E_{th}$")
plt.legend()    
plt.tight_layout()
plt.savefig("./" + ext + "/" + prefix+'Energies' + "." + ext)



# uncomment the following 2 lines if you want debug
#import pdb
#pbd.set_trace()
# plt.figure()
# plt.loglog(dt, error, '+-', linewidth=lw, color = "navy", markersize=ms)
# plt.xlabel(r'$\Delta$ t [s]', fontsize=fs)
# plt.ylabel(r'$\delta_x$ [m]', fontsize=fs)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig('Erreurs_position.png')

ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_x$ [J]')
ax.scatter(dt, error, label = "Erreurs position", color = "black", marker = "+", s=100)
ax.set_xscale('log')
ax.set_yscale('log')
a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error), color = "black", ax = ax, precisions = [4,4], plot = False)
x = np.linspace(0.1,13,1000)
ax.plot(x, (x**a)*np.exp(b), color = "navy", label = rf"$y = e^{{{b:.3f}}}x^{{{a:.3f}}}$") 
ax.legend()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "upper left",markers=["+","-"],ncol = 1)
u.save_pdf(fig, prefix + "Erreurs_position" + "." + ext)



ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_E$ [J]')
ax.scatter(dt, error_mec, label = "Erreurs énergie", color = "black", marker = "+", s=100)
ax.set_xscale('log')
ax.set_yscale('log')
# a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_mec), color = "red", ax = ax, precisions = [4,4], plot = False)
# x = np.linspace(0.1,13,1000)
# ax.plot(x, (x**a)*np.exp(b), color = "red", label = rf"$y = e^{{{b:.3f}}}x^{{{a:.3f}}}$") 
ax.legend()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "upper left",markers=["+","-"],ncol = 1)
u.save_pdf(fig,prefix + "Erreurs_energie" + "." + ext)


# plt.figure()
# plt.loglog(dt, error_mec, '+-', linewidth=lw, color = "navy", markersize=ms)
# plt.xlabel(r"$\Delta$ t [s]", fontsize=fs)
# plt.ylabel(r'$\delta_E$ [J]', fontsize=fs)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig('Erreurs_energie.png')

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, TODO: modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1  # TODO: Modify if needed

plt.figure()
plt.plot(dt**norder, vy_list, 'k+-', linewidth=lw)
plt.xlabel('\Delta t [s]', fontsize=fs)
plt.ylabel('v_y [m/s]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot3.png')

plt.show()
