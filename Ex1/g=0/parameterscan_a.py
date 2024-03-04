import numpy as np
import subprocess
import matplotlib.pyplot as plt
import utils_v2 as u

# Parameters
# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice1_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file

prefix = ["g0_Implicite_","g0_Semi_","g0_Explicite_"]
ext = "pdf"

nsteps =np.array([20,35,50,75,100,150,200,300,400,500,750,900,1000,1250,1500,1750,2000,5000]) # TODO change
toplot =[[50,2],[200,6],[1000,12],[2000,16],[5000,17]] # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 60 # TODO: Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps



pi = np.pi
g = 9.81  # m/s^2

# Analysis
# TODO insert the values
m = 0.056  

v0 = 5
omega   = 62.83185307179586
mu = 6
R = 0.033
rho = 1.2

gamma = mu * (R**3) * rho * omega/m

E_th = 0.5 * m * v0**2 + 0.2 * m * (R**2) * (omega**2)

T = 2*pi/gamma


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
param_alpha = [0,0.5,1] #implicite, semi-implicite, explicite
noms = ["Implicite","Semi-implicite","Explicite"]
colors = ["navy","darkgreen","magenta"]

# Simulations
outputs = [[],[],[]]  # List to store output file names
vy_list = []
for i in range(len(param_alpha)) : 
    for j in range(nsimul):
        output_file = f"{param_alpha[i]}_{paramstr}={param[j]}.out"
        outputs[i].append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[j]:.15g} output={output_file} alpha={param_alpha[i]}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')
    
    

error = [np.zeros(nsimul),np.zeros(nsimul),np.zeros(nsimul)]
error_mec = [np.zeros(nsimul),np.zeros(nsimul),np.zeros(nsimul)]
datas = [[],[],[]]



for i in range(len(param_alpha)) :
    for j in range(nsimul):  # Iterate through the results of all simulations
        data = np.loadtxt(outputs[i][j])  # Load the output file of the i-th simulation
        t = data[:, 0]
        
        datas[i].append(data)
        
        xx = data[-1, 1]  # final position, velocity, energy
        yy = data[-1, 2]
        vx = data[-1, 3]
        vy = data[-1, 4]
        En = data[-1, 5]
        
        En_min = np.min(data[:, 5])
        En_max = np.max(data[:, 5])
        
        # TODO compute the error for each simulation
        error[i][j] = np.sqrt((x_th-xx)**2 + (y_th-yy)**2)
        error_mec[i][j] = En_max-En_min


lw = 1.5 # line width. TODO: adjust if needed
fs = 20  # font size. TODO: adjust if needed
ms= 10 # marker size. TODO: adjust if needed

 
for i in range(len(param_alpha)) :
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_xlabel('x [m]', fontsize=fs)
    ax.set_ylabel('y [m]', fontsize=fs)

    for j in range(len(toplot)) : 
        if i == 2 and j == 0 : 
            continue
        ax.plot(datas[i][toplot[j][1]][:, 1], datas[i][toplot[j][1]][:, 2], label = r"$N_{steps}$" + f" = {toplot[j][0]}")

    t = np.linspace(0,T,1000)
    x = (v0/gamma)*np.sin(gamma*t)
    y = (v0/gamma)*(1- np.cos(gamma*t))

    ax.plot(x, y, label = "Solution analytique", color = "black", linestyle = "--")
    
   
        
    plt.axis('equal')
    if i == 0 : 
        ax.set_xlim(-40,20)
        ax.set_ylim(-10,40)
    if i == 2 : 
        ax.set_xlim(-55,40)
        ax.set_ylim(-20,60)
    plt.tight_layout()
    u.set_legend_properties(ax,colors=["black","navy"], fontsize = 12, loc = "best", markers=["+","-"],ncol = 1)
    
    plt.grid()
    plt.savefig("./" + ext + "/" + prefix[i]+'Positions' + "." + ext)



    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    ax.set_xlabel('t [s]', fontsize=fs)
    ax.set_yscale('log')
    ax.grid(which='both')
    for j in range(len(toplot)) : 
        ax.plot(datas[i][toplot[j][1]][:, 0], datas[i][toplot[j][1]][:, 5], label = r"$N_{steps}$" + f" = {toplot[j][0]}")
    ax.axhline(y=E_th, color='black', linestyle='--', label = r"$E_{th}$")
    u.set_legend_properties(ax,colors=["black","navy"], fontsize = 12, loc = "best", markers=["+","-"],ncol = 1) 
    plt.tight_layout()
    plt.savefig("./" + ext + "/" + prefix[i]+'Energies' + "." + ext)






ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_x$ [m]', grid_bool=False)
print("Erreur position : ")
for i in range(len(param_alpha)) :
    ax.loglog(dt, error[i], color = colors[i], marker = "+", linestyle = "-", markersize = 10, label = noms[i])

    # a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error[i]), color = colors[i], ax = ax, precisions = [4,4], plot = False)
    
    # ax.plot(x, (x**a)*np.exp(b), color = colors[i]) 
    # print(f"alpha = {param_alpha[i]} : a = {a:.3f} et b = {np.exp(b):.3f}")

x = np.linspace(np.min(dt),np.max(dt),1000)
ax.plot(x, 45.0*x, color = "black", label = r"~ $\Delta_t$", linestyle = "--")
ax.plot(x, 1*x**2, color = "red", label = r"~ $\Delta_t^2$", linestyle = "--")

ax.legend()    
u.set_legend_properties(ax,colors=["black","navy"], fontsize = fs, loc = "best",markers=["+","-"],ncol = 1)
plt.grid(True, which="both", linestyle='--')
u.save_pdf(fig, "g0_ALL_Erreurs_position" + "." + ext)





ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$\Delta$ t [s]", ylabel=r'$\delta_E$ [J]', grid_bool=False)
print("Erreur énergie : ")
for i in range(len(param_alpha)) :
    ax.loglog(dt, error_mec[i], color = colors[i], marker = "+", linestyle = "-", markersize = 10, label = noms[i])

    # if i != 1 : 
    #     a,a_error, b, b_error = u.linear_fit(np.log(dt),np.log(error_mec[i]), color = colors[i], ax = ax, precisions = [4,4], plot = False)
    #     x = np.linspace(0.03,13,1000)
    #     ax.plot(x, (x**a)*np.exp(b), color = colors[i]) 
    #     print(f"alpha = {param_alpha[i]} : a = {a:.3f} et b = {np.exp(b):.3f}")
     
x = np.linspace(np.min(dt),np.max(dt),1000)
ax.plot(x, 3.5*x, color = "black", label = r"~ $\Delta_t$", linestyle = "--")

ax.legend()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = fs, loc = "best",markers=["+","-"],ncol = 1)

plt.grid(True, which="both", linestyle='--')

u.save_pdf(fig, "g0_ALL_Erreurs_energie" + "." + ext)



#Nsteps = 75 !!!!!!!
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"t[s]", ylabel=r'E[J]', grid_bool=False)
for i in range(len(param_alpha)) :
    ax.plot(datas[i][3][:, 0], datas[i][3][:, 5], label = noms[i])
    
plt.tight_layout()
plt.legend()
ax.set_yscale('log')
plt.grid()
u.set_legend_properties(ax,colors=["black","navy"], fontsize = fs, loc = "best",markers=["+","-"],ncol = 1)
u.save_pdf(fig, "g0_ALL_Energie" + "." + ext)


