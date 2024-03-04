import numpy as np
import subprocess
import matplotlib.pyplot as plt
import utils_v2 as u

# Parameters
# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice1_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file

prefix = ["Ng0_Implicite_","Ng0_Semi_","Ng0_Explicite_"]
ext = "png"

nsteps1 = np.arange(20,40,1)
nsteps2 = np.arange(40,100,5)
nsteps3 = np.array([42,44,46,48,52,58])
nsteps4 = np.arange(100,200,20)
nsteps = np.concatenate((nsteps1,nsteps2,nsteps3,nsteps4),axis = 0)
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 21.642813194695652

dt = tfin / nsteps



pi = np.pi
g = 9.81  # m/s^2

# Analysis
# TODO insert the values
m = 0.056  

v0 = 0
omega   = 62.83185307179586
mu = 6
R = 0.033
rho = 1.2

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
vy_list = []
xy_list = []
for j in range(nsimul):
    output_file = f"{paramstr}={param[j]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[j]:.15g} output={output_file} alpha=0.5"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')
    
    

error = np.zeros(nsimul)
error_mec = np.zeros(nsimul)
datas = []

xfinal = np.zeros(nsimul)

for j in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[j])  # Load the output file of the i-th simulation
    t = data[:, 0]
    
    datas.append(data)
    

    xx = data[-1, 1]  # final position, velocity, energy
    yy = data[-1, 2]
    vx = data[-1, 3]
    vy = data[-1, 4]
    En = data[-1, 5]

    vy_list.append(vy)
    if j == 7 :
        xy_list.append(xx)
        xy_list.append(yy)
    xfinal[j] = xx


lw = 1.5 # line width. TODO: adjust if needed
fs = 16  # font size. TODO: adjust if needed
ms = 10 # marker size. TODO: adjust if needed



ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"x[m]", ylabel=r'y[m]', grid_bool=False)

N = 7
ax.plot(datas[N][:, 1], datas[N][:, 2], label = r"$N_{steps}$" + f" = 2000")

ax.scatter(xy_list[0],xy_list[1],label = "Final position",color = "black",marker = "o",s = 300)

ext = "png"
plt.tight_layout()
plt.grid()
plt.axis('equal')
u.save_pdf(fig, "Ng0_trainee_Position" + "." + ext)


print(vy_list)


norder = 2
ax,fig = u.create_figure_and_apply_format(figsize=(10, 6),xlabel=r"$dt^{2}$[$m^2$]", ylabel=r'$x[m]$', grid_bool=False)

ax.plot(dt**norder,xfinal,linestyle = "-", marker = "+",color = "blue",markersize = 15) 
plt.grid()
u.save_png(fig, "Ng0_trainee_convergence" + "." + ext)



print("x final position = ",xy_list[0])
print("y final position = ",xy_list[1])
print("tfin = ",tfin)

