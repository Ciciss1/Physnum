import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Ex2_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([100,200,500,1000,2000,5000,10000]) # TODO change as you need
nsimul = len(nsteps)  # Number of simulations to perform


dt = 1/nsteps # define it as you need


theta_0 = 1e-10
# Analysis
# TODO: Insert here the expressions for the exact final solution
final_theta_th = 1e-10
final_thetadot_th = 0

l = 0.200000000000
g = 9.8100000000000
omega_0 = np.sqrt(g/l)
T = 2*np.pi/omega_0

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  

for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
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
    
plt.figure(figsize=(10, 6))
plt.loglog(dt, final_theta_list) # this is an example, change it if you need
plt.xlabel("$dt$", fontsize=20)  # LaTeX for theta
plt.ylabel(r'${\theta}$', fontsize=20)  # Lating the y-axis label with fontsize 20
plt.xticks(fontsize=15)  
plt.yticks(fontsize=15)  
plt.grid(True)

# TODO add also the other plots that you need

print('error_theta', error_theta)
print(final_theta_list)
  
plt.savefig("error_theta.png")  # Save the figure to a file  
    
