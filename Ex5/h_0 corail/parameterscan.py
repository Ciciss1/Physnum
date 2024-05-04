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
    
    
#Iterating throught
for i in range(nsimul_Ns):
    data_uniform_T = np.loadtxt(outputs_Ns_uniform_T[i])
    data_uniform_heat = np.loadtxt(outputs_Ns_uniform_heat[i])
    data_nonuniform_T = np.loadtxt(outputs_Ns_nonuniform_T[i])
    data_nonuniform_heat = np.loadtxt(outputs_Ns_nonuniform_heat[i])
    
