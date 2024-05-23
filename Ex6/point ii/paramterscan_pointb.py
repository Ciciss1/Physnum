import sys
# Add the directory containing the module to sys.path
sys.path.append('/workspaces/Physnum')

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import math as m
import utils_v2 as u

from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation
from matplotlib import animation

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

ext = "pdf"


# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Exercice6_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


#for input file
tfin = 0.1
xL = -1
xR = 1
n_v = 2
x0 = -0.5
n = 16
sigma_norm = 0.04

dt = 1e-4

V0 = 1
Nintervals_array = np.array([250,300,350,450, 550, 650,750, 850,950])
Ncases = len(Nintervals_array)

#creating the outputs 
outputs = []
for i in range(Ncases) : 
    outputs.append(f"Nintervals_{Nintervals_array[i]}")
    
# # run the simulations
# for i in range(Ncases):
#     cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals_array[i]} output={outputs[i]}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')


prob_right_003_array = np.zeros(Ncases)
for i in range(Ncases):
    #extract data
    data_obs = np.loadtxt(f"./outputs/{outputs[i]}_obs.out")
    data_pot = np.loadtxt(f"./outputs/{outputs[i]}_pot.out")
    data_psi2 = np.loadtxt(f"./outputs/{outputs[i]}_psi2.out")

    t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T

    
    #extract avery 3 colomns
    psi_module = data_psi2[:,::3]

    #we want to get the probability prob_right at t=0.03, knowing that t starts at zero, gos to 0.1, and as dt of 1e-4,
    prob_right_003 = prob_right[301]

    prob_right_003_array[i] = prob_right_003
    
#-----------------convergence study-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$N_{intervals}$", ylabel=r'$P_{x>0}(t=0.03)$')

ax.plot(1/Nintervals_array, prob_right_003_array, label=r"$P_{x>0}(t=0.03)$")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"POINT_B_convergence", ext = ext)

    



