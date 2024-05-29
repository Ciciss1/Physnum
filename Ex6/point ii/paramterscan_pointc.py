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
tfin = 0.03
xL = -1
xR = 1
n_v = 2
x0 = -0.5
n = 16
sigma_norm = 0.04

V0 = 1700
Nintervals = 512

# dt_array = np.logspace(-5, -3, num=15, base=10, dtype=float)
dt_array = np.array([1e-5, 2.5e-5,5e-5, 1e-4,2e-4,2.5e-4,3.333333333333e-4,5e-4,6.6666666666666e-4,7.5e-4,1e-3])
Ncases = len(dt_array)

#creating the outputs 
outputs = []
for i in range(Ncases) : 
    outputs.append(f"dt_{dt_array[i]}")
    
# run the simulations
for i in range(Ncases):
    cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt_array[i]} Nintervals={Nintervals} output={outputs[i]}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


prob_right_003_array = np.zeros(Ncases)
for i in range(Ncases):
    #extract data
    data_obs = np.loadtxt(f"./outputs/{outputs[i]}_obs.out")
    data_pot = np.loadtxt(f"./outputs/{outputs[i]}_pot.out")
    data_psi2 = np.loadtxt(f"./outputs/{outputs[i]}_psi2.out")

    t,_,prob_right,_,_,_,_,_,_,_ = data_obs.T

    #we want to get the probability prob_right at t=0.03, so we need to find the index of the time closest to 0.03
    print("t", t[-1])
    prob_right_003_array[i] = prob_right[-1]
    
#-----------------convergence study-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\Delta t^2$[a.u.]", ylabel=r'$P_{x>0}(t=0.03)$')

ax.plot(dt_array**2, prob_right_003_array, label=r"$P_{x>0}(t=0.03)$", marker = "+", color = "navy")

#linear fit 
def f(x, a, b):
    return a*x + b
popt, pcov = curve_fit(f, dt_array**2, prob_right_003_array)

#plot the fit
# x = np.linspace(0, dt_array[-1], 100)
# ax.plot(x, f(x, popt[0],popt[1]+1e-5), label=r"~ $\Delta t^2$", color = "red", linestyle = "--")

plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, f"POINT_C_convergence_dt", ext = ext)

    



