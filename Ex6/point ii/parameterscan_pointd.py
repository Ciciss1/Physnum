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
Nintervals = 256

x = np.linspace(xL, xR, Nintervals+1)

V0s1 = np.linspace(300,2500,30)
V0s2 = np.linspace(2500,5000,10)
V0s = np.concatenate((V0s1,V0s2))
NV0s = len(V0s)

#create outputs
outputs=[]
for i in range(NV0s):
    outputs.append(f"pointd_V0={V0s[i]}")
    

E_V0 = np.zeros(NV0s)
Ptrans_array = np.zeros(NV0s)
for i in range(NV0s):
    V0 = V0s[i]
    output = outputs[i]

    # # run the simulation
    # cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

    # print(cmd)
    # subprocess.run(cmd, shell=True)
    # print('Done.')


    #extract data
    data_obs = np.loadtxt(f"./outputs/{output}_obs.out")
    data_pot = np.loadtxt(f"./outputs/{output}_pot.out")
    data_psi2 = np.loadtxt(f"./outputs/{output}_psi2.out")

    t,prob_left,prob_right,E,xmoy,x2moy,pmoy,p2moy,xincertitude,pincertitude = data_obs.T
    
    E_V0[i] = np.mean(E)/V0
    
    #prob at t=0.03
    index = np.argmin(np.abs(t-0.03))
    P_trans = prob_right[index]
    Ptrans_array[i] = P_trans
    
#-----------------plot Ptrans(E_V0)-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(8, 6),xlabel=r"$\frac{<E>}{V_0}$", ylabel=r'$P_{trans}$')

ax.plot(E_V0, Ptrans_array, marker = "+", color = "navy",label=r"$P_{trans}(E/V0)$")

#get the probabilitiy at E/V0 = 1, considering linear greps between all the points
i = 0
while E_V0[i] > 1:
    i += 1

#the grep has the form y = a*x + b, with a = (y2 - y1)/(x2 - x1) and b = y1 - a*x1
#we want to find the value of y at x = 1
#hence, a = (y2 - y1)/(x2 - x1) = (Ptrans_array[i] - Ptrans_array[i-1])/(E_V0[i] - E_V0[i-1])
#and b = y1 - a*x1 = Ptrans_array[i-1] - a*E_V0[i-1]
a  = (Ptrans_array[i-1] - Ptrans_array[i])/(E_V0[i-1] - E_V0[i])
b = Ptrans_array[i] - a*E_V0[i]

P_trans_1 = a*1 + b

ax.plot(1, P_trans_1, marker = "o", color = "red",label = r"$P_{trans}$ at $E=V_0$ :  " + f"{P_trans_1:.3f}")

ax.vlines(1, 0, P_trans_1, color = "red", linestyle = "--")
ax.hlines(P_trans_1, 0, 1, color = "red", linestyle = "--")

print("Probabilité avec E=V0 : ",P_trans_1)


plt.tight_layout()
u.set_legend_properties(ax, fontsize=18)
u.savefig(fig, "POINT_D_Ptrans(E_V0)", ext = ext)


