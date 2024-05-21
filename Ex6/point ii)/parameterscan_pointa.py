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
executable = 'Exercice6_2024_students'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'input_example'  # Name of the input file


#for input file
tfin = 0.1
xL = -1
xR = 1
V0 = 1
n_v = 2
x0 = -0.5
n = 16
sigma_norm = 0.04

dt = 1e-4
Nintervals = 256
output = "pointa"

# run the simulation

cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')



#extract
data = np.loadtxt(output + "_f")
t = data[:, 0]

print("SHAPE DATA : ",data.shape)


#-----------------plot la vague à différents t-----------------
ax,fig = u.create_figure_and_apply_format(figsize=(12, 6),xlabel=r"$x$ [m]", ylabel=r'Amplitude [m]')

x = np.linspace(xL, xR, N)

for i in [0,200,300,600,1100,1500]:
    ax.plot(x,data[i,1:], label=f"t = {t[i]:.0f} s")
        
# u.add_ticks("x",[xa,xb,xc,xc],[r'$x_a$',r'$x_b$',r'$x_c$',r'$x_d$'],ax=ax,divide =5)
# ax.vlines((xc+xd)/2,-1,4,linestyles="--",color="black",label=r"$Pic du corail$")
ax.vlines(xa,-1,4,linestyles="-.",color="red",label=r"$x_a$")
ax.vlines(xb,-1,4,linestyles="--",color="blue",label=r"$x_b$")
ax.vlines(xc,-1,4,linestyles="-.",color="purple",label=r"$x_c$")
ax.vlines(xd,-1,4,linestyles="--",color="black",label=r"$x_d$")


plt.tight_layout()
u.set_legend_properties(ax, fontsize=18,ncol=2)
u.savefig(fig, "Vague_difft", ext = ext)



