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
executable = 'Exercice5_students'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'input_example'  # Name of the input file

g = 9.81

hL =7000
hC = 35
hR = 200
xL = 0
xa = 300e3
xb = 700e3
xc = 720e3
xd = 850e3
xR = 1000e3

A = 1
x1 = 50e3
x2 = 250e3

tfin = 12000

#se propage de gauche à droite
initial_state = "right"
#pas mode propre
initialization = "pas mode"

cb_gauche = "sortie"
cb_droite = "sortie"

#output tous les 
n_stride = 500

#umax initial
umax = np.sqrt(hL*g)

#nombre de point maillage
N = 100000
dx = (xR - xL) / (N-1)


#beta = u * dt / dx
#on veut un dt tel que max(beta) = 1
dt = dx / umax
nsteps = int(tfin / dt)


cmd = f"{repertoire}{executable} {input_filename} N={ Ns[i]:.15g} output={output_file} IsUniform=true"

print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')

