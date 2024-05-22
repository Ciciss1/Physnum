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
Nintervals = 256

output = "V0_0"
    
# run the simulations
cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} xL={xL} xR={xR} V0={V0} nv={n_v} x0={x0} n={n} sigma_norm={sigma_norm} dt={dt} Nintervals={Nintervals} output={output}"

print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')
