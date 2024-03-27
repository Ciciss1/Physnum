import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import utils_v2 as u

from scipy.signal import find_peaks
from scipy.optimize import curve_fit



# TODO adapt to what you need (folder path executable input filename)

# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Ex2_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


nsteps = 20000

pi=3.1415926535897932384626433832795028841971e0

dt = 1/nsteps # define it as you need


theta_0 = [0,1e-6]
thetadot_0 = [15.0,15.0]

nsimul = len(theta_0) 

L = 0.2
g = 9.81
alpha = 0
d = 0.01

omega_0 = np.sqrt(g/L)
omega = 2*omega_0

T_excitation = 2*pi/omega

paramstr = 'theta0'  # Parameter name to scan
param = theta_0  # Parameter values to scan


# Simulations
outputs = []  

#periode d excitation
N = 40
tfin = N*T_excitation

print(tfin)

for i in range(nsimul):
    output_file = f"./outputs/{paramstr}={param[i]}_Nsteps={nsteps}.out"
    outputs.append(output_file)
    
# for i in range(nsimul):
#     output_file = outputs[i]
#     cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file} N={N}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')
    

data1 = np.loadtxt(outputs[0])
data2 = np.loadtxt(outputs[1])

theta1 = data1[:, 1]
thetadot1 = data1[:, 2]
theta2 = data2[:, 1]
thetadot2 = data2[:, 2]

t = data1[:, 0]

ext = "pdf"


#----------plot phase diagramm 1------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta_1$ [rad]", ylabel=r'$\dot{\theta}_1$ [rad/s]')

ax.plot(theta1,thetadot1,linestyle = "-", color = "navy") 

plt.tight_layout()
u.savefig(fig, "ExcitationLyapounov_PortraitPhase_theta1",ext)


#----------plot phase diagramm 2------------
ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$\theta_2$ [rad]", ylabel=r'$\dot{\theta}_2$ [rad/s]')

ax.plot(theta2,thetadot2,linestyle = "-", color = "navy") 

plt.tight_layout()
u.savefig(fig, "ExcitationLyapounov_PortraitPhase_theta2",ext)



#compute delta(t) in phase space
delta = np.sqrt(omega_0*(theta1-theta2)**2 + (thetadot1-thetadot2)**2)

ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\delta$ [s$^{-1}$]')

ax.plot(t,delta,linestyle = "-", color = "navy", label = r"$\delta(t)$")


def model(x, a, b):
    return b * np.exp(a * x)

params, covariance = curve_fit(model, t, delta)

# a,a_error,b,b_error = u.linear_fit(t,np.log(delta),"red",precisions=[3,10],ax = ax, plot = False)

# ax.plot(t,np.exp(a*t+b),linestyle = "--", color = "red", label = fr"fit $\delta(t)$ : $y(x) = ({a:.3f} \pm {a_error:.3f})x - ({abs(b):.3f} \pm {abs(b_error):.3f})$")
ax.plot(t,model(t,*params),linestyle = "--", color = "red", label = fr"fit $\delta(t)$ : $y(x) = {params[1]:.3f} \exp({params[0]:.3f}x)$")


u.set_legend_properties(ax,colors=["black","navy"], fontsize = 20, loc = "best",markers=["+","-"],ncol = 1)
plt.tight_layout()
u.savefig(fig, "ExcitationLyapounov_delta(t)",ext)





ax,fig = u.create_figure_and_apply_format(figsize=(10, 8),xlabel=r"$t$ [s]", ylabel=r'$\delta$ [s$^{-1}$]')

ax.plot(t,delta,linestyle = "-", color = "navy", label = r"$\delta(t)$")

ax.set_yscale("log")

plt.tight_layout()
u.savefig(fig, "ExcitationLyapounov_delta(t)_semilog",ext)
