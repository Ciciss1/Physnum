# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:22:14 2024

@author: Administrator
"""

import numpy as np
import subprocess
import matplotlib.pyplot as plt

# Parameters
# TODO adapt to what you need (folder path executable input filename)

repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'Ex2_2024_student'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file
ext = "pdf"


thetadot = np.array([1.0,2.0,5.0,9.0,10.5,11.0,11.025,11.3,11.5,11.6,11.8,11.9,12.0,16.5,17.5,17.6,17.9,20,22,-17.5,-17.6,-17.9,-20]) 
# thetadot = np.array([-17.5,-17.6,-17.9,-20,-22])
nsimul   = len(thetadot)        # Number of simulations to perform

nsteps   = 2000 #change according what you need


paramstr  = 'thetadot0'  # Parameter name to scan
paramstr2 = 'nsteps'  # Parameter name to scan
paramstr3 = 'sampling'  # Parameter name to scan
param     = thetadot  # Parameter values to scan

print(thetadot)

# Simulations
outputs = []  # List to store output file names
vy_list = []
for i in range(nsimul):
    output_file = f"./outputs/{paramstr}={param[i]}.out"
    outputs.append(output_file)

# for i in range(nsimul):  # Iterate through the simulations
#     output_file = outputs[i]
#     cmd = f"{repertoire}{executable} {input_filename} {paramstr3}={nsteps:.15g} {paramstr2}={nsteps:.15g} {paramstr}={param[i]:.15g} N={10000} output={output_file}"

#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')

error = np.zeros(nsimul)
plt.figure(figsize=(10, 6))

for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t = data[:, 0]

    xx = data[:, 1] 
    yy = data[:, 2]   
    plt.plot(np.mod(xx+np.pi,2*np.pi)-np.pi, yy, '*', markersize=1, label = f'$\dot{{\\theta}}$ = {thetadot[i]}')  
    
plt.xlabel(r'$\theta$', fontsize=20)  # LaTeX for theta
plt.ylabel(r'$\dot{\theta}$', fontsize=20)  # Lating the y-axis label with fontsize 20

plt.xticks(fontsize=15)  
plt.yticks(fontsize=15)  


rect_xlim = [-2.3,-1.2]
rect_ylim= [-5,8]
#add a rectangle on the plot, with a little text pointing to the rectangle with "1"
plt.gca().add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))
plt.text(rect_xlim[1]-0.2, rect_ylim[0]+0.7, '1', fontsize=15, color='black', zorder = 15)

rect_xlim = [-np.pi,-2.5]
rect_ylim= [8.5,12.5]
#add a rectangle on the plot, with a little text pointing to the rectangle with "2"
plt.gca().add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))
plt.text(rect_xlim[1]-0.2, rect_ylim[0]+0.7, '2', fontsize=15, color='black', zorder = 15)

rect_xlim = [2,np.pi]
rect_ylim= [-15,-8.5]
#add a rectangle on the plot, with a little text pointing to the rectangle with "3"
plt.gca().add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))
plt.text(rect_xlim[1]-0.2, rect_ylim[0]+0.7, '3', fontsize=15, color='black', zorder = 15)

rect_xlim = [-0.5,0.2]
rect_ylim= [15,19.5]
#add a rectangle on the plot, with a little text pointing to the rectangle with "4"
plt.gca().add_patch(plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1]-rect_xlim[0], rect_ylim[1]-rect_ylim[0], fill=None, edgecolor='black', lw=2, zorder = 15))
plt.text(rect_xlim[1]-0.2, rect_ylim[0]+0.7, '4', fontsize=15, color='black', zorder = 15)


plt.grid(True)  
# plt.legend()
plt.savefig(f"./png/poincare" + f".{ext}")  # Save the figure to a file


plt.figure(figsize=(10, 8))
for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t = data[:, 0]

    xx = data[:, 1] 
    yy = data[:, 2]   
    plt.plot(np.mod(xx+np.pi,2*np.pi)-np.pi, yy, '*', markersize=1, label = f'$\dot{{\\theta}}$ = {thetadot[i]}')  
    
plt.xlabel(r'$\theta$', fontsize=20)  # LaTeX for theta
plt.ylabel(r'$\dot{\theta}$', fontsize=20)  # Lating the y-axis label with fontsize 20
plt.xticks(fontsize=15)  
plt.yticks(fontsize=15)  

plt.grid(True)  
plt.xlim(-2.3,-1.2)
ylim_down = -5
ylim_up = 8
plt.ylim(ylim_down,ylim_up)

plt.savefig(f"./png/poincare_zoom1" + f".{ext}")  # Save the figure to a file



plt.xlim(-np.pi,-2.5)
ylim_down = 8.5
ylim_up = 12.5
plt.ylim(ylim_down,ylim_up)

plt.savefig(f"./png/poincare_zoom2" + f".{ext}")  # Save the figure to a file



plt.xlim(2,np.pi)
ylim_down = -15
ylim_up = -8.5
plt.ylim(ylim_down,ylim_up)

plt.savefig(f"./png/poincare_zoom3" + f".{ext}")  # Save the figure to a file




# plt.xlim(1.2,2.2)
# ylim_down = -8
# ylim_up = 5.5
# plt.ylim(ylim_down,ylim_up)

# plt.savefig(f"./png/poincare_zoom4" + f".{ext}")  # Save the figure to a file


# plt.xlim(-1,1)
# ylim_down = 9.8
# ylim_up = 13
# plt.ylim(ylim_down,ylim_up)

# plt.savefig(f"./png/poincare_1_05" + f".{ext}")  # Save the figure to a file


plt.xlim(-0.5,0.3)
ylim_down = 15
ylim_up = 20
plt.ylim(ylim_down,ylim_up)

plt.savefig(f"./png/poincare_zoom4" + f".{ext}")  # Save the figure to a file

