# -*- coding: utf-8 -*-
"""Created on Tue Nov 07 2023.

@author: August

This file contains examples on how to run DASE simulations
and save the partial results. The relevant simulation and 
material parameters can be set in the config.py file, and 
the profile function parameters can be altered in profiles.py
"""

#########################################
#      Import & Directory Creation      #
######################################### 

from numpy import (
    linspace,
    Infinity,
    meshgrid,
    array,
)

from matplotlib.pyplot import (
    subplots,
    grid,
    colorbar,
    show,
    savefig,
)

from csv import writer
from datetime import datetime
import os, sys
sys.path.append('../')

from DASE_Libraries import DASE, profiles, config, materials

now = datetime.now()
current = now.strftime("%Y-%m-%d-%H-%M-%S")
data_folder = f"../Data/SMQ_{current}"
os.makedirs(data_folder)

########################
#      Simulation      #
######################## 

# example of a parameter sweep for the SMQ values of a design

SMQ_list = []
cutoffradius_range = linspace(1.0e-6, 3.0e-3, 10)       # the range of dome cutoff radii
r_outer_range = linspace(180e-6, 600e-6, 10)            # the range of outer forcing ring radii
SMQ_opt = -Infinity

for cutoffradius in cutoffradius_range:
    profile_lamb = lambda r : profiles.cylinder(r, h=0.434e-6, r_cyl=cutoffradius)  # here you should set the profile of interest
    H, _ = DASE.Hamiltonian(profile_lamb, l=0, Dirichlet=True)
    eta_list, normalized_psi_list = DASE.SolveOSE(H)
    frequencies = eta_list*DASE.f_0
    SMQ_list_r = []
    for r in r_outer_range:
        overlaps = DASE.EM_ForcingOverlaps(normalized_psi_list, r_outer=r)
        SMQ_r = DASE.SMQ(overlaps, frequencies, g=250e+3)
        if SMQ_r > SMQ_opt:
            SMQ_opt = SMQ_r
            r_opt = r
            cutoffradius_opt = cutoffradius
        SMQ_list_r.append(SMQ_r)
    SMQ_list.append(SMQ_list_r)
opt_parameters = (SMQ_opt, r_opt, cutoffradius_opt)

with open(data_folder+"/SMQ_Map_values.csv", "w", newline='') as file:      # saving a file with SMQ values for further analysis
    write = writer(file) 
    write.writerows(SMQ_list)

fig,ax = subplots()
grid(which='major', axis='both')
X, Y = meshgrid(r_outer_range*1e6, cutoffradius_range*1e6)
Z = array(SMQ_list)
surf = ax.contourf(X, Y, Z, 100, cmap='coolwarm')
ax.plot(opt_parameters[1] * 1e6, opt_parameters[2] * 1e6, marker='o', markersize=6, color='k', 
        label=f'Optimum: SMQ = {opt_parameters[0]*100:0.3f}%, \n'
        f'r$_{{outer}}$ = {opt_parameters[1]*1e6:0.2f}$\mu$m, '
        f'r$_{{cutoff}}$ = {opt_parameters[2]*1e6:0.2f}$\mu$m')
fig.set_size_inches(8.7, 7)
ax.set_title('SMQ Map');
ax.set_xlabel('$r_{outer}$ $(\mu m)$')
ax.set_ylabel('Cutoff Radius $(\mu m)$')
ax.grid(visible=False)
ax.legend()
colorbar(surf, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
savefig(fname=data_folder+"/SMQ_Map.png")
show()                                                                      # Also showing the final figure to allow for further inspection
