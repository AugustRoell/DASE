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

from datetime import datetime
import os, sys
sys.path.append('../')

from DASE_Libraries import DASE, profiles, config

now = datetime.now()
current = now.strftime("%Y-%m-%d-%H-%M-%S")
data_folder = f"../Data/sim_{current}"
os.makedirs(data_folder)

########################
#      Simulation      #
######################## 

# simple example, simulating a spherical cap and saving 
# figures of the potential and the first few wave functions 

H, U_array = DASE.Hamiltonian(profiles.sphere)
DASE.PotentialSaveFig(U_array, path=data_folder+"/potential.png")

eta_list, normalized_psi_list = DASE.SolveOSE(H)
for sol in range(config.Nsol):
    DASE.WaveFunctionSaveFig(normalized_psi_list, path=data_folder+f"/sol_{sol}", sol=sol)