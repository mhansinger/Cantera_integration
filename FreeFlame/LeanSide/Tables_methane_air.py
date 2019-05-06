# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 17:19:24 2015
Premixed freely propagating flame
Adiabatic
Table generation TNF workshop (Methane-Air, 1 bar)
@author: julian & hagen
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from write_full_csv import write_full_csv

# Operating conditions
p = 1e5 # operating pressure [Pa]
Tin = 300 # inlet temperature unburnt gas [K]
yF1 = 1 # fuel
yO2 = 0.232 # oxidizer

# Parameters#
nf = 151 # number of flames to be computed
nu_st = 4 # stoichiometric 
LFL = 0.025 # lower flammability limit CH4 mass fraction
UFL = 0.055 # upper flammability limit CH4 mass fraction
yFuel = np.linspace(UFL, LFL, nf) # Define YFuel within flammability limits

# Initial grid
#initial_grid = [0.0,0.0001, 0.001, 0.01, 0.015, 0.02, 0.029,0.02999, 0.03]
#initial_grid = [x / 3 for x in initial_grid]
initial_grid = np.linspace(0,0.05, 200)
tol_ss = [1.0e-9, 1.0e-14]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-14]  # [rtol atol] for time stepping
loglevel  = 5               # amount of diagnostic output (0-5)
refine_grid = True          # True to enable refinement, False to

# Loop over all mixture fractions within flammability limits
for index, yF in enumerate(yFuel):
    yOx = (1 - yF)*0.232
    yN2 = (1 - yF)*0.768
    reactants = ('CH4:%f O2:%f N2:%f' % (yF, yOx, yN2))
    print('Reactants are %s' %reactants)
    mixtureFraction = (nu_st*yF-yOx+yO2)/(nu_st*yF1+yO2)
    # Calculate Bilger Mixture fraction
    yC = yF*0.748753
    yH = yF*0.251247
    yO = yOx
    mixtureFractionBilger = (2*yC/12.01 + yH/(2*1.01) - (yO - 0.232)/16)/(2*(0.75/12.01) + 0.25/(2*1.01) + 0.232/16)
    print('Mixture fraction = %f' %mixtureFraction)
    print('Mixture fraction Bilger = %f' %mixtureFractionBilger)
    
    # IdealGasMix object used to compute mixture properties
    gas = ct.Solution('gri30.xml', 'gri30_mix')
    gas.TPY = Tin, p, reactants
    
    # Flame object
    f = ct.FreeFlame(gas, initial_grid)
    f.restore('solution.xml')
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.inlet.Y = reactants
    f.inlet.T = Tin
    
    # Solve with the energy equation disabled
    print('Solving without energy equation. \n')
    f.energy_enabled = False
    f.set_max_jac_age(50, 50)
    f.set_time_step(1e-5, [2, 5, 10, 80])
    f.set_refine_criteria(ratio=3, slope=1, curve=1)
    #f.set_initial_guess()
    f.solve(loglevel=1, refine_grid=True)
    
    # Solve with the energy equation enabled
    print('Switching on energy equation. \n')
    f.energy_enabled = True
    #f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    f.set_refine_criteria(ratio=3, slope=1, curve=1)
    f.solve(loglevel, refine_grid)
    f.set_refine_criteria(ratio = 5.0, slope = 0.3, curve = 0.3)
    f.solve(loglevel, refine_grid)
    f.set_refine_criteria(ratio = 3.0, slope = 0.1, curve = 0.1)
    f.solve(loglevel, refine_grid)
    
    f.save('solution.xml')
    f.save('Tabellen/solution_yF%s_f%s.xml' %(str(yF), str(mixtureFractionBilger)))
    
    print 'mixture averaged flamespeed = ',f.u[0]
    write_full_csv(f,'Tabellen/p%s_freeflame_f%s.csv' %(p, str(mixtureFractionBilger)),species='Y', quiet=False) 
