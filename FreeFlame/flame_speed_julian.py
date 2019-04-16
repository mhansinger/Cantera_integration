# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:47:47 2017
Compare Mechanisms Methane/Oxygen at 20bar
Flame Speed
@author: julian
"""

import cantera as ct
import numpy as np
import csv
from matplotlib.pylab import *

import os

#filename = os.path.basename(__file__[0:-3])

#User Input
mechanism = 'gri30.cti' #'2S_CH4_BFER' #Lu19 #valorani # 1S_CH4_MP1 # NAJM # 2S_CH4_BFER # F1D_Jones_KT1 # 2S_CH4_CM2 #Lu13 #Lu30
p = 101325.0  # pressure [Pa]
Tin = 300  # unburned gas temperature [K]

# ----------------------------------------------------------------#
mechanism_file = mechanism
gas = ct.Solution(mechanism_file)
fuel_species = 'CH4'
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

# air composition
air_N2_O2_ratio = 3.76
stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')

# equivalence ratio range
phi = 0.9
sl = 0

initial_grid = [0.0,0.0001,0.001,0.01, 0.02,0.029,0.03]
initial_grid = 2*array([0.0, 0.001, 0.005, 0.01, 0.0149, 0.015],'d')/3 #[x / 4.0 for x in initial_grid]
tol_ss = [1.0e-5, 1.0e-8]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-8]  # [rtol atol] for time stepping

'''
##############################################################################

X = np.zeros(gas.n_species)
X[ifuel] = phi
X[io2] = stoich_O2
X[in2] = stoich_O2*air_N2_O2_ratio
# set the gas state
gas.TPX = Tin, p, X

'''
##############################################################################
## calculation with user provided mass fractions
Y = np.zeros(gas.n_species)
Y[ifuel] = 0.04977
Y[io2] = 0.2204
Y[in2] = 0.72983

# set the gas state
gas.TPY = Tin, p, Y

##############################################################################


# Create flame
gas()
f = ct.FreeFlame(gas, initial_grid)
#f.restore('tmp_Slavinskaya.xml')
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)
# Set boundaries
X=gas.X
f.inlet.X = X
f.inlet.T = Tin
f.P = p

# Solve with the energy equation disabled
print('Solving without energy equation. \n')

f.energy_enabled = False
f.set_max_jac_age(50, 50)
f.set_time_step(1e-6, [2, 5, 10, 20, 80, 120, 200])
f.set_refine_criteria(ratio=7.0, slope=0.1, curve=0.1)
f.solve(loglevel=5, refine_grid=True)
#f.set_refine_criteria(ratio = 3.0, slope = 0.8, curve = 0.8)
#f.solve(loglevel=1, refine_grid=True)

#With energy equation
f.energy_enabled = True
f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
f.solve(loglevel=5, refine_grid=True)
f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
f.solve(loglevel=1, refine_grid=True)
f.set_refine_criteria(ratio = 5.0, slope = 0.3, curve = 0.3)
f.solve(loglevel=1, refine_grid=True)
f.set_refine_criteria(ratio = 3.0, slope = 0.1, curve = 0.1)
f.solve(loglevel=1, refine_grid=True)
f.set_refine_criteria(ratio = 2.0, slope = 0.05, curve = 0.05, prune = 0.01)
f.solve(loglevel=1, refine_grid=True)
sl = f.u[0]
f.save(mechanism + '.xml')
print('\n mixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))
#print "\n Temperature = ", f.T

#################################################################
## Calculate thermal flame thickness
#################################################################

dTdx = np.zeros(f.flame.n_points-1)

for n in range(f.flame.n_points-1):
    dTdx[n] = (f.T[n+1]-f.T[n])/(f.flame.grid[n+1]-f.flame.grid[n])

slopeDT = max(dTdx)
thermalThickness = (max(f.T)-min(f.T))/slopeDT

#print('\nthermal Flamethickness = {:10f} m\n'.format(thermalThickness))

#################################################################

plt.plot(f.grid,f.T)
plt.show()
