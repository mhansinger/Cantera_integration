'''
Laminar flamespeed for a CH4/Air mixture at ambient condition as a function of ROF.
@author: Paola
Date: Nov-2018 modified for RAMEC
'''

import cantera as ct
import numpy as np
from matplotlib.pylab import *
import sys
import csv

##############################################################################
# Select input conditions
T = 300
P = 1.01325e05

######### Import gas phases from the chemical mechanism ######################
# Detailed Mechanisms: gri30, valorani, NAJM
# Analitically Reduced Mechanisms: Lu13(xml!), Lu19, Lu30
# Global mechanisms: F1D_JONES_KT1, 2S_CH4_BFER2, F1D_JONES_KT1, 2S_CH4_CM2, 1S_CH4_MP1
mechanism='utils/lu19.cti'

######## Select the mixture transport model according to the .cti file in use ###############
#transport='gri30_mix'
#transport='NAJM'
#transport= 'CH4_CM2_mix'
#transport= 'CH4_CM2_mix'
#transport= 'CH4_MP1_mix'
#transport='CH4_BFER_mix'
#transport='CH4_JONES_KT1_mix'
#transport='CH4_JONES_KT2_mix'
#transport='gas' # valoriani, Lu19.cti, Lu30.cti and and Lu13.xml

mechanism_file=mechanism #Change to .xml for Lu13, the rest is .cti
gas = ct.Solution(mechanism_file)
#gas.transport_model = 'gas'

# find fuel, nitrogen, and oxygen indexes
fuel_species = 'CH4'
m = gas.n_species
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

# air composition
air_N2_O2_ratio = 3.76
stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')

#Set tolerance properties for FreeFlame
tol_ss    = [1.0e-5, 1.0e-9]        # [rtol atol] for steady-state problem (-8)
tol_ts    = [1.0e-5, 1.0e-9]        # [rtol atol] for time stepping
loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)
refine_grid = True

#Set initial grid
initial_grid = 2*array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/3

##############################################################################
# create some arrays to hold the data
sl = 0
xeq = np.zeros(m)
yeq = np.zeros(m)

'''
##############################################################################
phi = 0.9

#set the stoichiometry
xeq[ifuel] = phi
xeq[io2] = stoich_O2
xeq[in2] = stoich_O2*air_N2_O2_ratio

# set the gas state
gas.TPX = T, P, xeq

'''
##############################################################################
## calculation with user provided mass fractions
Y = np.zeros(gas.n_species)
Y[ifuel] = 0.04977
Y[io2] = 0.2204
Y[in2] = 0.72983

# set the gas state
gas.TPY = T, P, Y


##############################################################################

xeq = gas.X
yeq = gas.Y

gas()

    #Create the free laminar premixed flame
f = ct.FreeFlame(gas,initial_grid)
    #Flame BC
f.inlet.X = xeq
f.inlet.T = T
    #Flame tolerances
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

    #First flame:
f.energy_enabled = False
f.set_refine_criteria(ratio = 10.0, slope = 1, curve = 1)
f.set_max_jac_age(50, 50)
f.set_time_step(1e-7, [2, 5, 10, 80])
f.solve(loglevel, refine_grid = False)


    #Second flame and so on..
f.energy_enabled = True
f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
f.solve(loglevel, refine_grid=True)


#################### Only for NAJM ####################### -> time convergence problem, changed grid refinement
#f.set_refine_criteria(ratio = 5.0, slope = 0.2, curve = 0.2) #parameters for NAJM
#f.solve(loglevel, refine_grid=True)
#sl[i]=f.u[0]
#print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))


##################### For the rest of the mechanisms #########
f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05)
f.solve(loglevel, refine_grid=True)
    #Not decreasing the ratio to 1.0 to avoid exceeding grid points #3000!
f.set_refine_criteria(ratio = 2.0, slope = 0.02, curve = 0.02)
f.solve(loglevel, refine_grid=True)
f.set_refine_criteria(ratio = 2.0, slope = 0.01, curve = 0.01)
f.solve(loglevel, refine_grid=True)
    #Solved enough flames to have the speed converging according to graph CERFACS Cantera
    #Save the flamespeed for step i
sl=f.u[0]
    #Print the sl on screen at the end of the loop
print('\nLewis 1 flamespeed = {:7f} m/s\n'.format(f.u[0]))

#################################################################
## Calculate thermal flame thickness
#################################################################

dTdx = np.zeros(f.flame.n_points-1)

for n in range(f.flame.n_points-1):
    dTdx[n] = (f.T[n+1]-f.T[n])/(f.flame.grid[n+1]-f.flame.grid[n])

slopeDT = max(dTdx)
thermalThickness = (max(f.T)-min(f.T))/slopeDT

# print('\nthermal Flamethickness = {:10f} m\n'.format(thermalThickness))

plt.plot(f.grid,f.T)
plt.show()

#################################################################

## write output CSV file for importing into Excel
#csv_file = 'flamespeed_'+ mechanism +'.csv'
#with open(csv_file, 'w') as outfile:
#    writer = csv.writer(outfile)
#    writer.writerow(['rof','sl', 'phi'])
#    for i in range(npoints):
#        writer.writerow([rof[i], sl[i], phi[i]])
#print "Output written to", "%s"%(csv_file)



