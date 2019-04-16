"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

from: CANTERA
"""

import cantera as ct
import numpy as np

mechanism = 'lu19.cti' #'gri30.cti'

out_filename = 'FreeFlame_lean_states.h5'

# Operating conditions
p = ct.one_atm  # pressure [Pa]
Tin = 300  # inlet temperature unburnt gas [K]
yF1 = 1  # fuel
yO2 = 0.232  # oxidizer

# Parameters#
width = 0.005  # m
nf = 300  # number of flames to be computed
nu_st = 4  # stoichiometric
LFL = 0.029  # lower flammability limit CH4 mass fraction
UFL = 0.055  # upper flammability limit CH4 mass fraction
#yFuel = np.linspace(LFL, UFL, nf)  # Define YFuel within flammability limits

yF = UFL

# # Tolerances, not used...
# tol_ss = [1.0e-9, 1.0e-10]  # [rtol atol] for steady-state problem
# tol_ts = [1.0e-5, 1.0e-10]  # [rtol atol] for time stepping
loglevel = 0 #3  # amount of diagnostic output (0-5)
refine_grid = True  # True to enable refinement, False to

yOx = (1 - yF) * 0.232
yN2 = (1 - yF) * 0.768
yH = 1e-10

# reactants in mass fraction
reactants = ('CH4:%f O2:%f N2:%f H:%f' % (yF, yOx, yN2, yH))

# Calculate Bilger Mixture fraction
yC = yF * 0.748753
yH = yF * 0.251247
yO = yOx
f_Bilger = (2 * yC / 12.01 + yH / (2 * 1.01) - (yO - 0.232) / 16) / \
           (2 * (0.75 / 12.01) + 0.25 / (2 * 1.01) + 0.232 / 16)

print('Bilger Mixture Fraction is: %f' % f_Bilger)

width = 0.01  # m
initial_grid = np.linspace(0, width, 150)

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution(mechanism)
gas.TPY = Tin, p, reactants

# Set up flame object
f = ct.FreeFlame(gas, initial_grid)
f.set_refine_criteria(ratio=3, slope=0.3, curve=0.3)
f.show_solution()

# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)

f.write_csv('lu19_freeFlame.csv', species='Y',quiet=False)


