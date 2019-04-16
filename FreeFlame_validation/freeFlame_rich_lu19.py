
"""
Premixed freely propagating flame
Adiabatic
Table generation TNF workshop (Methane-Air, 1 bar)
@author: julian & hagen, mod: hansinger

last change: 11.3.2019
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from write_full_csv import write_full_csv
import os
from os.path import join


# storage path
store_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/FreeFlame/FreeFlame_raw.h5'

mechanism_file =  'gri30.cti' #'utils/lu19.cti'#


# Operating conditions
p = 1e5  # operating pressure [Pa]
Tin = 300  # inlet temperature unburnt gas [K]
yF1 = 1  # fuel
yO2 = 0.232  # oxidizer

# Parameters#
nf = 1#751  # number of flames to be computed
nu_st = 4  # stoichiometric
LFL = 0.055  # lower flammability limit CH4 mass fraction
UFL = 0.205  # upper flammability limit CH4 mass fraction
yFuel = np.linspace(LFL, UFL, nf)  # Define YFuel within flammability limits

# Initial grid
# initial_grid = [0.0,0.0001, 0.001, 0.01, 0.015, 0.02, 0.029,0.02999, 0.03]
# initial_grid = [x / 3 for x in initial_grid]
initial_grid = np.linspace(0, 0.02, 200)
tol_ss = [1.0e-9, 1.0e-10]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-10]  # [rtol atol] for time stepping
loglevel = 5  # amount of diagnostic output (0-5)
refine_grid = True  # True to enable refinement, False to



# Loop over all mixture fractions within flammability limits
for index, yF in enumerate(yFuel):
    yOx = (1 - yF) * 0.232
    yN2 = (1 - yF) * 0.768
    reactants = ('CH4:%f O2:%f N2:%f' % (yF, yOx, yN2))
    print('Reactants are %s' % reactants)
    mixtureFraction = (nu_st * yF - yOx + yO2) / (nu_st * yF1 + yO2)
    # Calculate Bilger Mixture fraction
    yC = yF * 0.748753
    yH = yF * 0.251247
    yO = yOx
    mixtureFractionBilger = (2 * yC / 12.01 + yH / (2 * 1.01) - (yO - 0.232) / 16) / (
                2 * (0.75 / 12.01) + 0.25 / (2 * 1.01) + 0.232 / 16)
    #print('Mixture fraction = %f' % mixtureFraction)
    print('Mixture fraction Bilger = %f' % mixtureFractionBilger)

    # IdealGasMix object used to compute mixture properties
    gas = ct.Solution(mechanism_file)#('utils/lu19.cti')
    gas.TPY = Tin, p, reactants

    # Flame object
    f = ct.FreeFlame(gas, initial_grid)
    try:
        f.restore('solution.xml')
    except:
        print('No solution.xml')
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.inlet.Y = reactants
    f.inlet.T = Tin

    # Solve with the energy equation disabled
    print('Solving without energy equation. \n')
    f.energy_enabled = False
    f.set_max_jac_age(150, 150)
    f.set_time_step(1e-6, [2, 5, 10, 80])
    f.set_refine_criteria(ratio=3, slope=1, curve=1)
    # f.set_initial_guess()
    f.solve(loglevel=1, refine_grid=True)

    # Solve with the energy equation enabled
    print('Switching on energy equation. \n')
    f.energy_enabled = True
    # f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    f.set_refine_criteria(ratio=3, slope=1, curve=1)
    f.solve(loglevel, refine_grid)
  #  f.set_refine_criteria(ratio=5.0, slope=0.3, curve=0.3)
  #  f.solve(loglevel, refine_grid)
    f.set_refine_criteria(ratio=3.0, slope=0.5, curve=0.5)
    f.solve(loglevel, refine_grid)

    f.write_csv('GRI30.csv',species='Y')

    # hdf_database = pd.HDFStore(store_path)
    #
    # #update the hdf5 database
    # hdf_database.append('TNF_raw_data',data_df)
    # hdf_database.close()
