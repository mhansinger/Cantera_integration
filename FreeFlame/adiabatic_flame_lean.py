"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

from: CANTERA
"""

import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join

# Simulation parameters
# storage path
store_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/FreeFlame'

mechanism = 'utils/lu19.cti' #'gri30.cti'

out_filename = 'FreeFlame_lean_states.h5'

species_lu19 = ['C2H2','C2H4', 'C2H6', 'CH2CO', 'CH2O', 'CH3', 'CH3OH', 'CH4', 'CO', 'CO2', 'H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH']
thermo_fields = ['T','p','f_Bilger']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# Operating conditions
p = ct.one_atm  # pressure [Pa]
Tin = 600  # inlet temperature unburnt gas [K]
yF1 = 1  # fuel
yO2 = 0.232  # oxidizer

# Parameters#
width = 0.005  # m
nf = 300  # number of flames to be computed
nu_st = 4  # stoichiometric
LFL = 0.029  # lower flammability limit CH4 mass fraction
UFL = 0.0555  # upper flammability limit CH4 mass fraction
yFuel = np.linspace(LFL, UFL, nf)  # Define YFuel within flammability limits


# # Tolerances, not used...
# tol_ss = [1.0e-9, 1.0e-10]  # [rtol atol] for steady-state problem
# tol_ts = [1.0e-5, 1.0e-10]  # [rtol atol] for time stepping
loglevel = 0 #3  # amount of diagnostic output (0-5)
refine_grid = True  # True to enable refinement, False to


for yF in yFuel: #[0.055,0.06]:
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
    f.set_refine_criteria(ratio=3, slope=0.9, curve=0.9)
    f.show_solution()

    # Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)

    # store as data frame
    this_df = pd.DataFrame(columns=all_fields)

    # loop over species fields
    for s in species_lu19:
        s_index = f.gas.species_index(s)
        this_df[s] = f.Y[s_index,:]

    one_vector = f.T/f.T
    this_df['T'] = f.T
    this_df['p'] = one_vector * f.P
    this_df['f_Bilger'] = one_vector * f_Bilger

    # save the data to the database
    hdf_database = pd.HDFStore(join(store_path, out_filename))
    # update the hdf5 database
    hdf_database.append('FreeFlame_lean_data', this_df)
    hdf_database.close()
    print('Database updated')

    # plt.close('all')
    # plt.plot(f.grid, f.T)
    # plt.title('f_Bilger: %.3f' % f_Bilger)
    # plt.show(block=False)



