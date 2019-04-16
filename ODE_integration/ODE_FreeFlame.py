# ODE integration of FreeFlame data


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pickle
from os.path import join
import Ofpp     # reads OF data to numpy array
import h5py
from utils.clean_states import clean_states
import os

import cantera as ct

####################################
# time step to integrate over
dt = 1e-6




###################################


# path to FreeFlame data
# storage path
files_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/FreeFlame'
out_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/FreeFlame/ODE_integrated'

# species order is the same as in gas object. DO NOT CHANGE!!
species_lu19 = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'CH3', 'CH4', 'CO', 'CO2', 'CH2O', 'CH3OH', 'C2H2', 'C2H4', 'C2H6', 'CH2CO', 'N2']
thermo_fields = ['T','p','f_Bilger']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

files = ['FreeFlame_lean_states.h5','FreeFlame_rich_states.h5']

mechanism_file = 'utils/lu19.cti'

# construct output df
species_lu19_after = [f+'_after' for f in species_lu19]
species_lu19_difference = ['d'+f for f in species_lu19]
species_lu19_RR = ['RR_'+f for f in species_lu19]   # Reaction rate -> change / sec

out_columns = species_lu19 + thermo_fields + species_lu19_after + species_lu19_difference + species_lu19_RR
out_columns.append('T_after')
out_columns.append('dT')
out_columns.append('dt')

ODE_out_df = pd.DataFrame(columns=out_columns)

# create outpath directory if not available
if not os.path.isdir(out_path):
    os.mkdir(out_path)



# loop over lean and rich states
for f in files:
    this_path = join(files_path,f)
    print(this_path)
    # read in the thermo-chemical states as df

    raw_df=pd.read_hdf(this_path)
    mass_fraction_before_np = np.zeros((len(species_lu19)))
    gas = ct.Solution(mechanism_file)

    this_df = raw_df.iloc[455]

    for name in species_lu19:
        index = gas.species_index(name)
        mass_fraction_before_np[index] = this_df[name]


    gas.Y = mass_fraction_before_np
    print(gas.Y)
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])
    sim.advance(dt)

    print(r.thermo)
    #
    # # loop over rows in raw data
    # for i in range(raw_df.shape[0]):
    #     this_raw_df = raw_df.iloc[i]
    #
    #     # set up cantera reactor
    #     gas = ct.Solution(mechanism_file)
    #     T_before = this_raw_df['T']
    #
    #     if T_before > 800:
    #         # set the chemical state
    #         gas.TP = T_before, ct.one_atm
    #         #f_Bilger = raw_df['f_Bilger']
    #
            # mass_fraction_before_np = np.zeros((len(species_lu19)))
            # for name in species_lu19:
            #     index = gas.species_index(name)
            #     mass_fraction_before_np[index] = this_raw_df[name]
    #
    #         gas.Y = mass_fraction_before_np
    #
    #         try:
    #             assert gas.species_names == species_lu19
    #         except AssertionError:
    #             print('species order is not correct! Needs to be identical as in gas.species!')
    #             break
    #
    #         # create reactor
    #         r = ct.IdealGasConstPressureReactor(gas)
    #         sim = ct.ReactorNet([r])
    #         sim.advance(dt)
    #
    #         mass_fraction_after_np = r.thermo.Y
    #         Y_diff = mass_fraction_after_np - mass_fraction_before_np
    #         Y_RR = Y_diff / dt  # production rate of species mass fraction per second!
    #
    #         T_after = r.thermo.T
    #         diff_T = T_after - T_before
    #
    #         # print(T_before)
    #         # print()
    #         #
    #         this_out_df = pd.DataFrame(columns=out_columns)
    #         this_out_df[species_lu19] = mass_fraction_before_np
    #         this_out_df['T'] = T_before
    #         this_out_df['f_Bilger'] = f_Bilger
    #         this_out_df['p'] = ct.one_atm
    #         this_out_df['dt'] = dt
    #         this_out_df['dT'] = diff_T
    #         this_out_df['T_after'] = T_after
    #         this_out_df[species_lu19_after] = species_lu19_after
    #         this_out_df[species_lu19_difference] = Y_diff
    #         this_out_df[species_lu19_RR] = Y_RR
    #        # print(f_Bilger)
    #         print(Y_RR[gas.species_index('CH4')])
    #
    #
    #
    #
    #
    #










