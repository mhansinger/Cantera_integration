# read in the Sandia states from OF simulations and Lu19 mechanism
# @author: mhanisnger

# last change: 23.4.19


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pickle
from os.path import join
import Ofpp     # reads OF data to numpy array
import h5py
from utils.clean_states import clean_states_above, clean_states_below
from utils.compute_fBilger import compute_fBilger
import os
import cantera as ct

out_filename = 'Flamelet_states.h5'

gas = ct.Solution('utils/lu19.cti')

species_lu19 = gas.species_names

thermo_fields = ['T','p','f_Bilger']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# create DF for all the states in all Flames and time steps
Data_all_df = pd.DataFrame(columns=all_fields)

store_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/Diffusion_flamelets_database'

path_fine= '/home/max/HDD2_Data/OF4_Simulations/DiffusionFlames/diffusion_lu19_dataGen_fine'
path_normal= '/home/max/HDD2_Data/OF4_Simulations/DiffusionFlames/diffusion_lu19_dataGen'

fine_dirs = os.listdir(path_fine)
normal_dirs = os.listdir(path_normal)

#filter the directory names
fine_dirs = [d for d in fine_dirs if d.startswith('diff_lu19')]
normal_dirs = [d for d in normal_dirs if d.startswith('diff_lu19')]

fine_cases = [join(path_fine,c) for c in fine_dirs]
normal_cases = [join(path_normal,c) for c in normal_dirs]

All_cases = fine_cases+normal_cases
# print(All_cases)

for case in All_cases:
    thisFlame_path = case

    # get the time steps
    dirs = os.listdir(thisFlame_path)

    timesteps = [d for d in dirs if (d.startswith('0.') or d.endswith('e-06') or d.endswith('e-05'))]

    this_df = pd.DataFrame(columns=all_fields)

    # loop over the time step in all Flames
    for time in timesteps:
        thisTime_path = join(thisFlame_path, time)

        for column in this_df:
            # print(os.getcwd())
            print('Reading in %s from time %s of %s' % (column, time, case))

            current_path = join(thisTime_path, column)

            this_time_dir = os.listdir(thisTime_path)
            # # print(this_time_dir)
            # if 'PV' in this_time_dir:
            #     PV_FLAG = True
            # else:
            PV_FLAG = False

            # if 'f_Bilger' in this_time_dir:
            f_Bilger_FLAG = True
            # else:
            #     f_Bilger_FLAG = False

            if os.path.exists(current_path):
                this_df[column] = Ofpp.parse_internal_field(current_path)  # column is also the field name here
            else:
                #raise FileNotFoundError
                #print('%s Not Found!' % column)
                this_df[column] = 0

        # compute fBilger
        Y = this_df[species_lu19].values
        #print('Y.shape', Y.shape)
        f_Bilger = compute_fBilger(Y)
        #print('f_Bilger: ',f_Bilger)
        this_df['f_Bilger'] = f_Bilger


        # print(' ')
        # print(this_df.head())
        # print('  ')

        # if PV_FLAG:
        #     PV_max = max(this_df['PV'])

        # CLEAN THE DF!
        # remove all values above f=0.2 --> there is no reaction!
        if f_Bilger_FLAG:
            this_clean_df = clean_states_above(df=this_df, species='f_Bilger', threshold=0.20, sample_size=1)
            this_clean_df = clean_states_below(df=this_df, species='f_Bilger', threshold=0.005, sample_size=1)
        else:
            this_clean_df = this_df

        # if PV_FLAG:
        #     # Clean PV: 0 there is no reaction, PV_max: reaction is finished
        #     this_clean_df = clean_states_above(df=this_df, species='PV', threshold=0.98 * PV_max, sample_size=10000)
        #     this_clean_df = clean_states_below(df=this_df, species='PV', threshold=1e-8, sample_size=10000)

        print('\nData set contains %i entries\n' % len(this_clean_df))
        # append the data to Data_all_df
        Data_all_df = Data_all_df.append(this_clean_df, ignore_index=True)


hdf_database = pd.HDFStore(join(store_path,out_filename))

#update the hdf5 database
hdf_database.append('Flamelet_raw_data',Data_all_df)
hdf_database.close()









