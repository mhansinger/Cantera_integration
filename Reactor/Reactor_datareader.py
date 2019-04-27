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

out_filename = 'Reactor_states_2.h5'

gas = ct.Solution('utils/lu19.cti')

species_lu19 = gas.species_names

thermo_fields = ['T','p','f_Bilger']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# create DF for all the states in all Flames and time steps
Data_all_df = pd.DataFrame(columns=all_fields)

store_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/Reactor_database'

path= '/home/max/HDD2_Data/OF4_Simulations/CanteraIgnition/DataGeneration/IgniteLu19/Cases'

dirs = os.listdir(path)

#filter the directory names
dirs = [d for d in dirs if d.startswith('Lu19_f')]

cases = [join(path, c) for c in dirs]

for case in cases:
    thisFlame_path = case

    # get the time steps
    dirs = os.listdir(thisFlame_path)

    #print(dirs)

    timesteps = [d for d in dirs if (d.startswith('0.') or d.endswith('e-06') or d.endswith('e-05'))]

    # loop over the time step in all Flames
    for time in timesteps:
        thisTime_path = join(thisFlame_path, time)

        this_df = pd.DataFrame(data=np.zeros((1,len(all_fields))),columns=all_fields)

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
            f_Bilger_FLAG = False
            # else:
            #     f_Bilger_FLAG = False

            if os.path.exists(current_path):
                this_df[column] = Ofpp.parse_internal_field(current_path)  # column is also the field name here
                #print(this_df)
            else:
                this_df[column] = 0

        # compute fBilger
        Y = this_df[species_lu19].values
        #print('Y.shape', Y.shape)
        f_Bilger = compute_fBilger(Y)
        #print('f_Bilger: ',f_Bilger)
        this_df['f_Bilger'] = f_Bilger

        # append the data to Data_all_df
        Data_all_df = Data_all_df.append(this_df, ignore_index=True)


hdf_database = pd.HDFStore(join(store_path,out_filename))

#update the hdf5 database
hdf_database.append('Reactor_raw_data',Data_all_df)
hdf_database.close()









