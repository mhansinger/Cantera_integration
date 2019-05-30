# read in the TNF states from OF simulations and Lu19 mechanism
# @author: mhanisnger

# last change: 10.5.19


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pickle
from os.path import join
import Ofpp     # reads OF data to numpy array
import h5py
from utils.clean_states import clean_states_above, clean_states_below
import os
from dask.delayed import delayed
import dask.dataframe as dd

out_filename = 'TNF_states_DNS.h5'

species_lu19 = ['C2H2','C2H4', 'C2H6', 'CH2CO', 'CH2O', 'CH3', 'CH3OH', 'CH4', 'CO', 'CO2', 'H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH']

thermo_fields = ['T','f_Bilger']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# create DF for all the states in all Flames and time steps
Data_all_df = pd.DataFrame(columns=all_fields)

store_path = '/home/hansinger/DNS_database'

TNF_path_DNS= ['/home/hansinger/DNS_database']

for case in TNF_path_DNS:
    thisFlame_path = case #join(TNF_path,case)

    # get the time steps
    dirs = os.listdir(thisFlame_path)

    #print(dirs)

    timesteps = [d for d in dirs if d.startswith('0_')]

    this_df = pd.DataFrame(columns=all_fields)

    # loop over the time step in all Flames
    for time in timesteps:
        thisTime_path = join(thisFlame_path, time)

        for column in this_df:
            # print(os.getcwd())
            print('Reading in %s from time %s of %s' % (column, time, case))

            current_path = join(thisTime_path, column+'.gz')

            this_time_dir = os.listdir(thisTime_path)
            # print(this_time_dir)
            if 'PV.gz' in this_time_dir:
                PV_FLAG = True
            else:
                PV_FLAG = False

            if 'f_Bilger.gz' in this_time_dir:
                f_Bilger_FLAG = True
            else:
                f_Bilger_FLAG = False

            if os.path.exists(current_path):
                data = pd.read_csv(current_path, compression='gzip',names=[column])[column] #Ofpp.parse_internal_field(current_path)  # column is also the field name here
                this_df[column] = data.values
            else:
                #raise FileNotFoundError
                print('%s Not Found!' % column)
                this_df[column] = 0

        print(' ')
        print(this_df.head())
        print('  ')

        if PV_FLAG:
            PV_max = max(this_df['PV'])

        # CLEAN THE DF!
        # remove all values above f=0.2 --> there is no reaction!
        if f_Bilger_FLAG:
            this_clean_df = clean_states_above(df=this_df, species='f_Bilger', threshold=0.15, sample_size=10000)
            this_clean_df = clean_states_below(df=this_df, species='f_Bilger', threshold=0.005, sample_size=10000)

        if PV_FLAG:
            # Clean PV: 0 there is no reaction, PV_max: reaction is finished
            this_clean_df = clean_states_above(df=this_df, species='PV', threshold=0.98 * PV_max, sample_size=10000)
            this_clean_df = clean_states_below(df=this_df, species='PV', threshold=1e-8, sample_size=10000)

        print('\nData set contains %i entries\n' % len(this_clean_df))
        # append the data to Data_all_df
        Data_all_df = Data_all_df.append(this_clean_df, ignore_index=True)


hdf_database = pd.HDFStore(join(store_path,out_filename))

#update the hdf5 database
hdf_database.append('TNF_raw_DNS',Data_all_df)
hdf_database.close()









