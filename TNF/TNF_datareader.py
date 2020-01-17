# read in the TNF states from OF simulations and Lu19 mechanism
# @author: mhanisnger

# last change: 11.3.19

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


out_filename = 'TNF_states.h5'

species_lu19 = ['C2H2','C2H4', 'C2H6', 'CH2CO', 'CH2O', 'CH3', 'CH3OH', 'CH4', 'CO', 'CO2', 'H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH']

thermo_fields = ['T','p','f_Bilger','PV']#,'rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# create DF for all the states in all Flames and time steps
Data_all_df = pd.DataFrame(columns=all_fields)

store_path = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'

TNF_path_SuperMUC= '/home/max/HDD2_Data/OF4_Simulations/TNF_KIT/SuperMUC'
TNF_path_Cluster= '/home/max/HDD2_Data/OF4_Simulations/TNF_KIT/Cluster'


cases_SuperMUC = ['case-05_ESF_linear','case-05', 'case-05_lam_lowMa', 'case-05_lowT_lam','case-05_lam_lowMa_4Ord',
                  'case-44_react_lam_WALE','case-42_react_ESF_WALE_1Mio']
cases_Cluster=['case-27_lam_WALE_2.8Mio_linear','case-43_ESF_WALE_2.8Mio','case-43_ESF_WALE_2.8Mio_linear','case-46_lam_WALE_550k_Pr0.2']

MUC_paths = [join(TNF_path_SuperMUC,f) for f in cases_SuperMUC]
Cluster_paths = [join(TNF_path_Cluster,f) for f in cases_Cluster]

All_cases = MUC_paths + Cluster_paths

for case in All_cases:
    thisFlame_path = case #join(TNF_path,case)

    # get the time steps
    dirs = os.listdir(thisFlame_path)

    #print(dirs)

    timesteps = [d for d in dirs if d.startswith('0.')]

    this_df = pd.DataFrame(columns=all_fields)

    # loop over the time step in all Flames
    for time in timesteps:
        thisTime_path = join(thisFlame_path, time)

        for column in this_df:
            # print(os.getcwd())
            print('Reading in %s from time %s of %s' % (column, time, case))

            current_path = join(thisTime_path, column)

            this_time_dir = os.listdir(thisTime_path)
            # print(this_time_dir)
            if 'PV' in this_time_dir:
                PV_FLAG = True
            else:
                PV_FLAG = False

            if 'f_Bilger' in this_time_dir:
                f_Bilger_FLAG = True
            else:
                f_Bilger_FLAG = False

            if os.path.exists(current_path):
                this_df[column] = Ofpp.parse_internal_field(current_path)  # column is also the field name here
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
            this_clean_df = clean_states_above(df=this_df, species='f_Bilger', threshold=0.4, sample_size=10000)
            this_clean_df = clean_states_below(df=this_df, species='f_Bilger', threshold=0.001, sample_size=10000)

        if PV_FLAG:
            # Clean PV: 0 there is no reaction, PV_max: reaction is finished
            this_clean_df = clean_states_above(df=this_df, species='PV', threshold=0.98 * PV_max, sample_size=10000)
            this_clean_df = clean_states_below(df=this_df, species='PV', threshold=1e-8, sample_size=10000)

        print('\nData set contains %i entries\n' % len(this_clean_df))
        # append the data to Data_all_df
        Data_all_df = Data_all_df.append(this_clean_df, ignore_index=True)


hdf_database = pd.HDFStore(join(store_path,out_filename))

#update the hdf5 database
hdf_database.append('TNF_raw_data',Data_all_df)
hdf_database.close()









