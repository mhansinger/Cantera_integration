# read in the Sandia states from OF simulations and Lu19 mechanism
# @author: mhanisnger

# last change: 11.2.19


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pickle
from os.path import join
import Ofpp     # reads OF data to numpy array
import h5py
from clean_states import clean_states
import os


out_filename = 'Sandia_states.h5'

species_lu19 = ['C2H2','C2H4', 'C2H6', 'CH2CO', 'CH2O', 'CH3', 'CH3OH', 'CH4', 'CO', 'CO2', 'H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH']

thermo_fields = ['T','rho','p','heatRelease','thermo:mu','thermo:alpha']

all_fields = species_lu19 + thermo_fields

# create DF for all the states in all Flames and time steps
Data_all_df = pd.DataFrame(columns=all_fields)

Sandia_path = '/home/max/HDD2_Data/OF4_Simulations/Sandia'

Flames = ['FlameD_lu19_lam', 'Flame-E_lu19_lam', 'Flame-F_lu19_lam']

for Flame in Flames:
    thisFlame_path = join(Sandia_path,Flame)

    # get the time steps
    dirs = os.listdir(thisFlame_path)
    timesteps = [d for d in dirs if d.startswith('0.')]

    this_df = pd.DataFrame(columns=all_fields)

    # loop over the time step in all Flames
    for time in timesteps:
        try:
            thisTime_path = join(thisFlame_path, time)

            for column in this_df:
                print(os.getcwd())
                print('Reading in %s from time %s of %s' % (column, time, Flame))

                current_path = join(thisTime_path,column)

                if os.path.exists(current_path):
                    this_df[column] = Ofpp.parse_internal_field(current_path)  # column is also the field name here
                else:
                    raise FileNotFoundError

            print(' ')
            print(this_df.head())
            print('  ')

            # CLEAN THE DF!
            this_df = clean_states(df=this_df,species='T',threshold=295,sample_size=1e4)

            # append the data to Data_all_df
            Data_all_df = Data_all_df.append(this_df, ignore_index=True)

        except:
            pass


outPath = join(Sandia_path,'Sandia_Database')

hdf_database = pd.HDFStore(join(outPath,out_filename))

#update the hdf5 database
hdf_database.append('Sandia_Data',Data_all_df)
hdf_database.close()









