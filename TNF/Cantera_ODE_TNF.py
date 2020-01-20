#!/usr/bin/python

'''
This is to integrate chemical states based on given Y,T,P data

@author: mhansinger

last modified: 16.4.2019

PYTHON3.5 !!!
'''

import pandas as pd
import numpy as np
import dask.dataframe as dd
import cantera as ct
from os.path import join
import os
from tqdm import tqdm
#import matplotlib.pyplot as plt

from random import sample
from dask import compute, delayed

# import scipy.integrate
#
# from utils.ReactorOde import ReactorOde


#ignore cantera warnings
ct.suppress_thermo_warnings()

class Cantera_ODE_TNF(object):
    def __init__(self):
        print("Setting up the reactor")
        try:
            self.gas = ct.Solution('utils/lu19.cti')
            print('Lu19 mechanism is used')
        except:
            print('provide .cti file for lu19!')

        self.species_names = self.gas.species_names
        self.species_dict = {row : self.gas.species_index(row) for row in self.gas.species_names}
        self.T_ref = 298.15     # reference temperature for enthalpy
        #self.reactor = ct.IdealGasReactor(self.gas)
        #self.sim = ct.ReactorNet([self.reactor])

        self.data_integrated = None

        #dt = input('Chose a dt in ms for the integration:')
        self.dt = 1e-7
        self.P = 101325
        self.len_dataset = None
        self.TNF_database_org = None
        print('Time step is: %s' % str(self.dt))

    def read_data(self,name):

        self.TNF_data_path = join(self.path,name)
        self.TNF_database_org=pd.read_hdf(self.TNF_data_path)

        print('These are the data features:')
        print(self.TNF_database_org.columns)

    def read_data_dd(self, name):

        self.TNF_data_path = join(self.path, name)
        self.TNF_database_org=dd.read_hdf(self.TNF_data_path, key=self.key)

        print('These are the data features:')
        print(self.TNF_database_org.columns)

        # plot some histograms
        #self.plot_histograms('T')
        #self.plot_histograms('CO')
        #self.plot_histograms('f_Bilger')
        #self.plot_histograms('OH')
        #self.plot_histograms('CO2')
        #self.plot_histograms('CH4')

    def set_tables(self,name,path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database',key='TNF_raw_data'):
        print('reading in tables ...')
        
        self.path = path
        self.key = key
        #self.read_data(name=name)
        self.read_data_dd(name=name)

        print('setting up the output tables ...')

        self.columns_out = self.species_names.copy()
        self.columns_out.append('T')

        self.columns_out_after = [n+'_after' for n in self.species_names]
        self.columns_out_after.append('T_after')
        # reaction rate [1/s]
        self.columns_RR = ['RR_'+n for n in self.species_names]
        self.columns_out = self.columns_out + self.columns_out_after +self.columns_RR

        self.columns_out.append('f_Bilger')
        self.columns_out.append('rho')

        self.len_dataset = len(self.TNF_database_org)

        self.data_integrated_np = np.zeros((self.len_dataset,len(self.columns_out)))

        print('Output data set columns:')
        print(self.columns_out)


    def ODE_integrate(self,Y,T,p,steps):
        # do the ODE integration
        # Y is the species vector, T temperature, p pressure
        #self.gas = ct.Solution('utils/lu19.cti')
        self.gas.TPY = T, p, Y
        current_rho = self.gas.density
        r = ct.IdealGasConstPressureReactor(self.gas)
        sim = ct.ReactorNet([r])
        time = 0.0

        # do 10 small integration steps to get to dt! --> makes it less stiff
        for n in range(steps):
            time += self.dt/steps
            sim.advance(time)

        T_after = r.thermo.T
        Y_after = r.thermo.Y
        # mass fraction reaction rate dY/dt [1/s] same as in Cantera OF:
        RR_Y = (Y_after - Y) * current_rho / self.dt

        return T_after, Y_after, RR_Y, current_rho


    def loop_ODE(self,remove_T_below,steps):
        print(' ')
        # for row in tqdm(range(self.len_dataset)):
        for idx_fullset, this_set in tqdm(self.TNF_database_org.iterrows()):
            # print('Row number: ',idx_fullset)

            # this_set = self.TNF_database_org.iloc[row].calculate()      # calculate only if NOT dask!
            Y_vector = np.zeros(len(self.gas.species_names))

            # make sure species vector is in the correct order!
            for idx, sp in enumerate(self.gas.species_names):
                #print(idx,sp)
                Y_vector[idx] = this_set[sp]

            this_T = this_set['T']
            this_f_Bilger = this_set['f_Bilger']

            # CRITERIA TO REMOVE UNNECESSARY ENTRIES WHERE T IS AT T=300
            if this_T > remove_T_below:
                try:
                    ###############################
                    # ODE integration
                    T_after, Y_after, RR_Y, current_rho = self.ODE_integrate(Y=Y_vector, T=this_T, p=ct.one_atm,steps=steps)
                    ###############################

                    this_out_vector = np.hstack((Y_vector, this_T, Y_after, T_after, RR_Y, this_f_Bilger, current_rho))
                    self.data_integrated_np[idx_fullset, :] = this_out_vector

                except:
                    print('Error in integration')

        self.data_integrated=pd.DataFrame(data=self.data_integrated_np,columns=self.columns_out)

        # write database
        # self.write_hdf()

    # def filter_shuffle_data(self,condition='RR_CH4',threshold = 2):
    #     #remove entries where there is actually no reaction
    #
    #     # threshold is very sensitive... plot the reaction rates over Z and decide which thershold seems legitimit
    #
    #     all_indexes = self.data_integrated.index.tolist()
    #     print('index to_list done...')
    #
    #     remove_list = self.data_integrated.index[abs(self.data_integrated[condition]) > threshold].tolist()
    #     print('remove list done...')
    #
    #     ratio = len(remove_list) / len(all_indexes)
    #     print('ratio values to remove: ', ratio)
    #     remove_sample = sample(remove_list,int(len(remove_list)/10))
    #     print('sampling done ...')
    #
    #     remove_final = [f for f in remove_list if not f in remove_sample]
    #     print('remove final done ...')
    #
    #     indexes_to_keep = [f for f in all_indexes if not f in remove_final]
    #
    #     self.data_integrated = self.data_integrated.loc[indexes_to_keep]
    #     print('indexes to keep done ...')
    #
    #     #shuffle the data before writing to HDD
    #     #self.data_integrated = self.data_integrated.sample(frac=1).reset_index(drop=True)
    #     self.data_integrated_dd = dd.from_pandas(self.data_integrated,npartitions=5)
    #     self.data_integrated=self.data_integrated_dd.sample(frac=1).reset_index(drop=True).compute()
    #     print('shuffel is done ...')

    def write_hdf(self,nameDB,key='TNF_integrated'):

        # remove zero T values
        self.data_integrated = self.data_integrated[self.data_integrated['T']>0]

        # remove all the zero rows
        self.data_integrated=self.data_integrated.loc[(self.data_integrated!=0).all(1)]

        # reindex the dataset
        self.data_integrated = self.data_integrated.reset_index(drop=True)

        # hdf_database = pd.HDFStore(join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%s.h5' % str(self.dt)),
        #                            encoding='utf-8')
        #
        # # update the hdf5 database
        # hdf_database.append(nameDB, self.data_integrated)
        # hdf_database.close()

        #self.data_integrated.to_hdf(path_or_buf=join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%s' % str(self.dt)))

        path_name = join(self.path,nameDB+'_dt%s.h5' % str(self.dt))
        self.data_integrated.to_hdf(path_or_buf=path_name,key=key,format='table')

'''
    def plot_histograms(self,species):
        this_set = self.TNF_database_org[species].compute()

        if species=='T':
            this_set = this_set[this_set>310]
            plt.figure()
            plt.hist(this_set,bins=100)
            plt.title(species)
            plt.show(block=False)
        else:
            plt.figure()
            plt.hist(this_set,bins=100)
            plt.title(species)
            plt.show(block=False)
'''

# if __name__ == '__main__':
#     myReact = Cantera_ODE_TNF()
#     myReact.set_tables(name='TNF_states.h5')
#     myReact.loop_ODE(remove_T_below=310,steps=1)
#     #myReact.filter_shuffle_data(condition='RR_CH4',threshold=3)
#     myReact.write_hdf(nameDB='TNF_data_integrated')
#     # myReact.integrate_Ode(1000)
#
#     #reactor = myReact.integrate_cantera(iloc=0)


