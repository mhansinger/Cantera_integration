#!/usr/bin/python

'''
This is to integrate chemical states based on given Y,T,P data

@author: mhansinger

last modified: 10.4.2019

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
        self.dt = 1e-6
        self.P = 101325
        self.len_dataset = None
        self.TNF_database_org = None
        print('Time step is: %f' % self.dt)

    def read_data(self,path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'):

        self.TNF_data_path = join(path,'TNF_states.h5')
        self.TNF_database_org=pd.read_hdf(self.TNF_data_path)

        print('These are the data features:')
        print(self.TNF_database_org.columns)

    def read_data_dd(self,path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'):

        self.TNF_data_path = join(path,'TNF_states.h5',)
        self.TNF_database_org=dd.read_hdf(self.TNF_data_path,key='TNF_raw_data')

        print('These are the data features:')
        print(self.TNF_database_org.columns)

        # plot some histograms
        #self.plot_histograms('T')
        #self.plot_histograms('CO')
        #self.plot_histograms('f_Bilger')
        #self.plot_histograms('OH')
        #self.plot_histograms('CO2')
        #self.plot_histograms('CH4')

    def set_tables(self):
        print('reading in tables ...')
        self.read_data_dd()

        print('setting up the output tables ...')

        self.columns_out = self.species_names.copy()
        self.columns_out.append('T')

        self.columns_out_after = [n+'_after' for n in self.species_names]
        self.columns_out_after.append('T_after')
        # reaction rate [1/s]

        self.columns_R = ['R_'+n for n in self.species_names]
        self.columns_out = self.columns_out + self.columns_out_after +self.columns_R

        self.columns_out.append('f_Bilger')
        #self.columns_out.append('dt')

        self.len_dataset = len(self.TNF_database_org)

        self.data_integrated_np = np.zeros((self.len_dataset,len(self.columns_out)))

        #self.data_integrated = pd.DataFrame(columns=self.columns_out)
        #self.data_integrated = pd.concat([self.TNF_database_org,temp_pd],axis=1)
        print('Output data set columns:')
        print(self.columns_out)


    def ODE_integrate(self,Y,T,p):
        # do the ODE integration
        # Y is the species vector, T temperature, p pressure
        #self.gas = ct.Solution('utils/lu19.cti')
        self.gas.TPY = T, p, Y
        r = ct.IdealGasConstPressureReactor(self.gas)
        sim = ct.ReactorNet([r])
        time = 0.0

        # do 10 small integration steps to get to dt! --> makes it less stiff
        for n in range(10):
            time += self.dt/10
            sim.advance(time)

        T_after = r.thermo.T
        Y_after = r.thermo.Y
        R_Y = (Y_after - Y) * self.dt # mass fraction reaction rate dY/dt [1/s]

        return T_after, Y_after, R_Y


    def loop_ODE(self,remove_T_below):
        print(' ')
        #for row in tqdm(range(self.len_dataset)):
        for idx_fullset, this_set in tqdm(self.TNF_database_org.iterrows()):
            #print('Row number: ',idx_fullset)

            #this_set = self.TNF_database_org.iloc[row].calculate()      # calculate only if NOT dask!
            Y_vector = np.zeros(len(self.gas.species_names))

            # make sure species vector is in the correct order!
            for idx, sp in enumerate(self.gas.species_names):
                #print(idx,sp)
                Y_vector[idx] = this_set[sp]

            this_T = this_set['T']
            this_f_Bilger=this_set['f_Bilger']

            # CRITERIA TO REMOVE UNNECESSARY ENTRIES WHERE T IS AT T=300
            if this_T > remove_T_below:
                ###############################
                # ODE integration
                T_after, Y_after, R_Y = self.ODE_integrate(Y=Y_vector, T=this_T, p=ct.one_atm)
                ###############################

                this_out_vector = np.hstack((Y_vector, this_T, Y_after, T_after, R_Y, this_f_Bilger))

                self.data_integrated_np[idx_fullset, :] = this_out_vector

        self.data_integrated=pd.DataFrame(data=self.data_integrated_np,columns=self.columns_out)

        # write database
        self.write_hdf()

    def write_hdf(self):
        hdf_database = pd.HDFStore(join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%f' % self.dt))

        # update the hdf5 database
        hdf_database.append('TNF_data_integrated', self.data_integrated)
        hdf_database.close()

        # self.data_integrated.to_hdf(path_or_buf=join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%f' % self.dt))
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

if __name__ == '__main__':
    myReact = Cantera_ODE_TNF()
    myReact.set_tables()
    myReact.loop_ODE(remove_T_below=310)
    myReact.write_hdf()
    # myReact.integrate_Ode(1000)

    #reactor = myReact.integrate_cantera(iloc=0)

