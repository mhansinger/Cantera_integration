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

    def read_data(self,name,path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'):

        self.TNF_data_path = join(path,name)
        self.TNF_database_org=pd.read_hdf(self.TNF_data_path)

        print('These are the data features:')
        print(self.TNF_database_org.columns)

    def read_data_dd(self, name,path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'):

        self.TNF_data_path = join(path, name)
        self.TNF_database_org=dd.read_hdf(self.TNF_data_path, key='TNF_raw_data')


        print('These are the data features:')
        print(self.TNF_database_org.columns)

        # plot some histograms
        #self.plot_histograms('T')
        #self.plot_histograms('CO')
        #self.plot_histograms('f_Bilger')
        #self.plot_histograms('OH')
        #self.plot_histograms('CO2')
        #self.plot_histograms('CH4')

    def set_tables(self,name):
        print('reading in tables ...')
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
        RR_Y = (Y_after - Y) * current_rho / self.dt # mass fraction reaction rate dY/dt [1/s] same as in Cantera OF

        return T_after, Y_after, RR_Y


    def loop_ODE(self,remove_T_below,steps):
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
                try:
                    ###############################
                    # ODE integration
                    T_after, Y_after, RR_Y = self.ODE_integrate(Y=Y_vector, T=this_T, p=ct.one_atm,steps=steps)
                    ###############################

                    this_out_vector = np.hstack((Y_vector, this_T, Y_after, T_after, RR_Y, this_f_Bilger))
                    self.data_integrated_np[idx_fullset, :] = this_out_vector

                except:
                    print('Error in integration')

        self.data_integrated=pd.DataFrame(data=self.data_integrated_np,columns=self.columns_out)

        # write database
        self.write_hdf()

    def write_hdf(self,nameDB):
        hdf_database = pd.HDFStore(join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%s' % str(self.dt)))

        # update the hdf5 database
        hdf_database.append(nameDB, self.data_integrated)
        hdf_database.close()

        self.data_integrated.to_hdf(path_or_buf=join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_dt%s' % str(self.dt)))
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
    myReact.set_tables(name='TNF_states.h5')
    myReact.loop_ODE(remove_T_below=310,steps=1)
    myReact.write_hdf(nameDB='TNF_data_integrated.h5')
    # myReact.integrate_Ode(1000)

    #reactor = myReact.integrate_cantera(iloc=0)


