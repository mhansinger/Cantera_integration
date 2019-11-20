#!/usr/bin/python

'''
This file creates a data base based on reactor simulations

@author: mhansinger

last modified: July 2019

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
from numba import jit
from utils.compute_fBilger import compute_fBilger

from dask.delayed import delayed
import dask
import time

import matplotlib.pyplot as plt



class creacte_reactor_database(object):
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
        self.t_end = 1e-3

        self.CH4_range = np.linspace(0.005,0.25,1000) #np.linspace(0.03,0.08,1000)
        self.Y_O2_air = 0.233
        self.Y_N2_air = 1- self.Y_O2_air

        # indexes for cantera species
        self.CH4_INDEX = 9
        self.N2_INDEX = 18
        self.O2_INDEX = 3
        self.H2_INDEX = 0
        self.H_INDEX = 1
        self.OH_INDEX = 4
        self.CO_INDEX =10
        self.CO2_INDEX = 11

        self.PATH_DATABASE = '/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/Reactor_database'

        print('Time step is: %s' % str(self.dt))

    def set_tables(self):

        print('setting up the output tables ...')

        self.columns_input = self.species_names.copy()
        self.columns_input.append('T')
        self.columns_input.append('rho')

        self.columns_out = [n+'_after' for n in self.species_names]
        self.columns_out.append('T_after')
        # reaction rate [1/s]
        self.columns_RR = ['RR_'+n for n in self.species_names]
        self.columns_out = self.columns_input + self.columns_out + self.columns_RR

        self.columns_out.append('f_Bilger')
        #self.columns_out.append('rho')

        self.len_dataset = len(self.columns_out)

        self.data_integrated_np = np.zeros((1,len(self.columns_out)))

        print('Output data set columns:')
        print(self.columns_out)

    def get_index(self,name):
        # get the index of the species 'name' in self.columns_out
        return self.columns_out.index(name)

    def run_ignition(self,with_radicals=True):
        # set initial state of reactor to be ignited

        T_init = (1000 + np.random.random(1)[0]*1000)

        Y_CH4_init = self.CH4_range[np.random.randint(len(self.CH4_range))]

        if with_radicals:
            #print('NOT IMPLEMENTED!')
            Y_H_init = 1e-4 *np.random.random(1)[0]
            Y_OH_init = 0 #1e-3  * np.random.random(1)[0]
        else:
            Y_H_init = 0
            Y_OH_init = 0

        Y_O2_init = (1 - Y_CH4_init - Y_H_init - Y_OH_init) * self.Y_O2_air
        Y_N2_init = (1 - Y_CH4_init - Y_H_init - Y_OH_init) * self.Y_N2_air

        this_time = 0

        Y_initial = np.zeros(len(self.species_names))
        Y_initial[self.CH4_INDEX] = Y_CH4_init
        Y_initial[self.N2_INDEX] = Y_N2_init
        Y_initial[self.O2_INDEX] = Y_O2_init
        Y_initial[self.H_INDEX] = Y_H_init

        print(Y_CH4_init)
        print(T_init)
        print(' ')

        f_Bilger = compute_fBilger(Y_initial.T)

        data_integrated_np = np.zeros((1,len(self.columns_out)))

        while this_time<= self.t_end:
            T_after, Y_after, RR_Y, current_rho =self.ODE_integrate(Y=Y_initial, T=T_init, p=self.P)

            # update time
            this_time += self.dt

            update_vector = np.hstack([Y_initial,T_init, current_rho,Y_after,T_after,RR_Y,f_Bilger])

            T_init = T_after
            Y_initial = Y_after

            # update the data base
            self.data_integrated_np = np.vstack([self.data_integrated_np,update_vector])

            ## break criteria
            if abs(RR_Y[self.CH4_INDEX]) < 0.05:
                break


    def loop_reactor_ignition(self,loops):

        for i in range(loops):
            self.run_ignition()


    def ODE_integrate(self,Y,T,p):
        # do the ODE integration
        # Y is the species vector, T temperature, p pressure
        #self.gas = ct.Solution('utils/lu19.cti')
        self.gas.TPY = T, p, Y
        current_rho = self.gas.density
        r = ct.IdealGasConstPressureReactor(self.gas)
        sim = ct.ReactorNet([r])

        time = 0
        # do 10 small integration steps to get to dt! --> makes it less stiff
        #for n in range(steps):
        time += self.dt
        sim.advance(time)

        T_after = r.thermo.T
        Y_after = r.thermo.Y
        RR_Y = (Y_after - Y) * current_rho / self.dt # mass fraction reaction rate dY/dt [1/s]

        return T_after, Y_after, RR_Y, current_rho

    def write_database(self,key,name):
        file_name = name+'_dt%s.h5' % str(self.dt)

        data_base_pd = pd.DataFrame(self.data_integrated_np,columns=self.columns_out)
        data_base_pd.to_hdf(join(self.PATH_DATABASE,file_name),key=key,format='table')


if __name__ == '__main__':
    test = creacte_reactor_database()
    test.set_tables()
    test.loop_reactor_ignition(500)
    test.write_database(name='reactor_integrated',key='reactor_integrated')