#!/usr/bin/python

'''
This is to integrate chemical states based on given Y,T,P data

@author: mhansinger

last modified: 15.2.2019

PYTHON3.5 !!!
'''

import pandas as pd
import numpy as np
#import dask.dataframe as dd
import cantera as ct
from os.path import join
import os
import scipy.integrate

from utils.ReactorOde import ReactorOde

class Cantera_Reactor_Sandia(object):
    def __init__(self):
        print("Setting up the reactor")
        try:
            self.gas = ct.Solution('utils/lu19.xml')
            print('Lu19 mechanism is used')
        except:
            print('provide .cti file for lu19!')

        self.species_names = self.gas.species_names
        self.species_dict = { row : self.gas.species_index(row) for row in self.gas.species_names}
        self.T_ref = 298.15     # reference temperature for enthalpy
        #self.reactor = ct.IdealGasReactor(self.gas)
        #self.sim = ct.ReactorNet([self.reactor])

        self.data_integrated = None

        #dt = input('Chose a dt in ms for the integration:')
        self.dt = 1e-6
        self.P = 101325


    def read_data(self,path='/home/max/HDD2_Data/OF4_Simulations/Sandia/Sandia_Database/'):

        self.Sandia_data_path = join(path,'Sandia_states.h5')
        self.Sandia_database_org=pd.read_hdf(self.Sandia_data_path)

        print('These are the data features:')
        print(self.Sandia_database_org.columns)

    def set_tables(self):
        print('reading in tables ...')
        self.read_data()

        print('setting up the output tables ...')

        self.columns_out = [n+'_after' for n in self.species_names]
        self.columns_out.append('T_after')

        # difference to initial value
        columns_out2 = ['d'+n.split('_')[0] for n in self.columns_out]

        columns_out3 = ['RR_' + n.split('_')[0] for n in self.species_names]
        self.columns_out = self.columns_out + columns_out2 + columns_out3

        self.columns_out.append('dt')
        self.columns_out.append('heatRelease')

        temp_np = np.zeros((self.Sandia_database_org.shape[0],len(self.columns_out)))

        temp_pd = pd.DataFrame(temp_np,columns=self.columns_out)

        self.data_integrated = pd.concat([self.Sandia_database_org,temp_pd],axis=1)


    def integrate_Ode(self, iloc):
        # CVODE integration with cantera states

        initial_state = self.Sandia_database_org.iloc[iloc]

        # setting up the gas and solver
        self.gas.TPY = initial_state['T'], self.P, initial_state[self.species_names]

        y0 = np.hstack((self.gas.T,self.gas[self.gas.species_names].Y))

        ode = ReactorOde(self.gas)

        solver = scipy.integrate.ode(ode)
        solver.set_integrator('vode', method='bdf', with_jacobian=True)
        solver.set_initial_value(y0, 0.0)

        # do the forward integration
        solver.integrate(solver.t + self.dt)

        y1 = np.hstack([ self.gas.T,self.gas[self.gas.species_names].Y])

        print(y1)
        print(y0-y1)

        # Überlegen wie der output am besten wäre...


    # def integrate_cantera(self,iloc):
    #
    #     this_df = self.Sandia_database_org.iloc[iloc]
    #
    #     this_Y = []
    #     for n in self.species_names:
    #         this_Y.append(this_df[n])
    #
    #     this_TPY = this_df['T'], this_df['p'], this_Y
    #
    #     self.gas.TPY = this_TPY
    #
    #     # create reactor
    #     reactor = ct.IdealGasConstPressureReactor(self.gas)
    #     sim = ct.ReactorNet([reactor])
    #
    #     # advance in time
    #     sim.advance(self.dt)
    #
    #     # return reactor object
    #     return reactor





if __name__ == '__main__':
    myReact = Cantera_Reactor_Sandia()
    myReact.set_tables()

    myReact.integrate_Ode(1000)

    #reactor = myReact.integrate_cantera(iloc=0)


