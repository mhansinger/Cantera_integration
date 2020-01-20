'''
This is to shuffle the h5 data base

last change: AUGUST 2019

author: mhansinger
'''

import pandas as pd
import numpy as np
import dask.dataframe as dd
import cantera as ct
from os.path import join
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

from utils.compute_fBilger import compute_fBilger

from random import sample
from dask import compute, delayed

from sklearn.decomposition import PCA


class TNF_shuffle_filter(object):
    def __init__(self,augment_DB):

        self.columns = None
        self.data_integrated_dd = None
        self.nr_rows = None
        self.name= None
        self.augment_DB = augment_DB

        if self.augment_DB is True:
            self.gas = ct.Solution('utils/lu19.cti')
            print('Lu19 mechanism is used')

            self.dt = 1e-7
            self.p = 101325

            self.species_names = self.gas.species_names

    def read_data_dd(self, name, key='TNF_data_integrated', path='/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database'):

        self.name = name
        self.path = path
        # point to the hdf library
        self.data_integrated_dd = dd.read_hdf(join(path,name),key=key)

        self.data_integrated_dd = self.data_integrated_dd.reset_index(drop=True)
        self.columns = self.data_integrated_dd.columns
        self.nr_rows = self.data_integrated_dd.shape[0].compute()

    def shuffle_data(self):
        # randomize the order of the data set
        self.data_integrated_dd = self.data_integrated_dd.sample(frac=1).reset_index(drop=True)
        print('Database is shuffled!')

    def filter_data(self,main_RR,T_min,f_max,quantile):
        '''

        :param main_RR: species reaction rate
        :param threshold: below which value the data is thinned out
        :param T_min: remove values below T_min
        :return:
        '''

        self.shuffle_data()

        self.data_integrated_dd = self.data_integrated_dd[abs(self.data_integrated_dd['T']) > T_min]
        print('Removed all %s < %f' % ('T_min',T_min))
        #print('Median of %s: ', self.data_integrated_dd[main_RR].median().compute() % main_RR )

        self.data_integrated_dd = self.data_integrated_dd[(self.data_integrated_dd['f_Bilger']) > 1e-20]
        print('Removed all f_Bilger == 0')

        self.data_integrated_dd = self.data_integrated_dd[(self.data_integrated_dd['f_Bilger']) < f_max]
        print('Removed all f_Bilger < %f' % f_max)

        all_indexes = self.data_integrated_dd.index.compute()
        all_indexes = all_indexes.to_list()
        print('index to_list done...')

        # compute quantile before
        quantile_intially = self.data_integrated_dd[main_RR].quantile(quantile).compute()

        # compute quantile of main_RR
        print('Quantile %f of %s: %f' % (quantile, main_RR, quantile_intially))

        # create a dask.dataframe where the main_RR is fulfilled
        self.data_integrated_dd_filtered = self.data_integrated_dd[abs(self.data_integrated_dd[main_RR]) >= abs(quantile_intially)]
        removed_data = self.data_integrated_dd[abs(self.data_integrated_dd[main_RR]) < abs(quantile_intially)].reset_index(drop=True)

        # store the indexes of the removed values, a share of them needs to be appended
        remove_list = self.data_integrated_dd[abs(self.data_integrated_dd[main_RR]) < abs(quantile_intially)].index.compute()
        remove_list = remove_list.to_list()
        print('remove list done ...')
        
        # compute the ratio of removal to all indexes
        ratio = len(remove_list) / len(all_indexes)
        print('ratio values to remove: ', ratio)

        # # copmute quantile after removal
        self.quantile_after = self.data_integrated_dd_filtered[main_RR].quantile(quantile).compute()
        print('Quantile %f of %s after removal: %f' % (quantile, main_RR, self.quantile_after))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # augmenting the database!
        if self.augment_DB is True:
            self.augment_database()

        # clear the data base from unphysical values
        self.data_ODE_augmented_dd = self.data_ODE_augmented_dd[self.data_ODE_augmented_dd[main_RR] < 0].sample(frac=1).reset_index(drop=True)
        self.data_ODE_augmented_dd = self.data_ODE_augmented_dd[abs(self.data_ODE_augmented_dd[main_RR]) >= abs(quantile_intially)]
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        self.removed_data = removed_data.sample(frac=0.001).reset_index(drop=True)

        self.data_filtered_merged_dd = dd.concat([self.data_integrated_dd_filtered,self.data_ODE_augmented_dd]) #  dd.merge(self.data_integrated_dd_filtered, self.removed_data)
        print('Database is filtered!')

        self.data_filtered_merged_dd = dd.concat([self.data_filtered_merged_dd, self.removed_data])

        # copmute quantile after removal
        quantile_merge = self.data_filtered_merged_dd[main_RR].quantile(quantile).compute()
        print('Quantile %f of %s merge: %f' % (quantile, main_RR, quantile_merge))


    def create_subset(self,frac = 0.1):
        print('Creating subset ...')
        self.data_subset = self.data_filtered_merged_dd.sample(frac=frac).reset_index(drop=True).compute()


    def plot_subset(self,x='f_Bilger',y='RR_CH4',color_by='T'):
        # scatter a subset of data
        plt.scatter(self.data_subset[x],self.data_subset[y],c=self.data_subset[color_by],s=0.2)
        plt.set_cmap('hot')

        #plt.title('Subset of data, colored by: ',color_by)
        plt.xlabel(x)
        plt.ylabel(y)
        plt.colorbar()
        plt.show(block=False)


    # CANTERA integration if DB is to augment
    def augment_database(self):
        print('Augmenting database ... ')

        self.data_ODE_np = np.zeros((1, len(self.columns)))

        #loop over the 'high RR_CH4' data
        # ONCE over all high data values
        print('\nIntegrate all high values ... ')
        for idx_set, this_set in tqdm(self.data_integrated_dd_filtered.iterrows()):
            this_RR_CH4 = this_set['RR_CH4']

            this_y = this_set[self.species_names]
            this_T = this_set['T']
            this_f_Bilger = this_set['f_Bilger']

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            T_after, Y_after, RR_after, rho_after = self.ODE_integrate(this_y,this_T,self.p)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            update_vector_np = np.hstack((this_y, this_T, Y_after, T_after, RR_after, this_f_Bilger))

            self.data_ODE_np = np.vstack((self.data_ODE_np,update_vector_np))

        #self.data_ODE_augmented_dd = dd.from_array(x=self.data_ODE_np,columns=self.columns)
        #self.data_ODE_augmented_dd = dd.from_array(x=self.data_ODE_np,columns=self.columns)

        # compute the 0.1 quantile as selection criteria for the integration with random species changes...
        quantile_01 = self.data_integrated_dd['RR_CH4'].quantile(0.05).compute()

        # NOW loop over values where RR_CH4 > quantile_01
        print('\nIntegrate with random where RR_CH4 > quantile_01 ... ')
        for idx_set, this_set in tqdm(self.data_integrated_dd_filtered.iterrows()):

            this_RR_CH4 = this_set['RR_CH4']

            if this_RR_CH4 < quantile_01:

                for _ in range(10):
                    this_set_altered = this_set

                    # #Modify some species and Temperature
                    this_set_altered['CH4'] = this_set['CH4'] * (1 + np.random.randn() * 1e-3)
                    this_set_altered['CO2'] = this_set['CO2'] * (1 + np.random.randn() * 1e-3)
                    # this_set_altered['O2'] = this_set['O2'] * (1 + np.random.randn() * 1e-3)
                    # this_set_altered['H2O'] = this_set['H2O'] * (1 + np.random.randn() * 1e-3)
                    this_set_altered['OH'] = this_set['OH'] * (1 + np.random.randn() * 1e-4)
                    this_set_altered['H'] = this_set['H'] * (1 + np.random.randn() * 1e-4)

                    this_y = this_set_altered[self.species_names]
                    this_T = this_set_altered['T'] * (1 + np.random.randn() * 5e-3)

                    # if this_T > 2200:
                    #     this_T = this_set_altered['T']

                    # this_T = this_set_altered['T']

                    # recompute f_Bilger because species composition has changed!
                    this_f_Bilger = compute_fBilger(this_y.values)

                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    T_after, Y_after, RR_after, rho_after = self.ODE_integrate(this_y, this_T, self.p)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    update_vector_np = np.hstack((this_y, this_T, Y_after, T_after, RR_after, this_f_Bilger))

                    self.data_ODE_np = np.vstack((self.data_ODE_np, update_vector_np))

            self.data_ODE_augmented_dd = dd.from_array(x=self.data_ODE_np, columns=self.columns)



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



    def write_hdf(self,nameDB,dt='1e-7'):
        # not yet implemented!
        # try to use chunks and write new hdf in smaller tables, this is useful for the batch training then

        # # use partitions
        # npart = round(self.nr_rows/100)
        # parted_df = self.data_integrated_dd.repartition(npartitions=npart)

        # hdf_database = pd.HDFStore(join('/home/max/HDD2_Data/OF4_Simulations/ANN_Lu19_data/TNF_database','TNF_integrated_filtered_dt%s.h5' % str(dt)),
        #                            encoding='utf-8')
        #
        # # update the hdf5 database
        # hdf_database.append(nameDB, self.data_integrated)
        # hdf_database.close()

        new_name = join(self.path,'TNF_integrated_filtered_dt%s.h5' % str(dt))
        #self.data_integrated_dd.to_hdf(path_or_buf=new_name,key=nameDB)
        self.data_integrated_dd.to_hdf(path_or_buf=new_name, key=nameDB)
        print('Database is written ... ')



if __name__=='__main__':
    TNF = TNF_shuffle_filter(augment_DB=True)
    TNF.read_data_dd(name='TNF_integrated_sample_dt1e-07.h5')
    TNF.filter_data(main_RR='RR_CH4', T_min = 800,f_max=0.15,quantile=0.2)

    TNF.data_integrated_dd.compute()
    TNF.create_subset(frac=0.9)
    TNF.plot_subset(x='f_Bilger', y='RR_CH4', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_H2', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_CO', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_OH', color_by='T')
    TNF.write_hdf(nameDB='TNF_filtered',dt='1e-7')

    # works ... August, 2019

    plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['RR_CH4'].compute(),s=0.2,c=TNF.data_ODE_augmented_dd['T'].compute())
    plt.xlabel('f_Bilger')
    plt.ylabel('RR_CH4')
    plt.show()

    plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['RR_OH'].compute(),s=0.2,c=TNF.data_ODE_augmented_dd['T'].compute())
    plt.xlabel('f_Bilger')
    plt.ylabel('RR_OH')
    plt.show()

    # plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['T'].compute(),s=0.2)
    # plt.show()

# # %%
#     TNF.plot_subset(x='CH4', y='OH', color_by='RR_CO')
#
# # %%
#     plt.scatter(TNF.data_ODE_augmented_dd['CH4'].compute(),TNF.data_ODE_augmented_dd['OH'].compute(),s=0.2,c=TNF.data_ODE_augmented_dd['T'].compute())
#     plt.show()
#
    plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['T'].compute(),s=0.2)
    plt.show()

