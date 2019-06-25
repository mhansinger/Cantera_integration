'''
This is to shuffle the h5 data base

last change: June 2019

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

from random import sample
from dask import compute, delayed


class TNF_shuffle_filter(object):
    def __index__(self):

        self.columns = None
        self.data_integrated_dd = None
        self.nr_rows = None
        self.name= None

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

    def create_subset(self,frac = 0.001):
        print('Creating subset ...')
        self.data_subset = self.data_integrated_dd.sample(frac=frac).reset_index(drop=True).compute()


    def plot_subset(self,x='f_Bilger',y='RR_CH4',color_by='T'):
        # scatter a subset of data
        plt.scatter(self.data_subset[x],self.data_subset[y],c=self.data_subset[color_by],s=0.2)
        plt.set_cmap('hot')

        #plt.title('Subset of data, colored by: ',color_by)
        plt.xlabel(x)
        plt.ylabel(y)
        plt.colorbar()
        plt.show(block=False)


    def filter_data(self,condition,threshold):

        self.shuffle_data()

        self.data_integrated_dd = self.data_integrated_dd[abs(self.data_integrated_dd['T']) > 300]
        print('Removed all T < 300')

        self.data_integrated_dd = self.data_integrated_dd[(self.data_integrated_dd['f_Bilger']) > 1e-20]
        print('Removed all f_Bilger == 0')

        all_indexes = self.data_integrated_dd.index.compute()
        all_indexes = all_indexes.to_list()
        print('index to_list done...')


        # create a dask.dataframe where the condition is fulfilled
        self.data_integrated_dd_filtered = self.data_integrated_dd[abs(self.data_integrated_dd[condition]) >= threshold]
        removed_data = self.data_integrated_dd[abs(self.data_integrated_dd[condition]) < threshold]

        # store the indexes of the removed values, a share of them needs to be appended
        remove_list = self.data_integrated_dd[abs(self.data_integrated_dd[condition]) < threshold].index.compute()
        remove_list = remove_list.to_list()
        print('remove list done...')

        # compute the ratio of removal to all indexes
        ratio = len(remove_list) / len(all_indexes)
        print('ratio values to remove: ', ratio)

        self.removed_data = removed_data.sample(frac=0.01).reset_index(drop=True)

        self.data_filtered_merged_dd = dd.concat([self.data_integrated_dd_filtered,self.removed_data]) #  dd.merge(self.data_integrated_dd_filtered, self.removed_data)
        print('Database is filtered!')



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
    TNF = TNF_shuffle_filter()
    TNF.read_data_dd(name='TNF_integrated_sample_dt1e-07.h5')
    TNF.filter_data(condition='RR_CH4',threshold=5)

    TNF.data_integrated_dd.compute()
    TNF.create_subset(frac=0.5)
    TNF.plot_subset(x='f_Bilger',y='RR_CH4',color_by='T')
    TNF.write_hdf(nameDB='TNF_filtered',dt='1e-7')

    # works ... June, 2019
