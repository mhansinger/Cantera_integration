
# removes unnecessary thermodynamic-states from the given dataframe
# define trigger species and threshold

import numpy as np
import pandas as pd

def clean_states(df,species='f_Bilger',threshold=0.25,sample_size=2e5):

    print('\nCleaning data set')
    print('Removing all %s which are above %.2f' %(species,threshold))

    # get index values bigger than threshold
    index1 = list(df[species][df[species]<threshold].index)

    # sample data for
    index2 = list(df[species][df[species]>=threshold].sample(int(sample_size)).index)

    index_list = index1 + index2

    new_df = df.iloc[index_list]

    return new_df