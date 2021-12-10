import os,sys,re
import pandas as pd
import numpy as np

sa_data = pd.read_csv(sys.argv[1], sep='\t',header=None)
sa_data.columns = ['idx','sample','density','n_acc','n_tot']
sa_data['catalytic'] = sa_data['sample'].str.split('_',expand=True)[1]
#print(sa_data)
for state in np.unique(sa_data['catalytic'].values):
    sub_data = sa_data[sa_data['catalytic'] == state]
    for density in np.unique(sub_data['density'].values):
        for i in np.unique(sub_data['idx'].values):
            sub_data_des = sub_data[(sub_data['density'] == density) & (sub_data['idx'] == i)]
            n_acc = np.sum(sub_data_des['n_acc'])
            n_tot = np.sum(sub_data_des['n_tot'])
            print("%s\t%s\t%s\t%s\t%s" % (i, state, density, n_acc, n_tot))


