import scanpy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy
from pylab import *
from scipy import signal
import pickle
import seaborn as sns

def ccf(x, y):
    result = np.correlate(y - np.mean(y), x - np.mean(x), 'same') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    return result[length:]

def xcor(a,b, maxlen = 1000):
    foo = []
    for i in range(len(a)):
        res = ccf(a[i],b[i])
        foo.append(res)
    return foo

def process_xcors(new_mat, max_len=500):
    auto_cor = xcor(new_mat, new_mat, maxlen=max_len)
    return auto_cor

mat_total = np.load('length_filtered_molecules_tfs.npy')
mat_total = np.nan_to_num(pd.DataFrame(mat_total).rolling(33,axis=1, center=True, min_periods=1).mean())
nl = 500
auto_cor = process_xcors(mat_total, max_len=nl)
mat_total = None
auto_cor = np.vstack(auto_cor)
np.save('length_filtered_auto_cor_tfs', auto_cor)
# print(len(auto_cor))
# clusters, sp_obj = cluster_mats(auto_cor[:,:500],res=0.40,neighbors=15)
# #plotter(mat_total, auto_cor, clusters, 5)