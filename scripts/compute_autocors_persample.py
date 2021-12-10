#Goal: take processed IPD HMM pickles and compute per molecule autocorrelograms out to 500 nucleotides (filter out molecules < 1 kb)

import os,sys,re
import numpy as np
import scipy as sp
import pickle

def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def ccf(x, y):
    result = np.correlate(y - np.mean(y), x - np.mean(x), 'same') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    return result[length:]

def xcor(a,b, maxlen = 1000, lengths=False):
    foo = []
    for i in range(len(a)):
        res = ccf(a[i],b[i])
        foo.append(res)
    return foo

def process_xcors(tipds, min_len=1000, lengths=False):
    cors = []
    zmws = []
    for zmw in tipds:
        read = np.nan_to_num(tipds[zmw]) >= 0.5
        if len(read) < min_len: continue
        cor = ccf(read, read)[:500]
        cors.append(cor)
        zmws.append(zmw)
    return (np.vstack(cors), zmws)

def main():
    hmm_pickle = sys.argv[1]
    binary = eat_pickle_binary(hmm_pickle)
    cors, zmws = process_xcors(binary)
    np.save('%s.autocors.npy' % hmm_pickle, cors)
    print("File written.")
    fho = open('%s.zmw_ids.txt' % hmm_pickle, 'w')
    for zmw in zmws:
        print(zmw, file = fho)
    print('ZMWs written.')
    fho.close()

if __name__ == "__main__":
    main()
