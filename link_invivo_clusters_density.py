import os,sys,re
import scanpy
import numpy as np
import pandas as pd
import scipy as sp
import gzip as gz
import pickle

def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def density2dict(filename):
    pd_dens = pd.read_csv(filename,sep=',',warn_bad_lines=True)
    pd_dens['molecule_length'] = pd.to_numeric(pd_dens['moleculeLength'], errors='coerce')
    pd_dens.dropna(inplace=True)
    pd_dens['nucDens'] = (pd_dens['nucs'].astype(float) / pd_dens['moleculeLength'].astype(float)) * 1000
    lookup1 = dict(zip(pd_dens['zmw'], pd_dens['fracInaccessible']))
    lookup2 = dict(zip(pd_dens['zmw'], pd_dens['nucDens']))
    return (lookup1, lookup2)

def main():
    pwd_hmm = '/avicenna/vramani/analyses/pacbio/symlinks_for_hmm_calls/'
    pwd_aligned = '/avicenna/vramani/analyses/pacbio/symlinks_for_aligned_files/'
    pwd_density = '/avicenna/vramani/analyses/pacbio/symlinks_for_hmm_calls/density_links/'

    total_dataframe = pd.read_csv('./invivo_final_autocor_clusters.fix.data', sep='\t', 
    header=None, names=['mol_id', 'cluster_label','chromatin_type','sample_id'])
    total_dataframe[['samp_label','bio_rep','tech_rep']] = total_dataframe['sample_id'].str.split('_', expand=True,n=3) 

    file_labels = ['E14_rep1_1', 'E14_rep2_1','SNF2hKO_rep1_1','SNF2hKO_rep2_1',
                   'SNF2hWTAB_rep1_1','SNF2hWTAB_rep2_1','E14_rep1_2','E14_rep2_2',
                   'E14_rep3_1','SNF2hKO_rep3_1','SNF2hWTAB_rep3_1',
                   'E14_rep3_2','SNF2hKO_rep3_2','SNF2hWTAB_rep3_2',
                   'E14_rep3_3','SNF2hKO_rep3_3','SNF2hWTAB_rep3_3',
                   'SNF2hKO_rep1_2','SNF2hKO_rep2_2','SNF2hWTAB_rep1_2','SNF2hWTAB_rep2_2','E14_rep1_3',
                   'SNF2hKO_rep2_3','SNF2hKO_rep1_3','SNF2hKO_rep1_4','SNF2hWTAB_rep1_3']

    filelist_npy = ['210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        '210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c1_E14_plusM_A_rep1_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c1_SNF2hKO_plusM_A7_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c1_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c2_E14_plusM_A_rep1_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c2_SNF2hKO_plusM_A7_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c2_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c3_E14_plusM_A_rep1_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c3_SNF2hKO_plusM_A7_NNsingle_HMM.pickle',
        '210830_NA_SNF2H-pool_c3_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        'SAMv2_E14_mESCs_+m_preptest_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        'SAMv2_E14_mESCs_+m_preptest_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle',
        'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle']

    density_file_list = ['210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep1_density.csv',
        '210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep2_density.csv',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep1_density.csv',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep2_density.csv',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_density.csv',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep2_density.csv',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep1_density.csv',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep2_density.csv',
        '210830_NA_SNF2H-pool_c1_E14_plusM_A_rep1_density.csv',
        '210830_NA_SNF2H-pool_c1_SNF2hKO_plusM_A7_density.csv',
        '210830_NA_SNF2H-pool_c1_Snf2hWTAB_plusM_A3_density.csv',
        '210830_NA_SNF2H-pool_c2_E14_plusM_A_rep1_density.csv',
        '210830_NA_SNF2H-pool_c2_SNF2hKO_plusM_A7_density.csv',
        '210830_NA_SNF2H-pool_c2_Snf2hWTAB_plusM_A3_density.csv',
        '210830_NA_SNF2H-pool_c3_E14_plusM_A_rep1_density.csv',
        '210830_NA_SNF2H-pool_c3_SNF2hKO_plusM_A7_density.csv',
        '210830_NA_SNF2H-pool_c3_Snf2hWTAB_plusM_A3_density.csv',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep1_density.csv',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep2_density.csv',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep1_density.csv',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep2_density.csv',
        'SAMv2_E14_mESCs_+m_preptest_E14_mESCs_plusM_rep1_density.csv',                         
        'SAMv2_E14_mESCs_+m_preptest_SNF2hKO_mESCs_plusM_rep2_density.csv',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2_SNF2hKO_mESCs_plusM_rep1_density.csv',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep_SNF2hKO_mESCs_plusM_rep1_density.csv',
        'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_density.csv']  

    clusters = {}
    clusters_density = {}
    zipped = zip(filelist_npy, density_file_list, file_labels)

    for signal, density, rep_label in zipped:
        signal_open = "%s/%s" % (pwd_hmm, signal)
        rep = eat_pickle_binary(signal_open)
        density_open = "%s/%s" % (pwd_density, density) 
        frac, dens = density2dict(density_open)
        subSamp = total_dataframe[total_dataframe['sample_id'] == rep_label]
        for cluster in np.unique(subSamp['cluster_label']):
            if cluster not in clusters:
                clusters[cluster] = []
                clusters_density[cluster] = []
            subSamp_cluster = subSamp[subSamp['cluster_label'] == cluster]
            for mol_id in subSamp_cluster['mol_id'].values:
                zmw = int(mol_id.split('_')[-1])
                if zmw not in rep: continue
                read = np.nan_to_num(rep[zmw])[:1000] >= 0.5
                if len(read) < 1000: continue
                clusters[cluster].append(read)
                clusters_density[cluster].append((dens[zmw], rep_label))
    fho_lp = open('in_vivo_cluster_sig_averages_102821.data', 'w')
    fho_dens = open('in_vivo_cluster_densities_102821.data','w')
    for cluster in clusters:
        clusters[cluster] = np.vstack(clusters[cluster])
        cluster_mean = np.nanmean(clusters[cluster],axis=0)
        for i in range(len(cluster_mean)):
            print("%s\t%s\t%s" % (i, cluster_mean[i], cluster), file=fho_lp)
        for d, lab in clusters_density[cluster]:
            print("%s\t%s\t%s" % (cluster, d, lab), file=fho_dens)
    fho_lp.close()
    fho_dens.close()       
            
        
if __name__ == "__main__":
    main()