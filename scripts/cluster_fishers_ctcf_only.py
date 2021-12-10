import os,sys,re
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd
import scanpy

def sample_mols_cluster(mols, clusters, fho):
    for i in np.unique(clusters):
        temp = mols[clusters == i]
        if len(temp) < 1000: continue
        samp = temp[np.random.choice(temp.shape[0], 1000, replace=False), :]
        idx = 0
        for mol in samp:
            for j in range(len(mol)):
                print("%s\t%s\t%s\t%s" % (idx, j, int(mol[j]),i), file=fho)
            idx+=1

def cluster_cutoff(clusters):
    tot = len(clusters)
    cum = 0
    for i in range(len(np.unique(clusters))):
        frac = len(clusters[clusters == str(i)]) / tot
        cum += frac
        if cum >= 0.95:
            cutoff = i
            break
    return cutoff

def main():
    pwd = '/avicenna/vramani/analyses/pacbio/'
    ctcf_sites = np.load(pwd + 'ATAC_Ctcf_mols_500bp.npy')
    labs_full = pd.read_csv(pwd + 'KOvWT_access_ctcf_final.data')
    labs_short = pd.read_csv(pwd + 'KOvWT_access_ctcf_sep_labels.data')


    ctcf_sub = ctcf_sites[((labs_short['sample'] == 'SNF2hWTAB') | (labs_short['sample'] == 'SNF2hKO')) & (labs_short['factor'] == 'CTCF_mES.sites.flat.sorted.distances.filtered')]

    labs_full_sub = labs_full[((labs_short['sample'] == 'SNF2hWTAB') | (labs_short['sample'] == 'SNF2hKO')) & (labs_short['factor'] == 'CTCF_mES.sites.flat.sorted.distances.filtered')]

    labs_short_sub = labs_short[((labs_short['sample'] == 'SNF2hWTAB') | (labs_short['sample'] == 'SNF2hKO')) & (labs_short['factor'] == 'CTCF_mES.sites.flat.sorted.distances.filtered')]


    clust = scanpy.AnnData(X=np.nan_to_num(ctcf_sub))
    scanpy.tl.pca(clust)
    scanpy.pp.neighbors(clust,metric='correlation',n_neighbors=15)
    scanpy.tl.leiden(clust,resolution=0.25)
    clusters = np.array(clust.obs['leiden'])
    cut = cluster_cutoff(clusters)
    #print(len(np.unique(clusters)))

    fho = open('cluster_averages_CTCF_111021.txt','w')
    fho_samp = open('cluster_samps_CTCF_111021.txt','w')

    for i in np.unique(clusters):
        if int(i) > cut: continue
        process = np.nanmean(ctcf_sub[clusters == i], axis=0)
        for j in range(len(process)):
            print("%s\t%s\t%s" % (i, j - 250, process[j]), file=fho)

    fho.close()
    sample_mols_cluster(ctcf_sub, clusters, fho_samp)
    fho_samp.close()
    
    fho_glob = open('CTCF_fishers_allmols.data','w')
    fho_glob_rep = open('CTCF_fishers_byrep.data','w')
    labs_short_sub['cluster_label'] = clusters
    for i in np.unique(labs_short_sub['sample']):
        for j in np.unique(labs_short_sub['cluster_label']):
            subSamp = labs_short_sub
            num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['sample'] == i)].values)
            num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['sample'] != i)].values)
            num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['sample'] == i)].values)
            num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['sample'] != i)].values)
            odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                [num_clust_notlab, num_notclust_notlab]])
            print("%s\t%s\tOverall\t%s\t%s\t%s" % (i, j, num_clust_lab, odds_r, pval), file=fho_glob)
            for rep in np.unique(labs_short_sub['bio_rep']):
                subSamp = labs_short_sub[(labs_short_sub['bio_rep'] == rep)]
                num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['sample'] == i)].values)
                num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['sample'] != i)].values)
                num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['sample'] == i)].values)
                num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['sample'] != i)].values)
                odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                    [num_clust_notlab, num_notclust_notlab]])
                print("%s\t%s\t%s\tOverall\t%s\t%s\t%s" % (i, j, rep, num_clust_lab, odds_r, pval), file=fho_glob_rep)
    fho_glob.close()
    fho_glob_rep.close()
    
    labs_full_sub['clusters'] = clusters
    labs_full_sub.to_csv('CTCF_mols_sites_and_clusters.csv')

if __name__ == "__main__":
    main()
