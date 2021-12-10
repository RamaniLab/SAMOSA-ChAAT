import os,sys,re
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd

def fishers_test(cluster_df, label, fho):
#    cluster_df = pd.DataFrame({'labels' : labels_refined, 'clusters':clusters})
    tot = len(cluster_df)
    for i in np.unique(cluster_df['clusters']):
        for j in np.unique(cluster_df['density']):
            num_clust_lab = len(cluster_df[(cluster_df['clusters'] == i) & (cluster_df['density'] == j)])
            num_clust_notlab = len(cluster_df[(cluster_df['clusters'] == i) & (cluster_df['density'] != j)])
            num_notclust_lab = len(cluster_df[(cluster_df['clusters'] != i) & (cluster_df['density'] == j)])
            num_notclust_notlab = len(cluster_df[(cluster_df['clusters'] != i) & (cluster_df['density'] != j)])
            odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                [num_clust_notlab, num_notclust_notlab]])
            frac = num_clust_lab / tot
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i, j, num_clust_lab, odds_r, pval, label, frac), file=fho)
    
def cluster_viz(clusters, mols, fho):
    for i in np.unique(clusters):
        submols = mols[clusters == i]
        if len(submols) < 2500: continue
        submols_average = np.nanmean(submols, axis = 0)
        for j in range(len(submols_average)):
            print("%s\t%s\t%s" % (i, j, submols_average[j]), file = fho)
            
            
def main():
    total_dataframe = pd.read_csv('./invivo_final_autocor_clusters.fix.data', sep='\t', 
        header=None, names=['mol_id', 'cluster_label','chromatin_type','sample_id'])
    total_dataframe[['samp_label','bio_rep','tech_rep']] = total_dataframe['sample_id'].str.split('_', expand=True,n=3) 
    print(total_dataframe)
    fho = open('fishers_tests_smarca5_invivo_bydomain.txt','w')
    fho_rep = open('fishers_tests_smarca5_invivo_bydomain_replicates.txt','w')
    fho_glob = open('fishers_tests_smarca5_invivo_across_samp.txt','w')
    fho_glob_rep = open('fishers_tests_smarca5_invivo_across_samp_replicates.txt','w')
    
    for i in np.unique(total_dataframe['samp_label']):
        for j in np.unique(total_dataframe['cluster_label']):
            subSamp = total_dataframe
            num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['samp_label'] == i)].values)
            num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['samp_label'] != i)].values)
            num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['samp_label'] == i)].values)
            num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['samp_label'] != i)].values)
            odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                [num_clust_notlab, num_notclust_notlab]])
            print("%s\t%s\tOverall\t%s\t%s\t%s" % (i, j, num_clust_lab, odds_r, pval), file=fho_glob)
            for rep in np.unique(total_dataframe['bio_rep']):
                subSamp = total_dataframe[(total_dataframe['bio_rep'] == rep)]
                num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['samp_label'] == i)].values)
                num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['samp_label'] != i)].values)
                num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['samp_label'] == i)].values)
                num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['samp_label'] != i)].values)
                odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                    [num_clust_notlab, num_notclust_notlab]])
                print("%s\t%s\t%s\tOverall\t%s\t%s\t%s" % (i, j, rep, num_clust_lab, odds_r, pval), file=fho_glob_rep)
            for k in np.unique(total_dataframe['chromatin_type']):
                subSamp = total_dataframe[total_dataframe['samp_label'] == i]
                num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['chromatin_type'] == k)].values)
                num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['chromatin_type'] != k)].values)
                num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['chromatin_type'] == k)].values)
                num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['chromatin_type'] != k)].values)
                odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                    [num_clust_notlab, num_notclust_notlab]])
                print("%s\t%s\t%s\t%s\t%s\t%s" % (i, j, k, num_clust_lab, odds_r, pval), file=fho)
                for rep in np.unique(total_dataframe['bio_rep']):
                    subSamp = total_dataframe[(total_dataframe['samp_label'] == i) & (total_dataframe['bio_rep'] == rep)]
                    num_clust_lab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['chromatin_type'] == k)].values)
                    num_clust_notlab = len(subSamp[(subSamp['cluster_label'] == j) & (subSamp['chromatin_type'] != k)].values)
                    num_notclust_lab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['chromatin_type'] == k)].values)
                    num_notclust_notlab = len(subSamp[(subSamp['cluster_label'] != j) & (subSamp['chromatin_type'] != k)].values)
                    odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                        [num_clust_notlab, num_notclust_notlab]])
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i, j, k, rep, num_clust_lab, odds_r, pval), file=fho_rep)


    fho.close()
    fho_rep.close()
    fho_glob.close()
    fho_glob_rep.close()
    
if __name__ == "__main__":
    main()
