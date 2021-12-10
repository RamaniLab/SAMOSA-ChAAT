#The goal of this script is to draw a random set of molecules given a matrix of molecules with known cluster defintions and 
#arrange them in a text file for easy ggplotting.
import os,sys,re
import numpy as np
import pandas as pd
import scipy as sp
import pickle
import scanpy
from itertools import groupby
from Bio import Seq, SeqIO

import os,sys,re
import numpy as np
import pandas as pd
import pickle
import scanpy
from itertools import groupby
from Bio import Seq, SeqIO

def ccf(x, y):
    result = np.correlate(y - np.mean(y), x - np.mean(x), 'same') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    return result[length:]

def xcor(a,b, maxlen = 1000, lengths=False):
    foo = sp.signal.correlate(a,b, mode='same', method='fft') / (np.std(x, axis = 1) * np.std(y, axis = 1) * len(a[0]))
    length = (len(foo[0]) - 1) // 2
    return foo[:,length:]

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

def parseSampleLabels(string):
    split = string.split('_')
    #the initial replicates are all single-turnover
    if split[0] != 'CTCF':
        to = 'Single-turnover'
    footprinting_remodeling = {'plusM':('Footprinted','Native'), 'plusATP':('Footprinted','(+)ATP'),'plusADP':
                               ('Footprinted','(+)ADP'), 'minusATP':('Footprinted','(-)ATP'), 'minusM':('Negative','Negative')}
    density = {'15to1':'15:1','10to1':'10:1','5to1':'5:1','20to1':'20:1','naked':'Negative'}
    turnover = {'ST':'Single-turnover','MT':'Multi-turnover'}
    sequence_id = {'Ind':'S1','IndepArray':'S1','Indep':'S1','Dep':'S2','DepenArray':'S2','Depen':'S2',
                   'Independent':'S1', 'Dependent':'S2'}
    for i in split:
        if i in footprinting_remodeling:
            fp, cond = footprinting_remodeling[i]
        if i in density:
            dens = density[i]
        if i in turnover:
            to = turnover[i]
        if i in sequence_id:
            sid = sequence_id[i]
    return (fp, cond, dens, to, sid)

def findFootprints(tipds, valid_zmws=False):
    regions = {'zmw':[], 'length':[], 'start':[], 'end':[]}
    if not valid_zmws:
        valid_zmws = []
        for zmw in tipds:
            valid_zmws.append(zmw)
    for zmw in valid_zmws:
        hmm = tipds[zmw]
        inacregion = hmm
        inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5
        inacswitch = np.diff(inacregion)
        switchp = np.where(np.logical_or(inacswitch == 1, inacswitch == -1))[0]
        inInacReg = False
        regStart = -1
        regEnd = -1
        for point in switchp:
            if inacswitch[point] == -1 and not inInacReg:
                inInacReg = True
                regStart = point + 1
            if inacswitch[point] == 1 and inInacReg:
                inInacReg = False
                regEnd = point
                regions['zmw'].append(zmw)
                regions['length'].append(regEnd - regStart)
                regions['start'].append(regStart)
                regions['end'].append(regEnd)
    regionD = pd.DataFrame(regions)
    return(regionD)

def sample_mols_cluster(mols, cluster_filename, fho):
    cluster_fhi = open(cluster_filename)
    clusters = []
    for line in cluster_fhi:
        split = line.split('\t')
        cluster = split[1]
        clusters.append(cluster)
    clusters = np.array(clusters)
    for i in np.unique(clusters):
#       print(i)
        temp = mols[clusters == i]
        if len(temp) < 1000: continue
        samp = temp[np.random.choice(temp.shape[0], 1000, replace=False), :]
#        print(len(samp))
        idx = 0
        for mol in samp:
            for j in range(len(mol)):
                print("%s\t%s\t%s\t%s" % (idx, j, int(mol[j]),i), file=fho)
            idx+=1
    cluster_fhi.close()

def sample_mols_standard(mols, labels, fho):
    for i in np.unique(labels):
#       print(i)
        temp = mols[labels == i]
        if len(temp) < 1000: continue
        samp = temp[np.random.choice(temp.shape[0], 500, replace=False), :]
#        print(len(samp))
        idx = 0
        for mol in samp:
            for j in range(len(mol)):
                print("%s\t%s\t%s\t%s" % (idx, j, int(mol[j]),i), file=fho)
            idx+=1

def cluster_filter(cluster_df, lower, upper):
    valid_clusters = []
    for i in range(np.amax(cluster_df['clusters'])):
        clust_size = len(cluster_df[cluster_df['clusters'] == i])
        print(clust_size)
        if clust_size < 2500: 
            cutoff = i
            break
    cluster_df_filt = cluster_df[(cluster_df['clusters'] < cutoff) & (cluster_df['density'] > lower) & 
        (cluster_df['density'] <= upper)]
    print(len(cluster_df))
    print(len(cluster_df_filt))
    return cluster_df_filt        
        
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
        
    
##############################################
#####THIS IS ALL BOILERPLATE##################
#####CODE FOR PARSING OUR SAMPLE##############
#####NAMES.##################################
##############################################

#First, start with Snf2h independent site
refFile = '/avicenna/vramani/analyses/pacbio/pbrun10_CTCFpool_2/snf2h_independent_site_observed.fasta'
for ir, record in enumerate(SeqIO.parse(refFile, 'fasta')):
    if ir > 0:
        raise InputError('Reference fasta has multiple entries')
    irefseq = record.seq # reference sequence for independent CTCF site
#Next load up the Snf2h dependent site
refFile = '/avicenna/vramani/analyses/pacbio/pbrun10_CTCFpool_2/snf2h_dependent_site_observed.fasta'
for ir, record in enumerate(SeqIO.parse(refFile, 'fasta')):
    if ir > 0:
        raise InputError('Reference fasta has multiple entries')
    drefseq = record.seq # reference sequence for dependent CTCF site

#Load up the master list of all samples (CSV format)    
sampleRef1 = pd.read_csv('/avicenna/vramani/analyses/pacbio/pbrun10_CTCFpool_2/pbrun10_CTCFpool_2.sampleReference.csv')
sampleRef1['replicate'] = np.repeat('Rep1',len(sampleRef1['sampleName'].values))

sampleRef2 = pd.read_csv('/avicenna/vramani/analyses/pacbio/pbrun10_CTCFpool_1/pbrun10_CTCFpool_1.sampleReference.csv') 
sampleRef2['replicate'] = np.repeat('Rep1',len(sampleRef2['sampleName'].values))

sampleRef3 = pd.read_csv('/avicenna/vramani/analyses/pacbio/210516_NA_SNF2hCTCFarray_ST_rep2/210516_NA_SNF2hCTCFarray_ST_rep2.sampleReference.wynton.csv')
sampleRef3['replicate'] = np.repeat('Rep2',len(sampleRef3['sampleName'].values))

sampleRef4 = pd.read_csv('/avicenna/vramani/analyses/pacbio/210520_NA_SNF2hCTCFarray_MT_rep1/210520_NA_SNF2hCTCFarray_MT_rep1.sampleReference.wynton.csv')

ref4repname = []
for name in sampleRef4['sampleName']:
    if name.split('_')[-1] == 'rep1':
        ref4repname.append('Rep1')
    else:
        ref4repname.append('Rep2')
        
sampleRef4['replicate'] = ref4repname

sampleRef5 = pd.read_csv('/avicenna/vramani/analyses/pacbio/210608_NA_SNF2hCTCFarray_MT_rep2/210608_NA_SNF2hCTCFarray_MT_rep2.sampleReference.wynton.csv')
sampleRef5['replicate'] = np.repeat('Rep2',len(sampleRef5['sampleName'].values))

sampleRef = pd.concat([sampleRef1, sampleRef2, sampleRef3, sampleRef4, sampleRef5], ignore_index=True)

methylation_status = []
turnover_condition = []
sequence_id = []
catalytic_condition = []
density = []

for name in sampleRef['sampleName'].values:
    fp, cond, dens, to, sid = parseSampleLabels(name)
    methylation_status.append(fp)
    catalytic_condition.append(cond)
    turnover_condition.append(to)
    density.append(dens)
    sequence_id.append(sid)

sampleRef['meth_status'] = methylation_status
methylation_status = np.array(methylation_status)

sampleRef['catalytic'] = catalytic_condition  
catalytic_condition = np.array(catalytic_condition)

sampleRef['density'] = density
density = np.array(density)

sampleRef['turnover'] = turnover_condition
turnover_condition = np.array(turnover_condition)

sampleRef['seq'] = sequence_id
sequence_id = np.array(sequence_id)

sampleRef['sample_label'] = sampleRef[['density', 'catalytic', 'seq', 'turnover', 'replicate']].agg('_'.join, axis=1)

sampleRef.to_csv('all_samples_qc.txt',sep='\t')

####END BOILERPLATE#####

###############S1 Analyses###################
fho1 = open('S1_fishers_remodeling.data','w')
fho2 = open('S1_autocor_clusters.data','w')

cluster_file = open('S1_remodeled_leiden_autocor_labels.data')

subSample1 = sampleRef[(sampleRef['seq'] == 'S1') & 
    ((sampleRef['catalytic'] == 'Native') | (sampleRef['catalytic'] == '(+)ATP'))  & (sampleRef['density'] != 'Negative') & 
    (sampleRef['turnover'] == 'Single-turnover')]
cell_names = subSample1['cell'].values
samp_names = subSample1['sampleName'].values
samp_labels = subSample1['sample_label'].values

mols = []
labels = []

for i in range(len(cell_names)):
    cell_name = cell_names[i]
    samp_name = samp_names[i]
    samp_label = samp_labels[i]
    file_name = '/avicenna/vramani/analyses/pacbio/%s/processed/binarized/%s_%s_HMM.npy' % (cell_name, cell_name, samp_name)
    calls = np.load(file_name) > 0.5
    for mol in calls:
        crosscor = ccf(mol, mol)
        mols.append(crosscor)
        labels.append(samp_label)
mols = np.vstack(mols)
# crosscor = xcor(mols)
# print(crosscor)
labels = np.array(labels)
densities = []
labels_refined = []
clusters = []
for line in cluster_file:
    label, cluster, umap1, umap2, density = line.split()
    clusters.append(int(cluster))
    assembly, catalytic, fiber, turnover, replicate = label.split('_')
    labels_refined.append(catalytic)
    densities.append(int(density))

clusters = np.array(clusters)
clusters_df = pd.DataFrame({'clusters': clusters, 'catalytic': labels_refined, 'density': densities})
clusters_df_filt = cluster_filter(clusters_df, 1, 18)
clusters_df_filt_native = clusters_df_filt[clusters_df_filt['catalytic'] == 'Native']
clusters_df_filt_remodeled = clusters_df_filt[clusters_df_filt['catalytic'] == '(+)ATP']

fishers_test(clusters_df_filt_native, 'Native', fho1)
fishers_test(clusters_df_filt_remodeled, 'Remodeled', fho1)
fho1.close()

cluster_viz(clusters, mols, fho2)
fho2.close()

###############S2 Analyses###################

fho1 = open('S2_fishers_remodeling.data','w')
fho2 = open('S2_autocor_clusters.data','w')

cluster_file = open('S2_remodeled_leiden_autocor_labels.data')

subSample1 = sampleRef[(sampleRef['seq'] == 'S2') & 
    ((sampleRef['catalytic'] == 'Native') | (sampleRef['catalytic'] == '(+)ATP'))  & (sampleRef['density'] != 'Negative') & 
    (sampleRef['turnover'] == 'Single-turnover')]
cell_names = subSample1['cell'].values
samp_names = subSample1['sampleName'].values
samp_labels = subSample1['sample_label'].values

mols = []
labels = []

for i in range(len(cell_names)):
    cell_name = cell_names[i]
    samp_name = samp_names[i]
    samp_label = samp_labels[i]
    file_name = '/avicenna/vramani/analyses/pacbio/%s/processed/binarized/%s_%s_HMM.npy' % (cell_name, cell_name, samp_name)
    calls = np.load(file_name) > 0.5
    for mol in calls:
        crosscor = ccf(mol, mol)
        mols.append(crosscor)
        labels.append(samp_label)
mols = np.vstack(mols)
# crosscor = xcor(mols)
# print(crosscor)
labels = np.array(labels)
densities = []
labels_refined = []
clusters = []
for line in cluster_file:
    label, cluster, umap1, umap2, density = line.split()
    clusters.append(int(cluster))
    assembly, catalytic, fiber, turnover, replicate = label.split('_')
    labels_refined.append(catalytic)
    densities.append(int(density))

clusters = np.array(clusters)
clusters_df = pd.DataFrame({'clusters': clusters, 'catalytic': labels_refined, 'density': densities})
clusters_df_filt = cluster_filter(clusters_df, 1, 20)
clusters_df_filt_native = clusters_df_filt[clusters_df_filt['catalytic'] == 'Native']
clusters_df_filt_remodeled = clusters_df_filt[clusters_df_filt['catalytic'] == '(+)ATP']

fishers_test(clusters_df_filt_native, 'Native', fho1)
fishers_test(clusters_df_filt_remodeled, 'Remodeled', fho1)
fho1.close()

cluster_viz(clusters, mols, fho2)
fho2.close()






# fho1 = open('S2_fishers_remodeling.data','w')
# fho2 = open('S2_autocor_clusters.data','w')

# cluster_file = open('S2_remodeled_leiden_autocor_labels.data')

# subSample1 = sampleRef[(sampleRef['seq'] == 'S2') & 
#     ((sampleRef['catalytic'] == 'Native') | (sampleRef['catalytic'] == '(+)ATP'))  & (sampleRef['density'] != 'Negative') & 
#     (sampleRef['turnover'] == 'Single-turnover')]
# cell_names = subSample1['cell'].values
# samp_names = subSample1['sampleName'].values
# samp_labels = subSample1['sample_label'].values

# mols = []
# labels = []

# for i in range(len(cell_names)):
#     cell_name = cell_names[i]
#     samp_name = samp_names[i]
#     samp_label = samp_labels[i]
#     file_name = '/avicenna/vramani/analyses/pacbio/%s/processed/binarized/%s_%s_HMM.npy' % (cell_name, cell_name, samp_name)
#     calls = np.load(file_name) > 0.5
#     for mol in calls:
#         crosscor = ccf(mol, mol)
#         mols.append(crosscor)
#         labels.append(samp_label)
# mols = np.vstack(mols)
# labels = np.array(labels)
# densities = []
# labels_refined = []
# clusters = []
# for line in cluster_file:
#     label, cluster, umap1, umap2, density = line.split()
#     clusters.append(int(cluster))
#     assembly, catalytic, fiber, turnover, replicate = label.split('_')
#     labels_refined.append(catalytic)
#     densities.append(int(density))
# clusters = np.array(clusters)
# clusters_df = pd.DataFrame({'clusters':clusters, 'catalytic':labels_refined, 'density':densities})
# clusters_df_filt = cluster_filter(clusters_df, 1, 20)
# clusters_df_filt_native = clusters_df_filt[clusters_df_filt['catalytic'] == 'Native']
# clusters_df_filt_remodeled = clusters_df_filt[clusters_df_filt['catalytic'] == '(+)ATP']

# fishers_test(clusters_df_filt_native, 'Native', fho1)
# fishers_test(clusters_df_filt_remodeled, 'Remodeled', fho1)
# fho1.close()

# cluster_viz(clusters, mols, fho2)
# fho2.close()