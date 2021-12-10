#The goal here is to have one script that loads in all of the fibers,
#clusters them (dependent & independent), saves the ZMW IDs, cluster labels, info
#and UMAP coordinates. Can also extend this one stop shop script to compute
#footprint sizes too for each molecule (probably worth doing).
import os,sys,re
import numpy as np
import pandas as pd
import pickle
import scanpy
from itertools import groupby
from Bio import Seq, SeqIO


def cluster_mols(arr):
    np.nan_to_num(arr, copy=False)
    auto_cor = scanpy.AnnData(X=arr)
    scanpy.tl.pca(auto_cor)
    scanpy.pp.neighbors(auto_cor, metric='correlation', n_neighbors=10)
    scanpy.tl.leiden(auto_cor, resolution=0.4)
    scanpy.tl.umap(auto_cor)
    clusters = np.array(auto_cor.obs['leiden'])
    um_x = auto_cor.obsm['X_umap'][:,0]
    um_y = auto_cor.obsm['X_umap'][:,1]
    return(clusters, um_x, um_y)

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

#Create master lists that are replicate-labeled for unremodeled S1 and S2 fibers, cluster 
#with leiden clustering and save the leiden clusters plus UMAP coordinates. Begin with unremodeled
#S1 fibers:
subSample1 = sampleRef[(sampleRef['seq'] == 'S1') & (sampleRef['catalytic'] == 'Native')  & (sampleRef['density'] != 'Negative')]
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
        mols.append(mol)
        labels.append(samp_label)
mols = np.vstack(mols)
clusters, um_x, um_y = cluster_mols(mols)

fho = open('S1_unremodeled_leiden_umap.data', 'w')
for i in range(len(clusters)):
    print("%s\t%s\t%s\t%s" % (labels[i], clusters[i], um_x[i], um_y[i]), file = fho)
fho.close()

#S2 fibers:
subSample1 = sampleRef[(sampleRef['seq'] == 'S2') & (sampleRef['catalytic'] == 'Native')  & (sampleRef['density'] != 'Negative')]
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
        mols.append(mol)
        labels.append(samp_label)
mols = np.vstack(mols)
clusters, um_x, um_y = cluster_mols(mols)

fho = open('S2_unremodeled_leiden_umap.data', 'w')
for i in range(len(clusters)):
    print("%s\t%s\t%s\t%s" % (labels[i], clusters[i], um_x[i], um_y[i]), file = fho)
fho.close()