#The goal here is to calculate site exposure frequencies for every 15 mer across sequences
#S1 and S2, for every density
import os,sys,re
import numpy as np
import pandas as pd
import pickle
import scanpy
from itertools import groupby
from Bio import Seq, SeqIO

#log2(p_remodeled / p_unremodeled)

def ccf(x, y):
    result = np.correlate(y - np.mean(y), x - np.mean(x), 'same') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    return result[length:]

def process_xcors(arr):
    cors = []
    for mol in arr:
        cor = ccf(mol, mol)
        cors.append(cor)
    return (np.vstack(cors))

def site_acc(arr):
    site_acc = []
    cutoff = 0.9
    for midp in np.arange(8, len(arr[0]) - 8):
        site_acc.append(np.mean(arr[:,midp-8:midp+8],axis=1) > cutoff)
    site_acc = np.transpose(np.vstack(site_acc))
    return(site_acc)

def peak_cruncher(ac):
    peaks1 = []
    for i in range(len(ac)):
        mol1 = ac[i]
        peak1, data = signal.find_peaks(mol1, height = 0.10, width = 25)
        if len(peak1) == 0:
            ac1p = np.nan
        else:
            ac1p = int(peak1[0])
        peaks1.append(ac1p)
    peaks1 = np.array(peaks1)
    med1 = np.nanmedian(peaks1)
    mad1 = scipy.stats.median_absolute_deviation(peaks1,nan_policy='omit')
    mis1 = np.count_nonzero(np.isnan(peaks1))
    return (peaks1, med1, mad1,mis1)

def cluster_mols(arr, res=0.4):
    np.nan_to_num(arr, copy=False)
    auto_cor = scanpy.AnnData(X=arr)
    scanpy.tl.pca(auto_cor)
    scanpy.pp.neighbors(auto_cor, metric='correlation', n_neighbors=10)
    scanpy.tl.leiden(auto_cor, resolution=res)
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


#Iterate through each density, and calculate site accessibility for each position in each molecule
fho = open('S1_site_accesibility_new.data', 'w')
mols = []
labels = []
densities = []
for density in np.unique(sampleRef['density']):
    if density != 'Negative':
        subSample1 = sampleRef[(sampleRef['seq'] == 'S1') & (sampleRef['density'] == density) & ((sampleRef['catalytic'] == 'Native') | (sampleRef['catalytic'] == '(+)ATP'))]
        cell_names = subSample1['cell'].values
        samp_names = subSample1['sampleName'].values
        samp_labels = subSample1['sample_label'].values
        for i in range(len(cell_names)):
            cell_name = cell_names[i]
            samp_name = samp_names[i]
            samp_label = samp_labels[i]
            file_name = '/avicenna/vramani/analyses/pacbio/%s/processed/binarized/%s_%s_HMM.npy' % (cell_name, cell_name, samp_name)
            density_csv = pd.read_csv('/avicenna/vramani/analyses/pacbio/%s/processed/density/%s_%s_density.csv' % (cell_name, cell_name, samp_name))
            calls = np.load(file_name) > 0.5
            densities.append(density_csv['nucs'].values)
            if len(calls) != len(density_csv['nucs'].values):
                print("this happens for filename %s. calls length is %s; density length is %s" % (file_name, len(calls), len(density_csv['nucs'].values)))
            for mol in calls:
                mols.append(mol)
                labels.append(samp_label)
labels = np.array(labels)
mols = np.vstack(mols)
densities = np.concatenate(densities)
sarray = site_acc(mols)

#Calculate single-molecule NRLs & cluster
# xcor_total = process_xcors(mols)
# nrls, nrl_median, nrl_mad, nrl_miss = peak_cruncher(xcor_total)
# nrls = np.array(nrls)
# clusters, umap_x, umap_y = cluster_mols(xcor_total)

# fho2 = open('S1_nrls_clusters.data','w')
# for i in range(len(mols)):
#     print("%s\t%s\t%s\t%" % (nrls[i],clusters[i],densities[i], labels[i]), file = fho2)
# fho2.close()
# print("Missed NRLs = %s" % nrl_miss)


#Calculate site accesibility for each sample
for density in np.unique(densities):
    for label in np.unique(labels):
        subarray = sarray[(densities == density) & (labels == label)]
        subarray_acc = np.sum(subarray, axis=0)
        subarray_tot = len(subarray)
        for i in range(len(subarray_acc)):
            print("%s\t%s\t%s\t%s\t%s" % (i, label, density, subarray_acc[i], subarray_tot), file = fho) 
fho.close()

fho = open('S2_site_accesibility_new.data', 'w')
mols = []
labels = []
densities = []
for density in np.unique(sampleRef['density']):
    if density != 'Negative':
        subSample1 = sampleRef[(sampleRef['seq'] == 'S2') & (sampleRef['density'] == density) & ((sampleRef['catalytic'] == 'Native') | (sampleRef['catalytic'] == '(+)ATP'))]
        cell_names = subSample1['cell'].values
        samp_names = subSample1['sampleName'].values
        samp_labels = subSample1['sample_label'].values
        for i in range(len(cell_names)):
            cell_name = cell_names[i]
            samp_name = samp_names[i]
            samp_label = samp_labels[i]
            file_name = '/avicenna/vramani/analyses/pacbio/%s/processed/binarized/%s_%s_HMM.npy' % (cell_name, cell_name, samp_name)
            density_csv = pd.read_csv('/avicenna/vramani/analyses/pacbio/%s/processed/density/%s_%s_density.csv' % (cell_name, cell_name, samp_name))
            calls = np.load(file_name) > 0.5
            densities.append(density_csv['nucs'].values)
            if len(calls) != len(density_csv['nucs'].values):
                print("this happens for filename %s. calls length is %s; density length is %s" % (file_name, len(calls), len(density_csv['nucs'].values)))
            for mol in calls:
                mols.append(mol)
                labels.append(samp_label)
    
labels = np.array(labels)
mols = np.vstack(mols)
densities = np.concatenate(densities)
sarray = site_acc(mols)

#Calculate site accesibility for each sample
for density in np.unique(densities):
    for label in np.unique(labels):
        subarray = sarray[(densities == density) & (labels == label)]
        subarray_acc = np.sum(subarray, axis=0)
        subarray_tot = len(subarray)
        for i in range(len(subarray_acc)):
            print("%s\t%s\t%s\t%s\t%s" % (i, label, density, subarray_acc[i], subarray_tot), file = fho) 
fho.close()