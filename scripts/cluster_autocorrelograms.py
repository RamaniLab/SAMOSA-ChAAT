import os,sys,re
import scanpy
import numpy as np
import pandas as pd

def zmw_filter(valid_zmws, zmw_list, auto_cor, samp_label):
    '''We are only going to cluster ZMWs across our replicates that fall within the domains / regions we care about.'''
    labels = []
    samps = []
    bool_filter = []
    valid = {}
    zmws = []
    valid_fhi = open(valid_zmws)
    list_open = open(zmw_list)
    for line in valid_fhi:
        split = line.split()
        zmw = split[0]
        valid[zmw] = split[4]
    for line in list_open:
        split = line.split()
        zmw = split[0]
        if zmw in valid:
            bool_filter.append(True)
            labels.append(valid[zmw])
            samps.append(samp_label)
        else:
            bool_filter.append(False)
            labels.append('NA')
            samps.append('NA')
        zmws.append("%s_%s" % (samp_label, zmw))
    list_open.close()
    valid_fhi.close()
    bool_filter = np.array(bool_filter)
    labels = np.array(labels)
    samps = np.array(samps)
    zmws = np.array(zmws)
    submat = auto_cor[bool_filter]
    labels = labels[bool_filter]
    zmws  = zmws[bool_filter]
    samps = samps[bool_filter]
    return (submat, labels, zmws, samps)

    
def cluster_autocorrelograms(auto_cor):
    np.nan_to_num(auto_cor, copy=False)
    auto_cor = scanpy.AnnData(X=auto_cor)
    scanpy.tl.pca(auto_cor)
    scanpy.pp.neighbors(auto_cor, metric='correlation', n_neighbors=10)
    scanpy.tl.leiden(auto_cor, resolution=0.4)
    clusters = np.array(auto_cor.obs['leiden'])
    return(clusters)

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
    pwd_hmm = '/avicenna/vramani/analyses/pacbio/symlinks_for_hmm_calls/'
    pwd_aligned = '/avicenna/vramani/analyses/pacbio/symlinks_for_aligned_files/'
    
    filelist_npy = ['210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c1_E14_plusM_A_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c1_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c1_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c2_E14_plusM_A_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c2_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c2_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c3_E14_plusM_A_rep1_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c3_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.autocors.npy',
        '210830_NA_SNF2H-pool_c3_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.autocors.npy',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        'SAMv2_E14_mESCs_+m_preptest_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        'SAMv2_E14_mESCs_+m_preptest_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.autocors.npy',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy',
        'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.autocors.npy']
    
    filelist_txt = ['210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210421_NA_SAMv2_mESCs_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_1.0prep_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210427_NA_SAMv2mESCs_2.0prep_E14_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c1_E14_plusM_A_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c1_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c1_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c2_E14_plusM_A_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c2_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c2_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c3_E14_plusM_A_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c3_SNF2hKO_plusM_A7_NNsingle_HMM.pickle.zmw_ids.txt',
        '210830_NA_SNF2H-pool_c3_Snf2hWTAB_plusM_A3_NNsingle_HMM.pickle.zmw_ids.txt',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        'pbrun11_mESCs_SNF2h_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        'pbrun11_mESCs_SNF2h_SNF2hWTAB_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',                     
        'SAMv2_E14_mESCs_+m_preptest_E14_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        'SAMv2_E14_mESCs_+m_preptest_SNF2hKO_mESCs_plusM_rep2_NNsingle_HMM.pickle.zmw_ids.txt',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep_SNF2hKO_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt',
        'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep_SNF2hWTAB_mESCs_plusM_rep1_NNsingle_HMM.pickle.zmw_ids.txt']
 
    aligned_file_list = ['210421_NA_SAMv2_mESCs.split.E14_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        '210421_NA_SAMv2_mESCs.split.E14_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_1.0prep.split.SNF2hKO_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_1.0prep.split.SNF2hKO_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_1.0prep.split.SNF2hWTAB_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_1.0prep.split.SNF2hWTAB_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_2.0prep.split.E14_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        '210427_NA_SAMv2mESCs_2.0prep.split.E14_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c1.split.E14_plusM_A_rep1.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c1.split.SNF2hKO_plusM_A7.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c1.split.Snf2hWTAB_plusM_A3.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c2.split.E14_plusM_A_rep1.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c2.split.SNF2hKO_plusM_A7.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c2.split.Snf2hWTAB_plusM_A3.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c3.split.E14_plusM_A_rep1.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c3.split.SNF2hKO_plusM_A7.ccs.aligned.sorted.bam',
        '210830_NA_SNF2H-pool_c3.split.Snf2hWTAB_plusM_A3.ccs.aligned.sorted.bam',
        'pbrun11_mESCs_SNF2h.split.SNF2hKO_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        'pbrun11_mESCs_SNF2h.split.SNF2hKO_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        'pbrun11_mESCs_SNF2h.split.SNF2hWTAB_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        'pbrun11_mESCs_SNF2h.split.SNF2hWTAB_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        'SAMv2_E14_mESCs_+m_preptest.split.E14_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        'SAMv2_E14_mESCs_+m_preptest.split.SNF2hKO_mESCs_plusM_rep2.ccs.aligned.sorted.bam',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2.split.SNF2hKO_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        'SAMv2_SNF2hKO_mESCs_+m1_1.0prep.split.SNF2hKO_mESCs_plusM_rep1.ccs.aligned.sorted.bam',
        'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep.split.SNF2hWTAB_mESCs_plusM_rep1.ccs.aligned.sorted.bam']    

    file_labels = ['E14_rep1_1', 'E14_rep2_1', 'SNF2hKO_rep1_1', 'SNF2hKO_rep2_1', 'SNF2hWTAB_rep1_1', 'SNF2hWTAB_rep2_1', 
                   'E14_rep1_2', 'E14_rep2_2', 'E14_rep3_1', 'SNF2hKO_rep3_1', 'SNF2hWTAB_rep3_1', 'E14_rep3_2', 
                   'SNF2h_KO_rep3_2', 'SNF2hWTAB_rep3_2', 'E14_rep3_3', 'SNF2hKO_rep3_3', 'SNF2hWTAB_rep3_3', 
                   'SNF2hKO_rep1_2', 'SNF2hKO_rep2_2', 'SNF2hWTAB_rep1_2','SNF2hWTAB_rep2_2', 'E14_rep1_3', 
                   'SNF2hKO_rep2_3', 'SNF2hKO_rep1_3', 'SNF2hKO_rep1_4','SNF2hWTAB_rep1_3']

    epigenome_marks = ['H3K4me3_peaks.bed.gz_zmws', 'H3K4me1_peaks.bed.gz_zmws', 'H3K36me3_peaks.bed.gz_zmws',
                   'H3K27me3_peaks.bed.gz_zmws', 'H3K9me3_peaks.bed.gz_zmws', 'atac_deseq_res.peaks.closed.gz_zmws',
                    'atac_deseq_res.peaks.open.gz_zmws']

    majsat_prefix = 'blast_mm_cen_tel_'
    majsat_suffix = '.stats.txt.gz.satellite_zmws'

    file_label_lookup = zip(filelist_npy, filelist_txt, aligned_file_list, file_labels)
    zmw_tot = []
    labels_tot = []
    samps_tot = []
    arrs = []
    for a,b,c,d in file_label_lookup:
        domain_files = []
        file = pwd_hmm + a
        autocor = np.load(file)
        zmw_list = pwd_hmm + b
        for mark in epigenome_marks:
            domain_file = "%s/%s.%s" % (pwd_aligned, c, mark)
            domain_files.append(domain_file)
        domain_files.append("%s/%s%s%s" % (pwd_aligned, majsat_prefix, c, majsat_suffix)) 
        for filtered in domain_files:
            subarr, labels, zmws, samps = zmw_filter(filtered, zmw_list, autocor, d)
            arrs.append(subarr)
            labels_tot.append(labels)
            zmw_tot.append(zmws)
            samps_tot.append(samps)
    arrs = np.vstack(arrs)
    labels_tot = np.concatenate(labels_tot)
    zmw_tot = np.concatenate(zmw_tot)
    samps_tot = np.concatenate(samps_tot)
    clusters = cluster_autocorrelograms(arrs)
    #filter molecules (we'll say all clusters that collectively make up <5% of the data)
    cutoff = cluster_cutoff(clusters)
    fho1 = open('invivo_final_autocor_clusters.data', 'w')
    fho2 = open('invivo_final_autocor_averages.data', 'w')
    for i in range(len(zmw_tot)):
        if int(clusters[i]) > cutoff: continue #use cutoff to filter out molecules from one of the filtered clusters
        print("%s\t%s\t%s\t%s" % (zmw_tot[i], clusters[i], labels_tot[i], samps_tot[i]), file = fho1)
    for i in np.unique(clusters):
        if int(i) > cutoff: continue #use cutoff to filter out molecules from one of the filtered clusters
        avg = np.nanmean(arrs[clusters == i],axis=0)
        for idx in range(len(avg)):
            print("%s\t%s\t%s" % (idx, avg[idx], i), file = fho2)
    fho1.close()
    fho2.close()
            
      
if __name__ == "__main__":
    main()
    

