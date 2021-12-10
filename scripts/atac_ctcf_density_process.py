#Goal of this script is to process the in vivo single-fiber data for SNF2h-KO and SNF2h-WTAB samples
#in the contexts of Snf2h-dependent ATAC-seq peaks & Ctcf sites called by ENCODE.
import os,sys,re
import pandas as pd
import pickle
import numpy as np
import scipy as sp

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

def distill_tipds_flat_density(tipds, sitelist, density, frac, length, label, subsample=0):
    sites = open(sitelist,'r')
    hole_nos = {}
    for i in sites:
        split = i.split()
        #check to make sure the feature is in the middle of the read 
        #and enough on both sides
        if int(split[2]) <= int(split[4]) <= int(split[3]):
            index1 = int(split[4]) - int(split[2])
            index2 = int(split[3]) - int(split[4])
            strand = split[-1]
            r_strand = split[-2]
            if index1 > (length / 2) and index2 > (length / 2):
                hole_nos[split[0]] = (index1, strand, r_strand, (split[1],split[4]))
    sites.close()
    lengths = []
    labs = []
    reads = []
    densities = []
    fracs = []
    sites = []
    mno = 0
    for hole in hole_nos:
        if int(hole) not in tipds: 
#           print("This happens\t%s" % hole)
            continue
        if int(hole) not in density: continue
        hole_dens = density[int(hole)]
        hole_frac = frac[int(hole)]
        read = tipds[int(hole)][12:-12]
        index,strand,r_strand, site_info = hole_nos[hole]
        if r_strand == '-':
            read = read[::-1]
        if strand == "+":
            extract = tipds[int(hole)][int(index - (length / 2)): int(index + (length / 2))]
            if len(extract) != length: continue
            reads.append(extract)
            labs.append(label)
            lengths.append(len(tipds[int(hole)]))
            densities.append(hole_dens)
            fracs.append(hole_frac)
            sites.append(site_info)
        else:
            extract = tipds[int(hole)][::-1][int(index - (length / 2)):int(index + (length / 2))]
            if len(extract) != length: continue
            reads.append(extract)
            labs.append(label)
            lengths.append(len(tipds[int(hole)]))
            densities.append(hole_dens)
            fracs.append(hole_frac)
            sites.append(site_info)
        if subsample != 0:
            mno += 1
            if mno == subsample: break
    new_mat = np.vstack(reads)
    return (new_mat, np.array(lengths), np.array(labs), np.array(densities), np.array(fracs), sites)

def main():
    pwd_hmm = '/avicenna/vramani/analyses/pacbio/symlinks_for_hmm_calls/'
    pwd_aligned = '/avicenna/vramani/analyses/pacbio/symlinks_for_aligned_files/'
    pwd_density = '/avicenna/vramani/analyses/pacbio/symlinks_for_hmm_calls/density_links/'
        
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
        '210830_NA_SNF2H-pool_c1.split.E14_plusM_A_rep1.ccs.aligned.sorted.bam',
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
    
    zipped = zip(filelist_npy, density_file_list, aligned_file_list, file_labels)
    length = 500
    factors = ['atac_deseq_res.peaks.open.flat',
               'atac_deseq_res.peaks.closed.flat',
               'CTCF_mES.sites.flat.sorted.distances.filtered',
               'mm10.archetype_motifs.CTCF.v1.0.bed.sampled.sorted.flat']    
    mats = []
    labels_tot = []
    densities_tot = []
    fracs_tot = []
    sites_tot = []
    for signal, density, aligned, rep_label in zipped:
        signal_open = "%s/%s" % (pwd_hmm, signal)
        rep = eat_pickle_binary(signal_open)
        density_open = "%s/%s" % (pwd_density, density) 
        frac, dens = density2dict(density_open)
        for factor in factors:
            site_open = "%s/%s.%s.5kb.zmws" % (pwd_aligned, factor, aligned)
            mat_rep1, lengths, labels, densities, fracs, sites = distill_tipds_flat_density(rep, site_open, dens, frac, length, "%s_%s" % (rep_label, factor))
            mats.append(mat_rep1)
            densities_tot.append(densities)
            fracs_tot.append(fracs)
            labels_tot.append(labels)
            sites_tot.append(sites)
        rep = 0
    fracs_tot = np.concatenate(fracs_tot)
    densities_tot = np.concatenate(densities_tot)
    labels_tot = np.concatenate(labels_tot)
    sites_tot = np.concatenate(sites_tot)
    #print(sites_tot)
    mats = np.vstack(mats) > 0.5
    chrs = []
    sites = []
    for chrid, site in sites_tot:
        chrs.append(chrid)
        sites.append(site)
    #Save densities, fractions, and labels as dataframe
    mdata = pd.DataFrame({'frac_access': fracs_tot, 'density_est': densities_tot, 'labels_agg':labels_tot, 
                          'chrid':chrs, 'sites':sites})
    mdata.to_csv('KOvWT_access_ctcf_final.data',sep=',')
    mdata_labels = pd.DataFrame(columns=['sample','bio_rep','tech_rep','factor'])
    mdata_labels[['sample','bio_rep','tech_rep','factor']] = mdata['labels_agg'].str.split('_',expand=True,n=3)
    mdata_labels.to_csv('KOvWT_access_ctcf_sep_labels.data',sep=',')
    #Save mols matrix for plotting, clustering, etc.
    np.save('ATAC_Ctcf_mols_%sbp.npy' % length, mats)
    
if __name__ == "__main__":
    main()