import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import os
import sys
import socket

def reverse_complement(sequence):
    revbase = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    newSequence = ''
    for b in sequence[::-1]:
        newSequence += revbase[b]
    return newSequence

def complement(sequence):
    revbase = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    newSequence = ''
    for b in sequence:
        newSequence += revbase[b]
    return newSequence

def createNNinputs(samp, block, dataPBase, sampName, sampCell, zmwlimit=None):
    bases = ['A', 'C', 'G', 'T']
    checkbases = np.concatenate([np.arange(-3,0), np.arange(1,11)])
    contextColDic = {checkbases[x]:{bases[i]:(x * 4 + i) for i in range(4)} for x in np.arange(len(checkbases))}
    usepercs = np.arange(10,91,10)

    ipdArr = [] # the mean IPD at this Adenine
    zmwArr = [] # The ZMW hole number
    zmwPos = [] # The position of this adenine within the CCS
    sampleArr = [] # The sample number of this zmw
    numSubArr = [] # The number of subreads that contributed to this adenine IPD mean

    infile = os.path.join(dataPBase, sampCell, 'processed', 'full',
                          sampCell + '_' + sampName + '_full_block{:03d}.pickle'.format(block))
    with open(infile, 'rb') as fin:
        ipdfull = pickle.load(fin, encoding="latin1")
        
    zmwlist = list(ipdfull.keys())
    if zmwlimit != None:
        assert zmwlimit <= len(zmwlist), "There aren't that many ZMW in the sample!"
        zmwlist = zmwlist[0:zmwlimit]
        
    for zmw in tqdm(zmwlist, position=0, smoothing=0, desc='%s: part 1' % (sampName),
                    mininterval=10):
        ipdd = ipdfull[zmw]

        for pos in np.arange(3, len(ipdd['read'])-10):
            if ~np.isnan(ipdd['reverseM'][pos]) and ipdd['read'][pos] == 'A':
                ipdArr.append(ipdd['reverseM'][pos])
                zmwArr.append(zmw)
                zmwPos.append(pos)
                sampleArr.append(samp)
                numSubArr.append(ipdd['reverseNsub'][pos])
                
        for pos in np.arange(10, len(ipdd['read'])-3):
            if ~np.isnan(ipdd['forwardM'][pos]) and ipdd['read'][pos] == 'T':
                ipdArr.append(ipdd['forwardM'][pos])
                zmwArr.append(zmw)
                zmwPos.append(pos)
                sampleArr.append(samp)
                numSubArr.append(ipdd['forwardNSub'][pos])
                
    ipdArr = np.array(ipdArr, dtype=np.float32).reshape(-1,1)
    zmwArr = np.array(zmwArr, dtype=np.int32).reshape(-1,1)
    zmwPos = np.array(zmwPos, dtype=np.int32).reshape(-1,1)
    sampleArr = np.array(sampleArr, dtype=np.int16).reshape(-1,1)
    numSubArr = np.array(numSubArr, dtype=np.int16).reshape(-1,1)
    
    contextmat = np.full((ipdArr.shape[0], len(checkbases)*4), False, dtype=np.bool_)
    percsmat = np.full((ipdArr.shape[0], 4 * len(usepercs)), np.nan, dtype=np.float32)
    
    ic = 0

    for zmw in tqdm(zmwlist, position=0, smoothing=0, desc='%s: part 2' % (sampName),
                    mininterval=10):
        ipdd = ipdfull[zmw]
        
        for pos in np.arange(3, len(ipdd['read'])-10):
            if ~np.isnan(ipdd['reverseM'][pos]) and ipdd['read'][pos] == 'A':
                revcontext = ipdd['read'][pos-3:pos+11]
                for offs in checkbases:
                    contextmat[ic, contextColDic[offs][revcontext[3+offs]]] = True
                for i in range(4):
                    percsmat[ic,np.arange(len(usepercs)) + (i * len(usepercs))] = ipdd['percentiles'][bases[i]][usepercs]
                ic += 1

        for pos in np.arange(10, len(ipdd['read'])-3):
            if ~np.isnan(ipdd['forwardM'][pos]) and ipdd['read'][pos] == 'T':
                forcontext = reverse_complement(ipdd['read'][pos-10:pos+4])
                for offs in checkbases:
                    contextmat[ic, contextColDic[offs][forcontext[3+offs]]] = True
                for i in range(4):
                    percsmat[ic,np.arange(len(usepercs)) + (i * len(usepercs))] = ipdd['percentiles'][bases[i]][usepercs]
                ic += 1
    
    if not os.path.exists(dataPBase + '%s/processed/forNN' % (sampCell)):
        os.makedirs(dataPBase + '%s/processed/forNN' % (sampCell))
        

    np.savez(os.path.join(dataPBase, sampCell, 'processed', 'forNN',
                          sampCell + '_' + sampName + '_forNN_block{:03d}.npz'.format(block)), 
             contextmat = contextmat,
             ipdArr = ipdArr,
             contextColDic = contextColDic, 
             zmwArr = zmwArr, 
             zmwPos = zmwPos,
             sampleArr = sampleArr,
             percsmat = percsmat,
             numSubArr = numSubArr)

        
def main():
    samp = int(sys.argv[1])
    block = int(sys.argv[2])
    
    if 'biochem1' in socket.gethostname():
        dataPBase = '/avicenna/vramani/analyses/pacbio/'
        figPBase = '/avicenna/cmcnally/pbanalysis/'
    if 'titan' in socket.gethostname():
        dataPBase = '/data/users/goodarzilab/colin/results/pacbio/'
    if 'wynton' in socket.gethostname():
        dataPBase = '/wynton/group/goodarzilab/ramanilab/results/pacbio/'
    if 'rumi' in socket.gethostname():
        raise Exception('no pacbio results folder on rumi')
    
    sampleRef = pd.read_csv(dataPBase + 'sampleRef_K562_mESC.csv', sep=',', index_col=0)
    if 'titan' in socket.gethostname():
        sampleRef = pd.concat([sampleRef,
                           pd.read_csv(dataPBase + '210520_NA_K562Brdu_repeat/210520_NA_K562Brdu_repeat.sampleReference.wynton.csv',index_col=0),
                           pd.read_csv(dataPBase + '211014_MO_BrdU_invivo/211014_MO_BrdU_invivo.sampleReference.wynton.csv',index_col=0),
                           pd.read_csv(dataPBase + '210930_MO_E14_K562_BrdU/210930_MO_E14_K562_BrdU.sampleReference.wynton.csv',index_col=0)],
                          ignore_index=True)
    
    createNNinputs(samp, block, dataPBase, sampleRef['sampleName'][samp], sampleRef['cell'][samp])
    
    
if __name__ == "__main__":
    main()
    
    
