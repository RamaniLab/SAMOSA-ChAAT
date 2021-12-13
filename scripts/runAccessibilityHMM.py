import numpy as np
import pandas as pd
import socket
import sys
import os
import pickle
from pomegranate import HiddenMarkovModel, State, BernoulliDistribution

if 'biochem1' in socket.gethostname():
    dataPBase = '/avicenna/vramani/analyses/pacbio/'
    figPBase = '/avicenna/cmcnally/pbanalysis/'
elif 'titan' in socket.gethostname():
    dataPBase = '/data/users/goodarzilab/colin/results/pacbio/'
elif 'wynton' in socket.gethostname():
    dataPBase = '/wynton/group/goodarzilab/ramanilab/results/pacbio/'
elif 'rumi' in socket.gethostname():
    raise Exception('no pacbio results folder on rumi')
else:
    # If the name doesn't match the above, assume it's a compute node of wynton, which have all sorts of different names
    # If that's an incorrect assumption the script will crash
    dataPBase = '/wynton/group/goodarzilab/ramanilab/results/pacbio/'
    

samp = int(sys.argv[1])
pieceN = int(sys.argv[2])

if len(sys.argv) > 3:
    if sys.argv[3].lower() == 'brdu':
        sampleRef = pd.concat([pd.read_csv(dataPBase + '210520_NA_K562Brdu_repeat/210520_NA_K562Brdu_repeat.sampleReference.wynton.csv',index_col=0),
                               pd.read_csv(dataPBase + '210930_MO_E14_K562_BrdU/210930_MO_E14_K562_BrdU.sampleReference.wynton.csv',index_col=0),
                               pd.read_csv(dataPBase + '211014_MO_BrdU_invivo/211014_MO_BrdU_invivo.sampleReference.wynton.csv',index_col=0)],
                              ignore_index=True)
    else:
        raise ValueError("brdu is the only 3rd argument currently recognized")
else:
    sampleRef = pd.read_csv(dataPBase + 'sampleRef_K562_mESC.csv', sep=',')

with open(dataPBase + '{0}/processed/forHMM/{0}_{1}_forHMM_piece{2:05d}.pickle'.format(sampleRef['cell'][samp],
                                                                                       sampleRef['sampleName'][samp],
                                                                                       pieceN), 'rb') as fin:
    hmmInput = pickle.load(fin)
    
goodZMW = list(hmmInput.keys())
hmmRes = {}

for zmw in goodZMW:
    goodDat = hmmInput[zmw]['inDat']
    cclen = hmmInput[zmw]['cclen']
    
    if goodDat.shape[0] < 2:
        continue
        # Don't try to create HMM if the input data is for some reason super short
        
    # Construct the HMM model for this molecule
    AcAd = []
    InacAd = []

    for row in goodDat.itertuples(index=True):
        AcAd.append(State(BernoulliDistribution(row.posProb), name="Ac_{0}".format(row.Index)))
        InacAd.append(State(BernoulliDistribution(row.negProb), name="Inac_{0}".format(row.Index)))

    model = HiddenMarkovModel()
    model.add_states(AcAd)
    model.add_states(InacAd)
    model.add_transition(model.start, AcAd[0], 0.5)
    model.add_transition(model.start, InacAd[0], 0.5)

    leaveInacProb = 1/1000
    leaveAcProb = 1/1000
    for b in np.arange(goodDat.shape[0]-1):
        dist = goodDat['pos'][b+1] - goodDat['pos'][b]
        stayInacP = (1 - leaveInacProb)**dist
        model.add_transition(InacAd[b], InacAd[b+1], stayInacP)
        model.add_transition(InacAd[b], AcAd[b+1], 1 - stayInacP)
        stayAcP = (1 - leaveAcProb)**dist
        model.add_transition(AcAd[b], AcAd[b+1], stayAcP)
        model.add_transition(AcAd[b], InacAd[b+1], 1 - stayAcP)

    model.add_transition(AcAd[-1], model.end, 1)
    model.add_transition(InacAd[-1], model.end, 1)
    model.bake()

    # Use viterbi to find most likely path based on the observed methylation
    path = model.viterbi(goodDat['methPred'])
    refbases = np.arange(cclen)
    pathRes = np.full(goodDat.shape[0], np.nan)
    for p in path[1]:
        psplit = p[1].name.split('_')
        if len(psplit) > 1:
            if psplit[0] == 'Ac':
                pathRes[int(psplit[1])] = 1
            if psplit[0] == 'Inac':
                pathRes[int(psplit[1])] = 0

    hmmRes[zmw] = np.interp(refbases, goodDat['pos'], pathRes, left=np.nan, right=np.nan)
    
if not os.path.exists(dataPBase + '%s/processed/binarized/HMMout' % (sampleRef['cell'][samp])):
    os.makedirs(dataPBase + '%s/processed/binarized/HMMout' % (sampleRef['cell'][samp]))
    
with open(dataPBase + '{0}/processed/binarized/HMMout/{0}_{1}_HMMres_piece{2:05d}.pickle'.format(sampleRef['cell'][samp],
                                                                                                 sampleRef['sampleName'][samp],
                                                                                                 pieceN), 'wb') as fout:
    pickle.dump(hmmRes, fout, protocol=4)
