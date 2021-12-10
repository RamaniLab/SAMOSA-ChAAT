'''
extractIPD.py
Colin McNally
2020/03/26

For a given sample to process, extract the IPD value for each base, and save the T IPDs and the binarized versions
'''


from __future__ import print_function, division
import sys, os
import pandas as pd
import numpy as np
import re
import argparse
import pbcore.io as pb
from pbcore.sequence import reverseComplement
from Bio import Seq, SeqIO
from tqdm import tqdm
import edlib
import pickle
#from sklearn.mixture import GaussianMixture
#from sklearn.exceptions import ConvergenceWarning
import multiprocessing as mp
import glob
import warnings
#import Queue - seems unused, also causes failure with Python 3
import time

parser = argparse.ArgumentParser(description="Extract IPD values from an input sample")
parser.add_argument('referenceFile', nargs='?', default='/avicenna/vramani/analyses/pacbio/pbrun3-9_SampleReference.csv',
                   help='The sample reference file from which to []')
parser.add_argument('sample', 
                    help='Either an integer index into the reference file, or a string identifier in the format [cell].[samplename]')
parser.add_argument('-o', '--outputlocation', help='Location to save the outputs')
parser.add_argument('-j', '--threads', type=int, help='Number of threads to use. Defaults to 1')
parser.add_argument('-q', '--quiet', action='store_true', help='Add to supress the progress bar')
parser.add_argument('-f', '--filter', type=int, help='Filter ZMW with less than this number of subreads')

usepercentiles = range(10,41)
baserc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
p = re.compile(r'\d+[=XDI]')



# Extract IPD for each molecule from an array sample, return full results in a numpy array
def extractIPDfullArray(cbamFile, bamFile, refFile, args):
    bam = pb.IndexedBamReader(bamFile)
    
    for ir, record in enumerate(SeqIO.parse(refFile, 'fasta')):
        if ir > 0:
            raise InputError('Reference fasta has multiple entries')
        refseq = record.seq

    #Find bases that are A or T
    refisa = [b == 'A' for b in refseq]
    refist = [b == 'T' for b in refseq]
    
    #find all ZMW with ccs that pass filter
    
    outlist = []
    
    validzmw = []
    cbam = pb.IndexedBamReader(cbamFile)
    for cc in cbam:
        validzmw.append(cc.HoleNumber)
        
    
    for zind, zmw in enumerate(tqdm(validzmw, position=0, leave=True,
                                    desc='Getting IPDs', smoothing=0.001, mininterval=10)):
        subrs = bam.readsByHoleNumber(zmw)
        
        if len(subrs) < 10:
            continue
            
        useSubreads = []
        for sr in subrs:
            uc = sr.unrolledCigar(orientation='read')
            matchbase = np.sum(uc == 7)
            if (matchbase / len(refseq)) >= 0.81 and (matchbase / len(uc)) >= 0.85:
                useSubreads.append(sr)
        
        if len(useSubreads) < 10:
            continue
        
        allipds = np.empty((len(useSubreads), len(refseq)))
        allipds.fill(np.nan)
        ###allpws = np.empty((len(useSubreads), len(refseq)))
        ###allpws.fill(np.nan)
        subrOrient = np.empty(len(useSubreads))

        for index, sr in enumerate(useSubreads):
            subrOrient[index] = sr.isForwardStrand
            #using subread base calls that are matches or mismatches
            #may want to consider ignoring mismatches in the future, especially for non-synethetic source DNA
            #insertions = 1, gaps = 2, leaving both out of downstream analysis
            usebases = sr.unrolledCigar(orientation="genomic") == 7
            ipds = sr.baseFeature('Ipd',aligned=True, orientation="genomic")
            ###pws = sr.baseFeature('PulseWidth',aligned=True, orientation="genomic")
            refpos = sr.referencePositions(aligned=True, orientation="genomic")
            #use the reference positions to index the per base IPD values into their position in the reference sequence
            allipds[index, refpos[usebases]] = ipds[usebases]
            ###allpws[index, refpos[usebases]] = pws[usebases]

        
        
        # I expect to see RuntimeWarnings in this block
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            allipds[allipds < 1] = 1
            allipds = np.log10(allipds)
        
            percentiles = np.percentile(allipds[~np.isnan(allipds)],np.arange(0,101), interpolation='linear')
        
            formean = np.nanmean(allipds[subrOrient == True,], axis=0)
            revmean = np.nanmean(allipds[subrOrient == False,], axis=0)
            
        outlist.append({'fmean':formean.astype(np.float32), 'rmean':revmean.astype(np.float32),
                        'percentiles':percentiles.astype(np.float32), 'nsubr':len(useSubreads), 'zmw':zmw})
        
        #if len(outlist) >= 10000:
        #    break
        
    return(outlist)

# Extract the IPD from each zmw passed in the zmwqueue, return a dictionary containing onlyT IPDs, the binarized versions,
# and other informatin to include in the zmwinfo file
def extractIPDfullgenomic(cbamfile, alncbamfile, sbamfile, zmwqueue, outqueue): #filebase, nit):
    with pb.IndexedBamReader(alncbamfile) as alncbam, pb.IndexedBamReader(sbamfile) as sbam,\
         pb.IndexedBamReader(cbamfile) as cbam:
        zmw = zmwqueue.get()
        while zmw is not None:
            goodzmw = True
            cc = cbam.readsByHoleNumber(zmw)[0]
            ccread = cc.read(aligned=False, orientation='native')
            if (len(ccread) < 100 or
                ccread.count('A') < 10 or
                ccread.count('C') < 10 or
                ccread.count('G') < 10 or
                ccread.count('T') < 10):
                goodzmw = False
                
            if goodzmw:
                subrs = sbam.readsByHoleNumber(zmw)
                if len(subrs) < 10:
                    goodzmw = False
            
            if goodzmw:
                ccalns = alncbam.readsByHoleNumber(zmw)

                zmres = {}

                if len(ccalns) > 0:
                    alnlen = np.array([ccal.readLength for ccal in ccalns])
                    usealn = np.where(alnlen == np.max(alnlen))[0][0]
                elif len(ccalns) == 1:
                    usealn = 0
                elif len(ccalns) == 0:
                    usealn = None

                allipds = np.empty((len(subrs), len(ccread)), dtype='float32')
                allipds.fill(np.nan)
                subrOrient = np.empty(len(subrs), dtype='bool')

                for index, sr in enumerate(subrs):
                    # Test if this subread aligns to the forward or reverse of the CCS
                    forwardSread = sr.read(aligned=False, orientation='native')
                    reverseSread = reverseComplement(forwardSread)
                    faln = edlib.align(forwardSread, ccread, mode='NW', task='path')
                    raln = edlib.align(reverseSread, ccread, mode='NW', task='path')
                    if faln['editDistance'] < raln['editDistance']:
                        subrOrient[index] = True
                        alndir = faln
                        useread = forwardSread
                    else:
                        subrOrient[index] = False
                        alndir = raln
                        useread = reverseSread

                    # Use the alignment information to extract IPD at each base that aligns to the CCS
                    origb = np.empty(len(useread), dtype=np.int16 )
                    origb.fill(np.nan)
                    ccb = np.empty(len(useread), dtype=np.int16)
                    ccb.fill(np.nan)
                    subI = 0
                    ccI = 0
                    for m in p.finditer(alndir['cigar']):
                        lg = int(m.group()[-len(m.group()):-1])
                        mtype = m.group()[-1]
                        if mtype == '=':
                            origb[subI:(subI + lg)] = range(subI, subI + lg)
                            ccb[subI:(subI + lg)] = range(ccI, ccI + lg)
                            subI += lg
                            ccI += lg
                        elif mtype == 'X':
                            subI += lg
                            ccI += lg
                        elif mtype == 'I':
                            subI += lg
                        elif mtype == 'D':
                            ccI += lg

                    ccb = ccb[~np.isnan(ccb)]
                    origb = origb[~np.isnan(origb)]
                    if not subrOrient[index]:
                        for i in range(len(origb)):
                            origb[i] = -1 - origb[i]

                    ipds = sr.baseFeature('Ipd',aligned=False, orientation="native")
                    allipds[index, ccb] = ipds[origb]

                readisb = {refb:np.where([b == refb for b in ccread])[0] for refb in ['A','C','G','T']}

                # Take the mean IPD at each position, after taking log
                with warnings.catch_warnings(): # ignoring warnings from taking the mean of columns that are all NaN
                    warnings.simplefilter("ignore", category=RuntimeWarning)

                    allipds[allipds < 1] = 1 # so minimum log(IPD) is 0
                    allipds = np.log10(allipds)

                    percentiles = {}
                    percentiles['all'] = np.percentile(allipds[~np.isnan(allipds)],np.arange(0,101),
                                                       interpolation='linear')
                    for base in ['A', 'C', 'G', 'T']:
                        ipdsAtTemplate = np.concatenate([allipds[subrOrient == True, :][:,readisb[reverseComplement(base)]],
                                                         allipds[subrOrient == False,:][:,readisb[base]]],
                                                        axis=None)
                        percentiles[base] = np.percentile(ipdsAtTemplate[~np.isnan(ipdsAtTemplate)],
                                                          np.arange(0,101),
                                                          interpolation='linear').astype(np.float32)

                    forwardMean = np.nanmean(allipds[subrOrient == True,:], axis=0)
                    forwardNcontrib = np.sum(~np.isnan(allipds[subrOrient == True,:]), axis=0)
                    reverseMean = np.nanmean(allipds[subrOrient == False,:], axis=0)
                    reverseNcontrib = np.sum(~np.isnan(allipds[subrOrient == False,:]), axis=0)

                #tonlyMean = np.empty(len(ccread), dtype='float32')
                #tonlyMean.fill(np.nan)
                #tonlyMean[readisb['T']] = forwardMean[readisb['T']]
                #tonlyMean[readisb['A']] = reverseMean[readisb['A']]
                #tonlyMean[readisb['C']] = forwardMean[readisb['C']]
                #tonlyMean[readisb['G']] = reverseMean[readisb['G']]

                # Save useful information about this molecule
                zmres['zmw'] = zmw
                zmres['cclen'] = len(ccread)
                zmres['nsubr'] = len(subrs)
                zmres['naln'] = len(ccalns)
                if usealn is not None:
                    zmres['chr'] = ccalns[usealn].referenceName
                    zmres['refStart'] = ccalns[usealn].referenceStart
                    zmres['refEnd'] = ccalns[usealn].referenceEnd
                    zmres['alnStart'] = ccalns[usealn].aStart
                    zmres['alnEnd'] = ccalns[usealn].aEnd
                else:
                    zmres['chr'] = "noAlignment"
                    zmres['refStart'] = -1
                    zmres['refEnd'] = -1
                    zmres['alnStart'] = -1
                    zmres['alnEnd'] = -1
                with warnings.catch_warnings(): # ignoring warnings from taking the mean of all NaN
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    # Something is wrong with the read if these throw a warning, but still no need to print on command line
                    # The stored value will be NaN
                    zmres['basemeanA'] = np.nanmean(np.concatenate([forwardMean[readisb['A']],
                                                                    reverseMean[readisb['T']]],axis=None))
                    zmres['basemeanC'] = np.nanmean(np.concatenate([forwardMean[readisb['C']],
                                                                    reverseMean[readisb['G']]],axis=None))
                    zmres['basemeanG'] = np.nanmean(np.concatenate([forwardMean[readisb['G']],
                                                                    reverseMean[readisb['C']]],axis=None))
                    zmres['basemeanT'] = np.nanmean(np.concatenate([forwardMean[readisb['T']],
                                                                    reverseMean[readisb['A']]],axis=None))

                zmres['full'] = {'forwardM':forwardMean.astype(np.float32),
                                 'forwardNSub':forwardNcontrib.astype(np.int16),
                                 'reverseM':reverseMean.astype(np.float32),
                                 'reverseNsub':reverseNcontrib.astype(np.int16),
                                 'read':ccread, 'percentiles':percentiles}

                outqueue.put(zmres) # put the results in the output queue
            # even if the zmw wasn't used, mark it as done and move on
            zmwqueue.task_done()
            zmw = zmwqueue.get()
        zmwqueue.task_done()
        outqueue.put(None)


# Collect genomic results for each molecule, aggregate them in dictionaries and a data frame and save the results periodically
def listenerSaver(zmwqueue, outqueue, nNoneQueue, chunkSize, outbase, sampcn): #, totalccs):   
    MipdDic = {}
    ccdicdat = {}
    for key in ['zmw', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'cclen', 'nsubr', 'basemeanA', 'basemeanC', 
                'basemeanG', 'basemeanT']:
        ccdicdat[key] = []
    seenEnd = 0
    savepart = 0
    lastqs = zmwqueue.qsize()
    pbar = False
    while seenEnd < nNoneQueue:
        try:
            nextres = outqueue.get()
        except Queue.Empty:
            time.sleep(0.1)
            pass
        else:
            if nextres is None:
                seenEnd += 1
            else:
                if not pbar:
                    pbar = tqdm(desc='Processing zmws', total=(lastqs - nNoneQueue), smoothing=0.001,
                                mininterval=10, leave=True, position=0)
                pbar.update(1)
                MipdDic[nextres['zmw']] = nextres['full']
                
                for key in ccdicdat.keys():
                    ccdicdat[key].append(nextres[key])
                    
                if len(MipdDic.keys()) >= chunkSize:
                    ccdat = pd.DataFrame(ccdicdat)
                    #reordering zmw info column names
                    ccdat = ccdat[['zmw', 'cclen', 'nsubr', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'basemeanA',
                                   'basemeanC', 'basemeanG', 'basemeanT']] 
                    #need to convert strings to unicode or file can't be loaded in R
                    ccdat.chr = ccdat.chr.astype('unicode') 

                    #save files, one with zmw info, one with a dictionary of onlyT IPD along CCS
                    ccdat.to_pickle(os.path.join(outbase, 'processed', 'full', 'tmp.' + sampcn + '_part' + str(savepart) +
                                                 '_full_zmwinfo.pickle'))
                    with open(os.path.join(outbase, 'processed', 'full', 'tmp.' + sampcn + '_part' + str(savepart) + 
                                           '_full.pickle'), 'wb') as fout:
                        pickle.dump(MipdDic, fout, pickle.HIGHEST_PROTOCOL)

                    savepart += 1
                    MipdDic = {}
                    BingDic = {}
                    ccdicdat = {}
                    for key in ['zmw', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'cclen', 'nsubr', 'basemeanA',
                                'basemeanC', 'basemeanG', 'basemeanT']:
                        ccdicdat[key] = []
    pbar.close()
    ccdat = pd.DataFrame(ccdicdat)
    #reordering zmw info column names
    ccdat = ccdat[['zmw', 'cclen', 'nsubr', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'basemeanA', 'basemeanC',
                   'basemeanG', 'basemeanT']] 
    #need to convert strings to unicode or file can't be loaded in R
    ccdat.chr = ccdat.chr.astype('unicode') 
    
    #save files, one with zmw info, one with a dictionary of onlyT IPD along CCS
    nit = 0
    ccdat.to_pickle(os.path.join(outbase, 'processed', 'full', 'tmp.' + sampcn + '_part' + str(savepart) + '_full_zmwinfo.pickle'))
    with open(os.path.join(outbase, 'processed', 'full', 'tmp.' + sampcn + '_part' + str(savepart) + '_full.pickle'), 'wb') as fout:
        pickle.dump(MipdDic, fout, pickle.HIGHEST_PROTOCOL)

    
    
def main(sarg=None):
    if sarg == None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(sarg)
    sampleRef = pd.read_csv(args.referenceFile, index_col = 'index')
    strre = re.match('(.+)\.(.+)', args.sample)
    if re.match('[\d]+', args.sample) is not None:
        sampleInfo = sampleRef.loc[int(args.sample)]
    elif strre is not None:
        sampCell = strre.groups(0)[0]
        sampName = strre.groups(0)[1]
        sampleInfo = sampleRef.query('cell==@sampCell and sampleName==@sampName').squeeze()
    else:
        raise ValueError("Sample specification must be integer or '[cell].[samplename]'")
    if sampleInfo.shape[0] == 0:
        raise ValueError("No sample in the reference was matched by your sample specification")
    
    # Get output location, set up folders
    if args.outputlocation is not None:
        outBase = args.outputlocation
    else:
        outBase = os.path.join('/avicenna/vramani/analyses/pacbio', sampleInfo.cell)
        
    # make output folders if they don't exist
    if not os.path.exists(os.path.join(outBase, 'processed')):
        os.makedirs(os.path.join(outBase, 'processed'))
    if not os.path.exists(os.path.join(outBase, 'processed','full')):
        os.makedirs(os.path.join(outBase, 'processed','full'))
        
    if args.threads is not None:
        nthreads = args.threads
    else:
        nthreads = 1
    
    sampCN = sampleInfo.cell + '_' + sampleInfo.sampleName
    print(sampCN)
    
    if re.search('linker|array|observed|reference.fasta', sampleInfo.reference) is not None:
        if nthreads > 1:
            print('Multiple threads not implemented for array data')
        #sample is array / monosequence
        
        fullout = extractIPDfullArray(sampleInfo.ccsFile, sampleInfo.alignedSubreadsFile, sampleInfo.reference, args)
        with open(os.path.join(outBase, 'processed', 'full', sampCN + '_full.pickle'), 'wb') as fout:
            pickle.dump(fullout, fout, pickle.HIGHEST_PROTOCOL)
        
    else:
        # sample is genomic
        alttotal = None #8000 # for testing, only do a subset of molecules
        perchunk = 1000
        # Create a list of all zmw hole numbers
        
        validzmw = []
        
        with pb.IndexedBamReader(sampleInfo.ccsFile) as cbam:
            nccs = len(cbam)
            for ic, cc in enumerate(cbam):
                if alttotal and ic >= alttotal:
                    break
                validzmw.append(cc.HoleNumber)
        
        chunksize = perchunk * nthreads
        
        onchunk = 0
        validzmwc = []
        validzmwc.append([])
        for zmw in validzmw:
            validzmwc[onchunk].append(zmw)
            if len(validzmwc[onchunk]) == chunksize:
                onchunk += 1
                validzmwc.append([])
        del validzmw
        nchunks = len(validzmwc)
        
        validzmwQ = mp.JoinableQueue()
        for chunk in range(nchunks):
            for i in range(len(validzmwc[chunk])):
                validzmwQ.put(validzmwc[chunk][i])
            for thr in range(nthreads):
                validzmwQ.put(None)
        del validzmwc
        
        # Start a process to read output and save to disk
        procq = mp.Queue()
        listenp = mp.Process(target=listenerSaver, args=(validzmwQ, procq, nthreads*nchunks, chunksize, outBase, sampCN))
        listenp.daemon = True
        listenp.start()
        
        
        # start processes to extract IPD from each molecule
        
        for chunk in range(nchunks):
            workers = []
            for i in range(nthreads):
                x = mp.Process(target=extractIPDfullgenomic, args=(sampleInfo.ccsFile, sampleInfo.alignedCcsFile,
                                                                 sampleInfo.unalignedSubreadsFile, validzmwQ, procq))
                x.daemon = True
                workers.append(x)
                x.start()
            for thr in workers:
                thr.join()
        
        # wait until all extractIPDfullgenomic processes are finished
        listenp.join()

        # Results are now split between many tmp files. Identify them and join them together
        inffiles = [os.path.basename(x) for x in glob.glob(os.path.join(outBase, 'processed', 'full', 'tmp.' +  sampCN +
                                                                        '_part*_full_zmwinfo.pickle'))]
        onlytfiles = [os.path.basename(x) for x in glob.glob(os.path.join(outBase, 'processed', 'full', 'tmp.' + sampCN +
                                                                          '_part*_full.pickle'))]


        zminfoC = pd.DataFrame()
        zmtipdC = {}

        
        # combine zmw information files
        for inff in inffiles:
            zminfoC = pd.concat([zminfoC, pd.read_pickle(os.path.join(outBase,'processed','full',inff))], sort=False)
        # reset indices for combined zmwinfo dataframe
        zminfoC = zminfoC.sort_values('zmw')
        zminfoC.reset_index(drop=True, inplace=True)
        
        # combine onlyT and binarized files
        for onlytf in onlytfiles:
            with open(os.path.join(outBase,'processed','full',onlytf),'rb') as fopen:
                ottemp = pickle.load(fopen)
            zmtipdC.update(ottemp)
            
            
        # write combined files to pickle
        zminfoC.to_pickle(os.path.join(outBase, 'processed', 'full', sampCN + '_full_zmwinfo.pickle'))
        with open(os.path.join(outBase, 'processed', 'full', sampCN + '_full.pickle'), 'wb') as fout:
            pickle.dump(zmtipdC, fout, pickle.HIGHEST_PROTOCOL)

        # delete the temporary files
        os.system('rm ' + os.path.join(outBase, 'processed', 'full','tmp.' + sampCN + '_part*'))
    
if __name__ == "__main__":
    # Usage: extractIPD.py [referencefile] sampleIdentifier
    main()