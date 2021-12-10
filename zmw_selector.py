import os,sys,re
import scipy
import numpy as np
import pysam 
from collections import Counter
from intervaltree import Interval, IntervalTree

def readIterator(filenames,chrom, start, end):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom, start, end):
        yield read
      input_file.close()

def define_large_bins(chromsizes, resolutions):
    bins = {}
    valid_chroms = {}
    lines = chromsizes.readlines()
    for resolution in resolutions:
        bins[resolution] = {}
    for resolution in resolutions:
        hindex = 0
        for line in lines:
            chromname, length = line.split()
            valid_chroms[chromname] = True
            for i in range(0,int(length),resolution):
                bins[resolution][(chromname, i)] = hindex
                hindex += 1
    return bins, valid_chroms


def calculatePerBase(filenames, tss, valid_chroms, window):
    '''the goal here is to take a list of flat-file formatted features & 
    extract all ZMWs (from aligned CCS) where a portion of the read falls within n bases of the feature. '''
    chc = False
    zmws = []
    for line in tss:
        split = line.split()
        #check if there is an optional label (for certain flat files)
        if len(split) > 3:
            label = split[3]
        chrom = split[0]
        if chrom not in valid_chroms: continue
        t_start, strand = int(split[1]), split[2]
        start = t_start - window
        if start <= 0: continue
        end = t_start + window
        end_mat = np.zeros(5001)
        for read in readIterator(filenames, chrom, start, end):
            if read.is_reverse:
                rstrand = '-'
            else:
                rstrand = '+'
            rname = read.qname
            rstart = read.reference_start
            rend = read.reference_end
            zmw = rname.split('/')[1]
            if len(split) > 3:
                zmws.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (zmw, chrom, rstart, rend, t_start,label,rstrand, strand))
            else:
                zmws.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (zmw, chrom, rstart, rend, t_start,rstrand, strand))
            
    return zmws

def main():
    valid = {}
    tfile = open(sys.argv[1])
    valid_chroms = open(sys.argv[2])
    for line in valid_chroms:
        split = line.split()
        valid[split[0]] = True
    tss = tfile.readlines()
    filename = sys.argv[3]
    window = int(sys.argv[4])
    filenames = [filename]
    zmws = calculatePerBase(filenames, tss, valid, window)
    for zmw in zmws:
        print(zmw)

if __name__ == "__main__":
    main()
