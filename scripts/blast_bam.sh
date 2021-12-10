#!/bin/bash

# Run BLAST on BAM file containing long reads

# Overview:
# 1. Convert BAM -> FASTQ -> FASTA
#    - Retain reads with length above {MIN_READ_LEN}
# 2. Run BLAST using {BLASTDB}; convert output to BED format 
#    - Retain matches to database based on {EVALUE} threshold
#    - BED is gzipped and saved as *.raw.bed.gz
# 3. Resolve overlapping hits based on highest bitscore
#    (using dynamic programming) to create pruned BED
#    - Pruned BED is gzipped and saved as *.pruned.bed.gz
#    - Compute various statistics of interest that can be used
#      to filter reads based on matches
# * Formats of BED and statistics files specified in README.txt

if [[ -z $2 ]]; then
    OUTDIR=./
else
    OUTDIR=$2
    mkdir -p $OUTDIR
fi

set -euo pipefail

BAM=$1

# BLAST database prefix
BLASTDB="db/mm_cen_tel"

# E-value threshold
EVALUE=1e-3

BED_OUT_DIR=${OUTDIR}/blast_out
STATS_OUT_DIR=${OUTDIR}

mkdir -p $BED_OUT_DIR
mkdir -p $STATS_OUT_DIR

PREFIX=blast_$(basename $BLASTDB)_$(basename $BAM .bam)

# Run BLAST on reads that are at least MIN_READ_LEN long and retain
# hits that meet the maximum EVALUE threshold; conver to BED
RAW_BED=${BED_OUT_DIR}/${PREFIX}.raw.bed.gz

OUTFMT="6 qseqid qstart qend sseqid evalue sstrand pident sstart send length qlen slen bitscore"

samtools fasta ${BAM}                   \
| blastn -task blastn                   \
    -db ${BLASTDB}                      \
    -evalue ${EVALUE}                   \
    -dust no                            \
    -strand both                        \
    -outfmt "${OUTFMT}"                 \
| awk 'BEGIN {OFS="\t";} {              \
    if ($6=="plus") {                   \
        $6="+"; $8--;                   \
    } else {                            \
        $6="-"; $9--;                   \
        S=$8; $8=$9; $9=S;              \
    }                                   \
    $2--; print $0;                     \
}'                                      \
| gzip                                  \
> ${RAW_BED}

# Resolve overlapping BLAST hits based on bitscore to create a new BED file
# (with same format as above) and compute statistics and save output

PRUNED_BED=${BED_OUT_DIR}/${PREFIX}.pruned.bed.gz
STATS_FILE=${STATS_OUT_DIR}/${PREFIX}.stats.txt.gz

python - $RAW_BED $PRUNED_BED $STATS_FILE <<CODE
import gzip as gz
from itertools import groupby
import sys

def typecast(x):
    """Convert to int or float if possible"""
    try:
        return int(x)
    except ValueError:
        pass
    try:
        return float(x)
    except ValueError:
        return x

_bedhdr = [
    'qseqid', 'qstart', 'qend', 'sseqid', 'evalue', 'sstrand',
    'pident', 'sstart', 'send', 'length', 'qlen', 'slen', 'bitscore'
]

class BlastBedRow:
    def __init__(self, string):
        s = map(typecast, string.split())
        for k, v in zip(_bedhdr, s):
            setattr(self, k, v)
        self.qcovs = None
    
    def __str__(self):
        qcv = ['qcovs'] if (self.qcovs) else []
        fields = map(lambda x: str(self.__dict__.get(x)), _bedhdr + qcv)
        return '\t'.join(fields)
        
def binary_search(hits, i):
    """Binary search on alignment positions."""

    low,high = 0, i-1
    while low <= high:
        mid = (low+high)//2
        if hits[mid].qend <= hits[i].qstart:
            if hits[mid+1].qend <= hits[i].qstart:
                low = mid+1
            else:
                return mid
        else:
            high = mid-1
    return 0
    
def resolve_query_overlaps(hits):
    """Memoized dynamic programming to find best scoring set
    of non-overlapping alignments for a given query sequence"""
    sort = sorted(hits, key=lambda x: x.qend)
    n = len(sort)
    memo = [0]*n
    for i in range(1,n):
        j = binary_search(sort, i)
        memo[i] = max(sort[i].bitscore + memo[j], memo[i-1])
    
    best_hits = []
    i = n-1
    while i > 0:
        j = binary_search(sort, i)
        if sort[i].bitscore + memo[j] > memo[i-1]:
            best_hits.append(sort[i])
            i = j
        else:
            i = i-1
            
    return best_hits[::-1]

def calc_stats(hits):
    """Compute statistics BLAST hits including fraction of query
    covered by matches to subject on the per-subject and across all subject
    levels"""

    for qseqid, qgrp in groupby(hits, lambda x: x.qseqid):
        qgrp = sorted(qgrp, key= lambda x: x.sseqid)
        qlen = qgrp[0].qlen
        qcov = 0
        sids = []
        
        for sseqid, sgrp in groupby(qgrp, lambda x: x.sseqid):
            sgrp = list(sgrp)
            n = len(sgrp)
            scov = sum(map(lambda x: x.qend - x.qstart, sgrp))
            mean_alen = round(sum([x.length for x in sgrp])/n, 4)
            mean_pident = round(sum([x.pident for x in sgrp])/n, 4)
            strand = sorted(set([x.sstrand for x in sgrp]))
            
            # Update query-level statistics
            qcov += scov
            sids.append(sseqid)
            
            yield [
                qseqid, sseqid, n, round(100*scov/qlen, 4), mean_alen,
                mean_pident, strand
            ]
        
        yield [
            qseqid, '*', qlen, len(qgrp), round(100*qcov/qlen, 4),
            len(sids), sorted(set(sids))
        ]

def stats2string(stats):
    stats[-1] = ','.join(stats[-1])
    return '\t'.join(map(str, stats)) + '\n'

with gz.open(sys.argv[1], 'rt') as fh, gz.open(sys.argv[2], 'wt') as oh, gz.open(sys.argv[3], 'wt') as sh:
    bed = map(BlastBedRow, fh)
    for qseqid, hits in groupby(bed, lambda x: x.qseqid):
        pruned = resolve_query_overlaps(hits)
        stats = calc_stats(pruned)
        pruned_str = map(lambda x: str(x)+'\n', pruned)
        oh.writelines(pruned_str)
        sh.writelines(map(stats2string, stats))
CODE
