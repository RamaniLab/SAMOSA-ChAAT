#!/bin/bash

# Wrapper script that calls make_blastdb.sh

TOP_DIR='/avicenna/vramani/analyses/pacbio'

declare -a RUN_LIST=(
    'SAMv2_SNF2hWTAB_mESCs_+m_1.0prep'
    'SAMv2_SNF2hKO_mESCs_+m1_1.0prep2'
    'SAMv2_SNF2hKO_mESCs_+m1_1.0prep'
    'SAMv2_E14_mESCs_+m_preptest'
    '210830_NA_SNF2H-pool_c3'
    '210830_NA_SNF2H-pool_c2'
    '210830_NA_SNF2H-pool_c1'
)

LOG_DIR=logs_`date +%Y%m%d`
mkdir -p $LOG_DIR

for RUN in ${RUN_LIST[@]}
do
    echo "$RUN"
    mkdir -p $RUN
    for BAM in ${TOP_DIR}/${RUN}/aligned/*.bam
    do
        PFX=`basename $BAM .bam`
        echo "     $PFX"
        nohup ./blast_bam.sh $BAM $RUN > ${LOG_DIR}/nohup_${PFX}.log 2>&1 &
    done
done
