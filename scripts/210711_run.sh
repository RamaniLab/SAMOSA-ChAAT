#!/bin/bash

# Wrapper script that calls make_blastdb.sh

TOP_DIR='/avicenna/vramani/analyses/pacbio'

declare -a RUN_LIST=(
    '210427_NA_SAMv2mESCs_1.0prep'
    '210427_NA_SAMv2mESCs_2.0prep'
    '210421_NA_SAMv2_mESCs'
    'pbrun11_mESCs_SNF2h'
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