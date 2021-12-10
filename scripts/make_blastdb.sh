#!/bin/bash

# Set up BLAST databases for mouse repeat analysis

set -euo pipefail

DB_DIR=./db
DB_NAME='mm_cen_tel'

mkdir -p ${DB_DIR}

# Major satellite sequence
cat << EOF > db/major_satellite.fa
>major_satellite M32564.1
CTGAAAGAGGTGGAAAATTTAGAAATGTTCACTCTAGGACGCGAATATGG
CAAGAAAACTAAAAATCATGGAAAATGCGAAACATCCACTTGACGACTTG
AAAATGACGAAATCACTGGAAAACGTGAAAAATGAGAAATGCACACTGTA
GGACCTGGAATATGGCGAGAAAACTGAAAATCAAGTAAAATAAGAAATAT
ACACTTTGGGACGTGAA
EOF

# Minor satellite sequence 
cat << EOF > db/minor_satellite.fa
>minor_satellite X14462.1
GGAAAATGATAAAAACCTACACTGTAGAACATATTAGATGAGTGAGTTAC
ACTGAAAAACACATTCGTTGGAAACCGGCATTGTAGAACAGTGTATATCA
ATGAGTTACAATTAGAAACAT
EOF

# Telomere sequence: pentamerized repeat
cat << EOF > db/telomere.fa
>telomere TTAGGGx5
TTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
EOF

cat ${DB_DIR}/*.fa > ${DB_DIR}/${DB_NAME}.makeblastdb.log

( cat ${DB_DIR}/*.fa                                        \
  | makeblastdb -out ${DB_DIR}/${DB_NAME} -title ${DB_NAME} \
    -dbtype nucl                                            \
) >> ${DB_DIR}/${DB_NAME}.makeblastdb.log 2>&1

# =============================================================
# DFAM major and minor satellite consensus sequences (not used)
# cat << EOF > db/major_satellite.fa
# >major_satellite DF0003028.1 GSAT_MM
# gacctggaatatggcgagaaaactgaaaancgtggaaaatgagaaacgcacactgtagga
# cntggaatatggcgaggaaaactgaaaaacgtggaaaatcgagaaacgcacactgtagga
# cctggaatatggcgagaaaactgaaaancgtggaaaatgagaaacgcacactntaggacc
# tggaatatggcgaaancactgaaaaacgtggaaaatgagaaatgcacactgtaggacctg
# gaatatggcgagaaaactgaaaancgcggaaaatgagaaacgcacactgtaggacntgga
# atatggcgaggaaaactgaaaaacgtggaaaatcgagaaacgcacactgtaggacctgga
# atatggcgagaaaactgaaaancgtggaaaatgagaaacgcacactttaggacctggaat
# atggcgaaaacactgaaaaacgtggaaaatgagaaatgcacactgtaggac
# EOF
# cat << EOF > db/minor_satellite.fa
# >minor_satellite DF0004122.1 SYNREP_MM
# gaacatgttatatgagtgagttacantgaaaaacatanaaatnggaaacgngaattgtag
# aacatcgtatatnaangagttacaatgaaaaacatggaaaatgagaaaaaccacactgta
# gaacattgtgtatgagtgagttacaatgaaaaacatananatnggaaacgaganttgtag
# aacatcgtatataaaggagttacaatgaaaaacatggaaaatgagaaaaaccacactgta
# gaacattgtgtatgagtgagttaggg
# EOF