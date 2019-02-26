#!/bin/sh

## enhancer_bootstrap_bsub_wrapper.sh
## alex amlie-wolf 05/04/17
## a wrapper script to submit the enhancer bootstrapping script to bsub

if [ $# == 12 ]; then
    BS_SCRIPT=$1
    NUM_SAMPLES=$2
    MAF_BIN_SIZE=$3
    DIST_ROUND=$4
    DIST_THRESHOLD=$5
    LD_PARTNER_THRESHOLD=$6
    BG_SNP_INFOF=$7
    LD_SETS_DIR=$8
    REF_SUMMARY_DIR=$9
    LD_THRESH=${10}
    LD_AREA=${11}
    PARAMF=${12}
    
#    module load R/3.2.3

    Rscript "${BS_SCRIPT}" ${NUM_SAMPLES} ${MAF_BIN_SIZE} ${DIST_ROUND} \
    	${DIST_THRESHOLD} ${LD_PARTNER_THRESHOLD} ${BG_SNP_INFOF} ${LD_SETS_DIR} \
    	${REF_SUMMARY_DIR} ${LD_THRESH} ${LD_AREA} ${PARAMF}
elif [ $# == 13 ]; then
    ## running in shell mode, use tee
    BS_SCRIPT=$1
    NUM_SAMPLES=$2
    MAF_BIN_SIZE=$3
    DIST_ROUND=$4
    DIST_THRESHOLD=$5
    LD_PARTNER_THRESHOLD=$6
    BG_SNP_INFOF=$7
    LD_SETS_DIR=$8
    REF_SUMMARY_DIR=$9
    LD_THRESH=${10}
    LD_AREA=${11}
    PARAMF=${12}
    LOGFILE=${13}
    
#    module load R/3.2.3

    Rscript "${BS_SCRIPT}" ${NUM_SAMPLES} ${MAF_BIN_SIZE} ${DIST_ROUND} \
    	${DIST_THRESHOLD} ${LD_PARTNER_THRESHOLD} ${BG_SNP_INFOF} ${LD_SETS_DIR} \
    	${REF_SUMMARY_DIR} ${LD_THRESH} ${LD_AREA} ${PARAMF} 2>&1 | tee ${LOGFILE}        
else
    echo "Usage: $0 <bootstrap script> <number of background samples> <MAF bin size> <distance rounding> <LD distance threshold> <LD partner number threshold> <precomputed background SNP info file> <directory with precomputed LD sets> <reference summary directory> <LD threshold> <LD distance threshold> <INFERNO parameter file>"
fi
