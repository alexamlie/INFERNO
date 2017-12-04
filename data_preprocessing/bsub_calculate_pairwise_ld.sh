#!/bin/bash

## bsub_calculate_pairwise_ld.sh
## alex amlie-wolf 04/19/2016
## a script that wraps the LD calculation script in bsub calls

if [ $# == 3 ]; then
    DATADIR=$1
    FPREFIX=$2
    OUTDIR=$3

    DATESTRING=`date +%m_%d_%y_%R:%S`
    
    mkdir -p ${OUTDIR}/logs/${DATESTRING}/
    
    for F in ${DATADIR}/${FPREFIX}*; do
	FNAME=`basename $F`
	OUTPREFIX=${FNAME%.vcf*}
	## submit this to the pairwise LD calculation (assumed to be in this directory) the
	## trick is that the full file name is used as the prefix, so that only that file gets
	## analyzed. also, need to ask for a lot of memory (40 Gb!)
	bsub -J ${OUTPREFIX}_ld_calc -e ${OUTDIR}/logs/${DATESTRING}/${OUTPREFIX}.err \
	    -o ${OUTDIR}/logs/${DATESTRING}/${OUTPREFIX}.out -M 45000 \
	    -q wang_normal \
	    ./calculate_pairwise_ld.sh ${DATADIR} ${FNAME} ${OUTDIR}
    done
else
    ## run files directly, assume that we get a bunch of prefixes
    DATADIR=$1
    OUTDIR=$2
    shift 2
    
    DATESTRING=`date +%m_%d_%y_%R:%S`
    
    mkdir -p ${OUTDIR}/logs/${DATESTRING}/

    for FPREFIX in "$@"; do
	for F in ${DATADIR}/${FPREFIX}*; do 
    	    FNAME=`basename $F`
    	    OUTPREFIX=${FNAME%.vcf*}
    	    ./calculate_pairwise_ld.sh ${DATADIR} ${FNAME} ${OUTDIR} 2>&1 | tee ${OUTDIR}/logs/${DATESTRING}/${OUTPREFIX}_full_log.txt 
	done
    done
fi

