#!/bin/sh

## calculate_tss_distance.sh
## alex amlie-wolf 04/15/2016
## a script that calculates the distance of vcf files to annotations in a bed file
## it is meant to be used to calculate the distance of SNPs to their nearest TSS
## takes in the data directory, the prefix of the vcf files to be analyzed, the bed file,
## the path to bedtools (or a direction to load it with modules), and the output directory

if [ $# == 5 ]; then
    DATADIR=$1
    FPREFIX=$2
    BEDF=$3
    BEDTOOLS_BIN_DIR=$4
    OUTDIR=$5

    mkdir -p ${OUTDIR}

    ## just use the bedtools bin path if we have it
    ${BEDTOOLS_BIN_DIR}/bedtools --version > ${OUTDIR}/BEDTOOLS_VERSION.TXT
    
    for F in ${DATADIR}/${FPREFIX}*; do
	FNAME=`basename $F`
	OUTPREFIX=${FNAME%.vcf*}
	echo "Analyzing file ${FNAME}"
	## TODO: check if it's gzipped
	## we want to find the closest gene, but we need it to be in (sorted) bed format
	zcat ${F} | awk -v OFS=$'\t' '{if (substr($1, 1, 1)!="#") print "chr"$1, $2, $2+1, $3}' \
	    | sort -k1,1 -k2,2n | \
	    ${BEDTOOLS_BIN_DIR}/bedtools closest -d -t first -a stdin -b ${BEDF} \
	    > ${OUTDIR}/${OUTPREFIX}_closest_tss.txt	
    done
elif [ $# == 4 ]; then
    ## otherwise, we don't have a bedtools bin, so try to load it as module
    DATADIR=$1
    FPREFIX=$2
    BEDF=$3
    OUTDIR=$4

    mkdir -p ${OUTDIR}    

    ## if we don't have bedtools, then try to load it as a module
    ## first try just normal bedtools
    if ! type "bedtools" > /dev/null 2>&1; then
	## try loading bedtools or bedtools 2
	if type "module" > /dev/null 2>&1; then
	    module load bedtools > /dev/null 2>&1
	    ## if this still didn't work
	    if ! type "bedtools" > /dev/null 2>&1; then
		module load bedtools2 > /dev/null 2>&1
		if ! type "bedtools" > /dev/null 2>&1; then
		    echo "No bedtools found! Please insall bedtools and put it in your \$PATH, install it as a module, or provide a link to a bin/ directory for bedtools."
		    exit 1
		fi
	    fi
	fi
    fi
    ## if we reach this stage, then some bedtools module loading should have worked
    bedtools --version > ${OUTDIR}/BEDTOOLS_VERSION.TXT
    
    for F in ${DATADIR}/${FPREFIX}*; do
	FNAME=`basename $F`
	OUTPREFIX=${FNAME%.vcf*}
	echo "Analyzing file ${FNAME}"
	## TODO: check if it's gzipped
	## we want to find the closest gene, but we need it to be in (sorted) bed format
	zcat ${F} | awk -v OFS=$'\t' '{if (substr($1, 1, 1)!="#") print "chr"$1, $2, $2+1, $3}' \
	    | sort -k1,1 -k2,2n | bedtools closest -d -t first -a stdin -b ${BEDF} \
	    > ${OUTDIR}/${OUTPREFIX}_closest_tss.txt	
    done
fi
    
