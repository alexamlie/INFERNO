#!/bin/bash

## calculate_maf.sh
## alex amlie-wolf 04/19/2016
## a script that calculates minor allele frequencies for vcf files in a data directory
## takes in the data directory, the prefix of the vcf files to be analyzed, and the output directory

if [ $# == 3 ]; then
    DATADIR=$1
    FPREFIX=$2
    OUTDIR=$3

    mkdir -p ${OUTDIR}
    
    for F in ${DATADIR}/${FPREFIX}*; do
	FNAME=`basename $F`
	OUTPREFIX=${FNAME%.vcf*}
	## now run plink, but specifically to calculate MAF
	echo "Analyzing file ${FNAME}"		
	plink --vcf ${F} --freq --out ${OUTDIR}/${OUTPREFIX}
	## now convert the output file into tab separated format:
	sed -r 's/^\s+|\s+$//g' ${OUTDIR}/${OUTPREFIX}.frq | sed -r 's/\s+/\t/g' > ${OUTDIR}/${OUTPREFIX}.frq.parsed
	## overwrite the original PLINK output file
	mv ${OUTDIR}/${OUTPREFIX}.frq.parsed ${OUTDIR}/${OUTPREFIX}.frq

	## now, annotate the frequency file with the positions and write it out
	printf "chr\tpos\trsID\tA1\tA2\tMAF\tNCHROBs\n" > ${OUTDIR}/${OUTPREFIX}.frq.full
	paste <(zcat ${F} | grep -v "#" | cut -f2) <(tail -n +2 ${OUTDIR}/${OUTPREFIX}.frq) | \
	    awk -v OFS=$'\t' '{print "chr"$2, $1, $3, $4, $5, $6, $7}' >> ${OUTDIR}/${OUTPREFIX}.frq.full
	mv ${OUTDIR}/${OUTPREFIX}.frq.full ${OUTDIR}/${OUTPREFIX}.frq
    done
fi
