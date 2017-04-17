#!/bin/bash

## calculate_pairwise_ld.sh
## alex amlie-wolf 04/19/2016
## a script that calculates pairwise LD information for files in some input directory
## takes in the data directory, the prefix of the vcf files to be analyzed, and the output directory

if [ $# == 3 ]; then
    DATADIR=$1
    FPREFIX=$2
    OUTDIR=$3

    mkdir -p ${OUTDIR}
    
    for F in ${DATADIR}/${FPREFIX}*; do
	FNAME=`basename $F`
	OUTPREFIX=${FNAME%.vcf*}
	## now run PLINK
	echo "Analyzing file ${FNAME}"
	## find all SNPs within a megabase
	plink --vcf ${F} --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.1 --out ${OUTDIR}/${OUTPREFIX}
	## now convert the output file into tab separated format:
	sed -r 's/^\s+|\s+$//g' ${OUTDIR}/${OUTPREFIX}.ld | sed -r 's/\s+/\t/g' > ${OUTDIR}/${OUTPREFIX}.ld.parsed
	## overwrite the original PLINK output file
	mv ${OUTDIR}/${OUTPREFIX}.ld.parsed ${OUTDIR}/${OUTPREFIX}.ld
	## compress this result to save space!
	gzip ${OUTDIR}/${OUTPREFIX}.ld
    done
fi
