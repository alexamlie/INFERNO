#!/bin/bash

## compare_FANTOM5_roadmap.sh
## alex amlie-wolf 03/14/18
## a script to compare FANTOM5 and Roadmap enhancer states within tissue categories

if [ $# == 5 ]; then
    OUTDIR=$1
    F5_CLASSES=$2
    F5_DIR=$3
    HMM_CLASSES=$4
    HMM_DIR=$5

    mkdir -p ${OUTDIR}/FANTOM5_merged/ ${OUTDIR}/roadmap_merged/ ${OUTDIR}/class_comparisons/
    ## clear any previous files
    rm -rf ${OUTDIR}/FANTOM5_merged/*
    rm -rf ${OUTDIR}/roadmap_merged/*

    echo "Parsing FANTOM5 data"
    while IFS=$'\t' read -r -a LINE_DATA; do
    	THIS_CLASS=`echo ${LINE_DATA[0]} | sed -e 's/ /_/g'`
    	THIS_FILE=${LINE_DATA[1]}
    	echo "${THIS_CLASS}, ${THIS_FILE}"
	
    	cat ${F5_DIR}/${THIS_FILE} >> ${OUTDIR}/FANTOM5_merged/${THIS_CLASS}.all_intervals.bed
    done < <(tail -n +2 ${F5_CLASSES})

    ## similar idea for Roadmap data
    echo "Parsing Roadmap data"
    while IFS=$'\t' read -r -a LINE_DATA; do
    	THIS_CLASS=`echo ${LINE_DATA[0]} | sed -e 's/ /_/g'`
    	THIS_EID=${LINE_DATA[2]}
    	echo "${THIS_CLASS}, ${THIS_EID}"
	
    	awk -F$'\t' '$4=="6_EnhG" || $4=="7_Enh" || $4=="12_EnhBiv"' \
    	    ${HMM_DIR}/${THIS_EID}_15_coreMarks_mnemonics.bed \
    	    >> ${OUTDIR}/roadmap_merged/${THIS_CLASS}.all_intervals.bed
    done < <(tail -n +2 ${HMM_CLASSES})

    module load bedtools2
    ## now merge the files in each class
    echo "Merging class files"
    for BEDF in ${OUTDIR}/FANTOM5_merged/*all_intervals.bed; do
    	THIS_DIR=`dirname ${BEDF}`
    	THIS_CLASS=`basename ${BEDF%.all_intervals.bed}`
    	echo "Merging ${THIS_CLASS} data in FANTOM5"
    	sort -k1,1 -k2,2n ${BEDF} | bedtools merge -i stdin > ${THIS_DIR}/${THIS_CLASS}.merged_intervals.bed
    done

    for BEDF in ${OUTDIR}/roadmap_merged/*all_intervals.bed; do
    	THIS_DIR=`dirname ${BEDF}`
    	THIS_CLASS=`basename ${BEDF%.all_intervals.bed}`
    	echo "Merging ${THIS_CLASS} data in roadmap"
    	sort -k1,1 -k2,2n ${BEDF} | bedtools merge -i stdin > ${THIS_DIR}/${THIS_CLASS}.merged_intervals.bed
    done
    
    ## and finally for each class, compare the two data sources to generate files to read into R
    ## get all the unique classes
    ## this formulation will let me figure out the unique numbers easily
    echo -e "class\tFANTOM5_total\troadmap_total\tshared_intervals" > \
	${OUTDIR}/FANTOM5_roadmap_overlap_summary.txt
    
    for CLASS in `find ${OUTDIR}/*_merged/ -name '*merged_intervals.bed' -printf '%f\n' | sed -e 's/.merged_intervals.bed//g' | sort -u`; do
	echo "Comparing results for $CLASS"

	F5_FILE=${OUTDIR}/FANTOM5_merged/${CLASS}.merged_intervals.bed
	HMM_FILE=${OUTDIR}/roadmap_merged/${CLASS}.merged_intervals.bed
	
	## if they both have files for this class, intersect them and summarize the results
	if [[ -e "${F5_FILE}" && -e "${HMM_FILE}" ]]; then
	    echo "Computing intersection"
	    ## just use the normal intersect call, to get the shared intervals only
	    bedtools intersect -a ${F5_FILE} -b ${HMM_FILE} > ${OUTDIR}/class_comparisons/${CLASS}_shared_intervals.bed
	    F5_TOT=`wc -l ${F5_FILE} | awk '{print $1}'`
	    HMM_TOT=`wc -l ${HMM_FILE} | awk '{print $1}'`
	    SHARED_COUNT=`wc -l ${OUTDIR}/class_comparisons/${CLASS}_shared_intervals.bed | awk '{print $1}'`

	    echo -e "${CLASS}\t${F5_TOT}\t${HMM_TOT}\t${SHARED_COUNT}" >> ${OUTDIR}/FANTOM5_roadmap_overlap_summary.txt
	elif [ -e "${F5_FILE}" ]; then
	    echo "Only FANTOM5 found"
	    F5_TOT=`wc -l ${F5_FILE} | awk '{print $1}'`
	    echo -e "${CLASS}\t${F5_TOT}\t0\t0" >> ${OUTDIR}/FANTOM5_roadmap_overlap_summary.txt
	else
	    echo "Only Roadmap found"
	    HMM_TOT=`wc -l ${HMM_FILE} | awk '{print $1}'`
	    echo -e "${CLASS}\t0\t${HMM_TOT}\t0" >> ${OUTDIR}/FANTOM5_roadmap_overlap_summary.txt
	fi	    
    done
    
fi
