#!/bin/sh

## run_ld_score_regression.sh
## alex amlie-wolf 06/18/18
## a script to run LD score regression

if [ $# == 16 ]; then
    LDSC_CODE_DIR=${1}
    OUTDIR=${2}
    SUMMARY_STATF=${3}
    MUNGE_SNPLIST=${4}
    RSID_COL=${5}
    PVAL_COL=${6}
    ALLELE1_COL=${7}
    ALLELE2_COL=${8}
    MAF_COL=${9}
    BETA_COL=${10}    
    HAS_HEADER=${11}
    LDSC_BASELINE_DIR=${12}
    LDSC_WEIGHTS_DIR=${13}
    LDSC_FRQ_DIR=${14}
    OUTPREFIX=${15}
    N=${16}
    
    mkdir -p ${OUTDIR}
    
    cd ${LDSC_CODE_DIR}
    ## comment this out if not on a module-based system
    module load anaconda/3
    source activate ldsc    
    
    if [ ${HAS_HEADER} == "True" ] ; then
	## need to get the header names
	RSID_NAME=`head -1 ${SUMMARY_STATF} | cut -f${RSID_COL}`
	PVAL_NAME=`head -1 ${SUMMARY_STATF} | cut -f${PVAL_COL}`
	ALLELE1_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE1_COL}`
	ALLELE2_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE2_COL}`
	MAF_NAME=`head -1 ${SUMMARY_STATF} | cut -f${MAF_COL}`	
	BETA_NAME=`head -1 ${SUMMARY_STATF} | cut -f${BETA_COL}`

	python ./munge_sumstats.py --sumstats ${SUMMARY_STATF} --merge-alleles ${MUNGE_SNPLIST} \
	    --out ${OUTDIR}/${OUTPREFIX}.munged_data --snp "${RSID_NAME}" --a1 "${ALLELE1_NAME}" \
	    --a2 "${ALLELE2_NAME}" --p "${PVAL_NAME}" --frq "${MAF_NAME}" \
	    --signed-sumstats "${BETA_NAME},0" --N ${N}

	python ./ldsc.py --h2 ${OUTDIR}/${OUTPREFIX}.munged_data.sumstats.gz --ref-ld-chr \
	    ${LDSC_BASELINE_DIR}/baseline. --w-ld-chr ${LDSC_WEIGHTS_DIR}/weights. \
	    --overlap-annot --frqfile-chr ${LDSC_FRQ_DIR}/1000G.mac5eur. \
	    --out ${OUTDIR}/${OUTPREFIX}.LDSC
    else
	NUM_COLS=`head -1 ${SUMMARY_STATF} | awk -F$'\t' '{print NF}'`
	## reformat summary stat file to include header
	awk -F$'\t' -v RSID=${RSID_COL} -v PVAL=${PVAL_COL} -v A1=${ALLELE1_COL} \
	    -v A2=${ALLELE2_COL} -v MAF=${MAF_COL} -v BETA=${BETA_COL} -v NUM_COLS=${NUM_COLS} \
	    'BEGIN{OFS=FS; for(i=1; i<=NUM_COLS; i++) {
                   if(i==RSID) {printf "rsID"} else if(i==PVAL) {printf "pval"} 
                   else if(i==A1) {printf "allele1"} else if(i==A2) {printf "allele2"}
                   else if(i==MAF) {printf "MAF"} else if(i==BETA) {printf "beta"}
                   else {printf i}; if(i!=NUM_COLS) {printf "\t"} else {printf "\n"}}} {print $0}' \
	    ${SUMMARY_STATF} > ${OUTDIR}/${OUTPREFIX}.summary_stats.header_added.txt
	
	python ./munge_sumstats.py --sumstats ${OUTDIR}/${OUTPREFIX}.summary_stats.header_added.txt \
	    --merge-alleles ${MUNGE_SNPLIST} --N ${N} \
	    --out ${OUTDIR}/${OUTPREFIX}.munged_data --snp "rsID" --a1 "allele1" \
	    --a2 "allele2" --p "pval" --frq "MAF" --signed-sumstats "beta,0"

	python ./ldsc.py --h2 ${OUTDIR}/${OUTPREFIX}.munged_data.sumstats.gz --ref-ld-chr \
	    ${LDSC_BASELINE_DIR}/baseline. --w-ld-chr ${LDSC_WEIGHTS_DIR}/weights. \
	    --overlap-annot --frqfile-chr ${LDSC_FRQ_DIR}/1000G.mac5eur. \
	    --out ${OUTDIR}/${OUTPREFIX}.LDSC	
    fi

    source deactivate ldsc
elif if [ $# == 17 ]; then
    LDSC_CODE_DIR=${1}
    OUTDIR=${2}
    SUMMARY_STATF=${3}
    MUNGE_SNPLIST=${4}
    RSID_COL=${5}
    PVAL_COL=${6}
    ALLELE1_COL=${7}
    ALLELE2_COL=${8}
    MAF_COL=${9}
    BETA_COL=${10}    
    HAS_HEADER=${11}
    LDSC_BASELINE_DIR=${12}
    LDSC_WEIGHTS_DIR=${13}
    LDSC_FRQ_DIR=${14}
    OUTPREFIX=${15}
    N=${16}
    LOGFILE=${17}
    
    mkdir -p ${OUTDIR}
    
    cd ${LDSC_CODE_DIR}
    ## comment this out if not on a module-based system
    module load anaconda/3
    source activate ldsc    
    
    if [ ${HAS_HEADER} == "True" ] ; then
	## need to get the header names
	RSID_NAME=`head -1 ${SUMMARY_STATF} | cut -f${RSID_COL}`
	PVAL_NAME=`head -1 ${SUMMARY_STATF} | cut -f${PVAL_COL}`
	ALLELE1_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE1_COL}`
	ALLELE2_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE2_COL}`
	MAF_NAME=`head -1 ${SUMMARY_STATF} | cut -f${MAF_COL}`	
	BETA_NAME=`head -1 ${SUMMARY_STATF} | cut -f${BETA_COL}`

	python ./munge_sumstats.py --sumstats ${SUMMARY_STATF} --merge-alleles ${MUNGE_SNPLIST} \
	    --out ${OUTDIR}/${OUTPREFIX}.munged_data --snp "${RSID_NAME}" --a1 "${ALLELE1_NAME}" \
	    --a2 "${ALLELE2_NAME}" --p "${PVAL_NAME}" --frq "${MAF_NAME}" \
	    --signed-sumstats "${BETA_NAME},0" --N ${N} 2>&1 | tee ${LOGFILE}

	python ./ldsc.py --h2 ${OUTDIR}/${OUTPREFIX}.munged_data.sumstats.gz --ref-ld-chr \
	    ${LDSC_BASELINE_DIR}/baseline. --w-ld-chr ${LDSC_WEIGHTS_DIR}/weights. \
	    --overlap-annot --frqfile-chr ${LDSC_FRQ_DIR}/1000G.mac5eur. \
	    --out ${OUTDIR}/${OUTPREFIX}.LDSC 2>&1 | tee -a ${LOGFILE}
    else
	NUM_COLS=`head -1 ${SUMMARY_STATF} | awk -F$'\t' '{print NF}'`
	## reformat summary stat file to include header
	awk -F$'\t' -v RSID=${RSID_COL} -v PVAL=${PVAL_COL} -v A1=${ALLELE1_COL} \
	    -v A2=${ALLELE2_COL} -v MAF=${MAF_COL} -v BETA=${BETA_COL} -v NUM_COLS=${NUM_COLS} \
	    'BEGIN{OFS=FS; for(i=1; i<=NUM_COLS; i++) {
                   if(i==RSID) {printf "rsID"} else if(i==PVAL) {printf "pval"} 
                   else if(i==A1) {printf "allele1"} else if(i==A2) {printf "allele2"}
                   else if(i==MAF) {printf "MAF"} else if(i==BETA) {printf "beta"}
                   else {printf i}; if(i!=NUM_COLS) {printf "\t"} else {printf "\n"}}} {print $0}' \
	    ${SUMMARY_STATF} > ${OUTDIR}/${OUTPREFIX}.summary_stats.header_added.txt
	
	python ./munge_sumstats.py --sumstats ${OUTDIR}/${OUTPREFIX}.summary_stats.header_added.txt \
	    --merge-alleles ${MUNGE_SNPLIST} --N ${N} \
	    --out ${OUTDIR}/${OUTPREFIX}.munged_data --snp "rsID" --a1 "allele1" \
	    --a2 "allele2" --p "pval" --frq "MAF" --signed-sumstats "beta,0" 2>&1 | tee ${LOGFILE}

	python ./ldsc.py --h2 ${OUTDIR}/${OUTPREFIX}.munged_data.sumstats.gz --ref-ld-chr \
	    ${LDSC_BASELINE_DIR}/baseline. --w-ld-chr ${LDSC_WEIGHTS_DIR}/weights. \
	    --overlap-annot --frqfile-chr ${LDSC_FRQ_DIR}/1000G.mac5eur. \
	    --out ${OUTDIR}/${OUTPREFIX}.LDSC 2>&1 | tee -a ${LOGFILE}	
    fi

    source deactivate ldsc
else
    echo "Incorrect number of arguments"
fi
