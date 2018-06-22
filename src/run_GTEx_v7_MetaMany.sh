#!/bin/sh

## run_GTEx_v7_MetaMany.sh
## alex amlie-wolf 06/13/18
## a script as part of INFERNO, to run MetaXcan on all the GTEx databases

if [ $# == 14 ]; then
    METAXCAN_DIR=${1}
    OUTDIR=${2}
    SUMMARY_STATF=${3}
    GTEX_V7_DBDIR=${4}
    RSID_COL=${5}
    POS_COL=${6}
    PVAL_COL=${7}
    CHR_COL=${8}
    ALLELE1_COL=${9}
    ALLELE2_COL=${10}
    MAF_COL=${11}
    EFFECT_COL=${12}
    HAS_HEADER=${13}
    USE_BETA=${14}
    
    mkdir -p ${OUTDIR}
    
    if [ ${HAS_HEADER} == "True" ] ; then
	## need to get the header names
	RSID_NAME=`head -1 ${SUMMARY_STATF} | cut -f${RSID_COL}`
	POS_NAME=`head -1 ${SUMMARY_STATF} | cut -f${POS_COL}`
	PVAL_NAME=`head -1 ${SUMMARY_STATF} | cut -f${PVAL_COL}`
	CHR_NAME=`head -1 ${SUMMARY_STATF} | cut -f${CHR_COL}`
	ALLELE1_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE1_COL}`
	ALLELE2_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE2_COL}`
	MAF_NAME=`head -1 ${SUMMARY_STATF} | cut -f${MAF_COL}`	
	EFFECT_NAME=`head -1 ${SUMMARY_STATF} | cut -f${EFFECT_COL}`

	if [ ${USE_BETA} == "True" ]; then
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${SUMMARY_STATF} --snp_column ${RSID_NAME} \
		--pvalue_column ${PVAL_NAME} --effect_allele_column ${ALLELE1_NAME} \
		--non_effect_allele_column ${ALLELE2_NAME} --chromosome_column ${CHR_NAME} \
		--position_column ${POS_NAME} --beta_column ${EFFECT_NAME} --freq_column ${MAF_NAME} \
		--separator $'\t' --covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db
	else
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${SUMMARY_STATF} --snp_column ${RSID_NAME} \
		--pvalue_column ${PVAL_NAME} --effect_allele_column ${ALLELE1_NAME} \
		--non_effect_allele_column ${ALLELE2_NAME} --chromosome_column ${CHR_NAME} \
		--position_column ${POS_NAME} --or_column ${EFFECT_NAME} --freq_column ${MAF_NAME} \
		--separator $'\t' --covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db
	fi
    else 
	NUM_COLS=`head -1 ${SUMMARY_STATF} | awk -F$'\t' '{print NF}'`
	## reformat summary stat file to include header
	awk -F$'\t' -v RSID=${RSID_COL} -v POS=${POS_COL} -v PVAL=${PVAL_COL} \
	    -v CHR=${CHR_COL} -v A1=${ALLELE1_COL} -v A2=${ALLELE2_COL} -v MAF=${MAF_COL} \
	    -v EFFECT=${EFFECT_COL} -v NUM_COLS=${NUM_COLS} \
	    'BEGIN{OFS=FS; for(i=1; i<=NUM_COLS; i++) {
                   if(i==RSID) {printf "rsID"} else if(i==POS) {printf "pos"}
                   else if(i==PVAL) {printf "pval"} else if(i==CHR) {printf "chr"}
                   else if(i==A1) {printf "allele1"} else if(i==A2) {printf "allele2"}
                   else if(i==MAF) {printf "MAF"} else if(i==EFFECT) {printf "effect"}
                   else {printf i}; if(i!=NUM_COLS) {printf "\t"} else {printf "\n"}}} {print $0}' \
	    ${SUMMARY_STATF} > ${OUTDIR}/summary_stats.header_added.txt

	if [ ${USE_BETA} == "True" ]; then
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${OUTDIR}/summary_stats.header_added.txt \
		--snp_column "rsID" --pvalue_column "pval" \
		--effect_allele_column "allele1" --non_effect_allele_column "allele2" \
		--chromosome_column "chr" --position_column "pos" --beta_column \
		"effect" --freq_column "MAF" --separator $'\t' \
		--covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db
	else
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${OUTDIR}/summary_stats.header_added.txt \
		--snp_column "rsID" --pvalue_column "pval" \
		--effect_allele_column "allele1" --non_effect_allele_column "allele2" \
		--chromosome_column "chr" --position_column "pos" --or_column \
		"effect" --freq_column "MAF" --separator $'\t' \
		--covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db
	fi
    fi

    ## after it's done, run the analysis R script
    Rscript ./src/analyze_GTEx_v7_metaXcan.R ${OUTDIR}
elif [ $# == 15 ]; then
    METAXCAN_DIR=${1}
    OUTDIR=${2}
    SUMMARY_STATF=${3}
    GTEX_V7_DBDIR=${4}
    RSID_COL=${5}
    POS_COL=${6}
    PVAL_COL=${7}
    CHR_COL=${8}
    ALLELE1_COL=${9}
    ALLELE2_COL=${10}
    MAF_COL=${11}
    EFFECT_COL=${12}
    HAS_HEADER=${13}
    USE_BETA=${14}
    LOGFILE=${15}
    
    mkdir -p ${OUTDIR}
    
    if [ ${HAS_HEADER} == "True" ] ; then
	## need to get the header names
	RSID_NAME=`head -1 ${SUMMARY_STATF} | cut -f${RSID_COL}`
	POS_NAME=`head -1 ${SUMMARY_STATF} | cut -f${POS_COL}`
	PVAL_NAME=`head -1 ${SUMMARY_STATF} | cut -f${PVAL_COL}`
	CHR_NAME=`head -1 ${SUMMARY_STATF} | cut -f${CHR_COL}`
	ALLELE1_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE1_COL}`
	ALLELE2_NAME=`head -1 ${SUMMARY_STATF} | cut -f${ALLELE2_COL}`
	MAF_NAME=`head -1 ${SUMMARY_STATF} | cut -f${MAF_COL}`	
	EFFECT_NAME=`head -1 ${SUMMARY_STATF} | cut -f${EFFECT_COL}`

	if [ ${USE_BETA} == "True" ]; then
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${SUMMARY_STATF} --snp_column ${RSID_NAME} \
		--pvalue_column ${PVAL_NAME} --effect_allele_column ${ALLELE1_NAME} \
		--non_effect_allele_column ${ALLELE2_NAME} --chromosome_column ${CHR_NAME} \
		--position_column ${POS_NAME} --beta_column ${EFFECT_NAME} --freq_column ${MAF_NAME} \
		--separator $'\t' --covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db 2>&1 | tee ${LOGFILE}
	else
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${SUMMARY_STATF} --snp_column ${RSID_NAME} \
		--pvalue_column ${PVAL_NAME} --effect_allele_column ${ALLELE1_NAME} \
		--non_effect_allele_column ${ALLELE2_NAME} --chromosome_column ${CHR_NAME} \
		--position_column ${POS_NAME} --or_column ${EFFECT_NAME} --freq_column ${MAF_NAME} \
		--separator $'\t' --covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db 2>&1 | tee ${LOGFILE}
	fi
    else 
	NUM_COLS=`head -1 ${SUMMARY_STATF} | awk -F$'\t' '{print NF}'`
	## reformat summary stat file to include header
	awk -F$'\t' -v RSID=${RSID_COL} -v POS=${POS_COL} -v PVAL=${PVAL_COL} \
	    -v CHR=${CHR_COL} -v A1=${ALLELE1_COL} -v A2=${ALLELE2_COL} -v MAF=${MAF_COL} \
	    -v EFFECT=${EFFECT_COL} -v NUM_COLS=${NUM_COLS} \
	    'BEGIN{OFS=FS; for(i=1; i<=NUM_COLS; i++) {
                   if(i==RSID) {printf "rsID"} else if(i==POS) {printf "pos"}
                   else if(i==PVAL) {printf "pval"} else if(i==CHR) {printf "chr"}
                   else if(i==A1) {printf "allele1"} else if(i==A2) {printf "allele2"}
                   else if(i==MAF) {printf "MAF"} else if(i==EFFECT) {printf "effect"}
                   else {printf i}; if(i!=NUM_COLS) {printf "\t"} else {printf "\n"}}} {print $0}' \
	    ${SUMMARY_STATF} > ${OUTDIR}/summary_stats.header_added.txt

	if [ ${USE_BETA} == "True" ]; then
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${OUTDIR}/summary_stats.header_added.txt \
		--snp_column "rsID" --pvalue_column "pval" \
		--effect_allele_column "allele1" --non_effect_allele_column "allele2" \
		--chromosome_column "chr" --position_column "pos" --beta_column \
		"effect" --freq_column "MAF" --separator $'\t' \
		--covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db 2>&1 | tee ${LOGFILE}
	else
	    ${METAXCAN_DIR}/MetaMany.py --gwas_file ${OUTDIR}/summary_stats.header_added.txt \
		--snp_column "rsID" --pvalue_column "pval" \
		--effect_allele_column "allele1" --non_effect_allele_column "allele2" \
		--chromosome_column "chr" --position_column "pos" --or_column \
		"effect" --freq_column "MAF" --separator $'\t' \
		--covariance_suffix _covariances.txt.gz..opeans_tw_0.5_signif.db \
		--output_directory ${OUTDIR} ${GTEX_V7_DBDIR}/*db 2>&1 | tee ${LOGFILE}
	fi
    fi

    ## after it's done, run the analysis R script
    Rscript ./src/analyze_GTEx_v7_metaXcan.R ${OUTDIR} 2>&1 | tee -a ${LOGFILE}
    
else
    echo "Usage: $0 <metaXcan dir> <output directory> <summary stat file> <GTEx v7 PredictDB directory> <rsid column num> <position column num> <p-value column num> <chromsoome column num> <allele 1 column num> <allele 2 column num> <MAF column num> <beta / OR column num> <true/false for summary stat header> <true/false for beta (true) or OR (false)>" 
fi
	
