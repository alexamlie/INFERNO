#!/bin/sh

## gtex_lncRNA_corr_bash_wrapper.sh
## alex amlie-wolf 06/20/2018
## a wrapper script to perform lncRNA correlation analysis
## this one is specific for bash analyses due to the complexity of the possible arguments

#    module load R/3.2.3

if [ $# == 12 ]; then
    CORR_SCRIPT=${1}
    LOGFILE=${2}
    OUTDIR=${3}
    COLOC_RESULT_FILE=${4}
    GTEX_EXPR_DIR=${5}
    SAMPLE_INFO_FILE=${6}
    GENCODE_LNCRNA_FILE=${7}
    FANTOM5_CLASS_FILE=${8}
    GTEX_CLASS_FILE=${9}
    ROADMAP_CLASS_FILE=${10}
    COLOC_H4_THRESH=${11} 
    COR_THRESH=${12}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${COR_THRESH} 2>&1 | tee ${LOGFILE}
elif [ $# == 13 ]; then
    CORR_SCRIPT=${1}
    LOGFILE=${2}
    OUTDIR=${3}
    COLOC_RESULT_FILE=${4}
    GTEX_EXPR_DIR=${5}
    SAMPLE_INFO_FILE=${6}
    GENCODE_LNCRNA_FILE=${7}
    FANTOM5_CLASS_FILE=${8}
    GTEX_CLASS_FILE=${9}
    ROADMAP_CLASS_FILE=${10}
    COLOC_H4_THRESH=${11} 
    PEARSON_THRESH=${12}
    SPEARMAN_THRESH=${13}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${PEARSON_THRESH} ${SPEARMAN_THRESH} 2>&1 | tee ${LOGFILE}
elif [ $# == 14 ]; then
    CORR_SCRIPT=${1}
    LOGFILE=${2}
    OUTDIR=${3}
    COLOC_RESULT_FILE=${4}
    GTEX_EXPR_DIR=${5}
    SAMPLE_INFO_FILE=${6}
    GENCODE_LNCRNA_FILE=${7}
    FANTOM5_CLASS_FILE=${8}
    GTEX_CLASS_FILE=${9}
    ROADMAP_CLASS_FILE=${10}
    COLOC_H4_THRESH=${11} 
    PEARSON_THRESH=${12}
    SPEARMAN_THRESH=${13}
    NUM_PCS=${14}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${PEARSON_THRESH} ${SPEARMAN_THRESH} ${NUM_PCS} 2>&1 | tee ${LOGFILE}
else
    echo "Usage: $0 <correlation script> <log file> <outdir> <coloc result file> <GTEx expression dir> <GTEx sample info file> <gencode lncRNA file> <FANTOM5 class file> <GTEx class file> <Roadmap class file> <COLOC P(H_4) threshold> .."
    echo "Last arguments can be threshold on both correlation measures, Pearson and Spearman thresholds (in that order), or Pearson and Spearman thresholds and number of PCs to use for regression"
fi
