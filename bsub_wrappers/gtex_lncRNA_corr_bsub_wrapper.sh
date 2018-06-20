#!/bin/sh

## gtex_lncRNA_corr_bsub_wrapper.sh
## alex amlie-wolf 06/22/2017
## a wrapper script to perform lncRNA correlation analysis

#    module load R/3.2.3

if [ $# == 11 ]; then
    CORR_SCRIPT=${1}
    OUTDIR=${2}
    COLOC_RESULT_FILE=${3}
    GTEX_EXPR_DIR=${4}
    SAMPLE_INFO_FILE=${5}
    GENCODE_LNCRNA_FILE=${6}
    FANTOM5_CLASS_FILE=${7}
    GTEX_CLASS_FILE=${8}
    ROADMAP_CLASS_FILE=${9}
    COLOC_H4_THRESH=${10} 
    COR_THRESH=${11}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${COR_THRESH}
elif [ $# == 12 ]; then
    CORR_SCRIPT=${1}
    OUTDIR=${2}
    COLOC_RESULT_FILE=${3}
    GTEX_EXPR_DIR=${4}
    SAMPLE_INFO_FILE=${5}
    GENCODE_LNCRNA_FILE=${6}
    FANTOM5_CLASS_FILE=${7}
    GTEX_CLASS_FILE=${8}
    ROADMAP_CLASS_FILE=${9}
    COLOC_H4_THRESH=${10} 
    PEARSON_THRESH=${11}
    SPEARMAN_THRESH=${12}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${PEARSON_THRESH} ${SPEARMAN_THRESH}
elif [ $# == 13 ]; then
    CORR_SCRIPT=${1}
    OUTDIR=${2}
    COLOC_RESULT_FILE=${3}
    GTEX_EXPR_DIR=${4}
    SAMPLE_INFO_FILE=${5}
    GENCODE_LNCRNA_FILE=${6}
    FANTOM5_CLASS_FILE=${7}
    GTEX_CLASS_FILE=${8}
    ROADMAP_CLASS_FILE=${9}
    COLOC_H4_THRESH=${10} 
    PEARSON_THRESH=${11}
    SPEARMAN_THRESH=${12}
    NUM_PCS=${13}

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${PEARSON_THRESH} ${SPEARMAN_THRESH} ${NUM_PCS}
else
    echo "Usage: $0 <correlation script> <outdir> <coloc result file> <GTEx expression dir> <GTEx sample info file> <gencode lncRNA file> <FANTOM5 class file> <GTEx class file> <Roadmap class file> <COLOC P(H_4) threshold> .."
    echo "Last arguments can be threshold on both correlation measures, Pearson and Spearman thresholds (in that order), or Pearson and Spearman thresholds and number of PCs to use for regression"
fi
