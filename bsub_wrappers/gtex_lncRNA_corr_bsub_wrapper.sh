#!/bin/sh

## gtex_lncRNA_corr_bsub_wrapper.sh
## alex amlie-wolf 06/22/2017
## a wrapper script to perform lncRNA correlation analysis

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

#    module load R/3.2.3

    Rscript "${CORR_SCRIPT}" ${OUTDIR} ${COLOC_RESULT_FILE} ${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} \
	${GENCODE_LNCRNA_FILE} ${FANTOM5_CLASS_FILE} ${GTEX_CLASS_FILE} ${ROADMAP_CLASS_FILE} \
	${COLOC_H4_THRESH} ${COR_THRESH}
else
    echo "Usage: $0 <ARGUMENTS!>"
fi
