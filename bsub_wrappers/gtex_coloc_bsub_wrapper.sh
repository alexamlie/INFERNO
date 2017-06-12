#!/bin/sh

## gtex_coloc_bsub_wrapper.sh
## alex amlie-wolf 06/02/17
## a wrapper script to submit the GTEx colocalization script to bsub

if [ $# == 22 ]; then
    COLOC_SCRIPT=${1}
    OUTDIR=${2}
    INFERNO_PARAM_F=${3}
    COLOC_H4_THRESH=${4}
    COLOC_ABF_THRESH=${5}    
    TOP_SNPF=${6}
    GWAS_SUMMARY_FILE=${7}
    GTEX_DIR=${8}
    GTEX_SAMPLE_SIZEF=${9}
    GTEX_CLASS_FILE=${10}
    GTEX_RSID_MATCH=${11}
    HG19_ENSEMBL_REF_FILE=${12}
    RELEVANT_CLASSES=${13}
    RSID_COL=${14}
    POS_COL=${15}
    PVAL_COL=${16}
    CHR_COL=${17}
    ALLELE1_COL=${18}
    ALLELE2_COL=${19}
    MAF_COL=${20}
    CASE_PROP=${21}
    SAMPLE_SIZE=${22}
        
    module load R/3.2.3
    
    Rscript "${COLOC_SCRIPT}" ${OUTDIR} ${INFERNO_PARAM_F} ${COLOC_H4_THRESH} ${COLOC_ABF_THRESH} \
	${TOP_SNPF} ${GWAS_SUMMARY_FILE} ${GTEX_DIR} ${GTEX_SAMPLE_SIZEF} ${GTEX_CLASS_FILE} \
	${GTEX_RSID_MATCH} ${HG19_ENSEMBL_REF_FILE} "${RELEVANT_CLASSES}" ${RSID_COL} ${POS_COL} \
	${PVAL_COL} ${CHR_COL} ${ALLELE1_COL} ${ALLELE2_COL} ${MAF_COL} ${CASE_PROP} ${SAMPLE_SIZE}
else
    echo "Usage: $0 <FILL ME IN!>"
fi
