#!/bin/bash

## INFERNO.sh
## Alex Amlie-Wolf
## a master script to run INFERNO analysis on a GWAS dataset

if [ $# == 4 ]; then
    TOP_SNPF=$1
    CFG_FILE=$2
    OUTDIR=$3
    OUTPREFIX=$4
    
    source ${CFG_FILE}

    module load python-2.7.9
    python src/expand_and_annotate_snps.py --loglevel full --kg_pop ${KG_POP} \
	--ld_threshold ${LD_THRESH} --ld_check_area ${LD_AREA} --gene_bed_file ${GENE_BED_FILE} \
	--kgxref_file ${KGXREF_FILE} --unstranded_partition_dir ${UNSTRANDED_PARTITION_DIR} \
	--fantom5_dir ${FANTOM5_DIR} --enhancer_locus_window ${ENH_LOCUS_WINDOW} \
	--fantom5_correlation_file ${F5_CORR_FILE} --gtex_dir ${GTEX_DIR} --factorbook_file \
	${FACTORBOOK_FILE} --roadmap_chromhmm_dir ${ROADMAP_HMM_DIR} --homer_motif_bed_file \
	${HOMER_BED_FILE} --homer_motif_pwm_file ${HOMER_PWM_FILE} --homer_motif_seq_file \
	${HOMER_SEQ_FILE} ${KG_DIR} ${TOP_SNPF} ${OUTDIR} ${OUTPREFIX}
else
    echo "Usage: $0 <INFERNO-formatted top SNP file> <Config file> <Output directory> <Output file prefix"
fi
    
