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

    mkdir -p ${OUTDIR}/
    
    module load python/2.7.9
    module load bedtools2
    
    ## first run the annotation script
    time python ./src/expand_and_annotate_snps.py --loglevel full --kg_pop ${KG_POP} \
    	--ld_threshold ${LD_THRESH} --ld_check_area ${LD_AREA} --gene_bed_file ${GENE_BED_FILE} \
    	--kgxref_file ${KGXREF_FILE} --unstranded_partition_dir ${UNSTRANDED_PARTITION_DIR} \
    	--fantom5_dir ${FANTOM5_DIR} --enhancer_locus_window ${ENH_LOCUS_WINDOW} \
    	--fantom5_correlation_file ${F5_CORR_FILE} --gtex_dir ${GTEX_DIR} --factorbook_file \
    	${FACTORBOOK_FILE} --roadmap_chromhmm_dir ${ROADMAP_HMM_DIR} --homer_motif_bed_file \
    	${HOMER_BED_FILE} --homer_motif_pwm_file ${HOMER_PWM_FILE} --homer_motif_seq_file \
    	${HOMER_SEQ_FILE} --dashr_locus_file ${DASHR_LOCUS_F} --targetscan_dir ${TARGETSCAN_DIR} \
    	${KG_DIR} ${TOP_SNPF} ${OUTDIR} ${OUTPREFIX}

    ## run the summarization script
    time ./analysis_scripts/count_annotation_overlaps.sh -l ${ENH_LOCUS_WINDOW} \
    	-f ${F5_CLASSES} -g ${GTEX_CLASSES} -r ${ROADMAP_CLASSES} \
    	${OUTDIR} ${LD_THRESH} ${LD_AREA} ${OUTPREFIX} ${OUTDIR}/summaries/
    
    ## next run the full analysis script
    module load R/3.2.3
    ## need to get the parameter file from this most recent run
    PARAMF=`ls -t ${OUTDIR}/parameters/*parameters* | head -1`
    ## skip subtitles
    time Rscript ./analysis_scripts/Rscript_run_full_analysis.R ./analysis_scripts/ ${PARAMF} ${OUTPREFIX} TRUE
    
    ## submit a job to run bootstrapping
    bsub -M 40000 -J ${OUTPREFIX}_enh_bootstrapping -o ${OUTDIR}/logs/enh_bootstrapping.o%J \
    	-e ${OUTDIR}/logs/enh_bootstrapping.e%J ./bsub_wrappers/enhancer_bootstrap_bsub_wrapper.sh \
    ./src/enhancer_only_sample_and_expand_matched_input_variants.R \
    	${NUM_SAMPLES} ${MAF_BIN_SIZE} ${DIST_ROUND} ${DIST_THRESHOLD} ${LD_PARTNER_THRESHOLD} \
    	${BG_SNP_INFOF} ${LD_SETS_DIR} ${REF_SUMMARY_DIR} "${PARAMF}"

    ## submit a job to run generic colocalization
#    bsub -M 40000 ./src/generic_gtex_colocalization_analysis.R
else
    echo "Usage: $0 <INFERNO-formatted top SNP file> <Config file> <Output directory> <Output file prefix"
fi
    
