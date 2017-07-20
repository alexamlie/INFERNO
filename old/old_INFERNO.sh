
## a master script to run INFERNO analysis on a GWAS dataset
## runs in two modes:
## 1. directly running input variants through the pipeline (does not require summary stats)
## 2. using summary statistics to do p-value expansion as well as colocalization

if [ $# == 4 ]; then
    ## no p-value expansion, LD pruning, or colocalization
    TOP_SNPF=$1
    CFG_FILE=$2
    OUTDIR=$3
    OUTPREFIX=$4
    
    source ${CFG_FILE}

    mkdir -p ${OUTDIR}/logs/
    
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
    time Rscript ./analysis_scripts/Rscript_run_full_analysis.R ./analysis_scripts/ ${PARAMF} "${OUTPREFIX}" TRUE
    
    ## submit a job to run bootstrapping
    bsub -M 40000 -J ${OUTPREFIX}_enh_bootstrapping -o ${OUTDIR}/logs/${OUTPREFIX}_enh_bootstrapping.o%J \
    	-e ${OUTDIR}/logs/${OUTPREFIX}_enh_bootstrapping.e%J ./bsub_wrappers/enhancer_bootstrap_bsub_wrapper.sh \
    	./src/enhancer_only_sample_and_expand_matched_input_variants.R \
    	${NUM_SAMPLES} ${MAF_BIN_SIZE} ${DIST_ROUND} ${DIST_THRESHOLD} ${LD_PARTNER_THRESHOLD} \
    	${BG_SNP_INFOF} ${LD_SETS_DIR} ${REF_SUMMARY_DIR} "${PARAMF}"

elif [[ $# == 15 ]]; then
    RSID_COL=$1
    POS_COL=$2
    PVAL_COL=$3
    CHR_COL=$4
    ALLELE1_COL=$5
    ALLELE2_COL=$6
    MAF_COL=$7
    SIG_MULT=$8
    CASE_PROP=$9
    SAMPLE_SIZE=${10}
    SUMMARY_FILE=${11}
    TOP_SNPF=${12}
    CFG_FILE=${13}
    OUTDIR=${14}
    OUTPREFIX=${15}

    source ${CFG_FILE}

    mkdir -p ${OUTDIR}/logs/

    module load python/2.7.9
    module load bedtools2
    
    ## p-value expansion from the top variant file
    echo "Performing p-value expansion"
    time python data_preprocessing/pval_expand_tagsnp_set.py --single_file --sig_multiplier \
    	${SIG_MULT} --dist_threshold ${LD_AREA} --rsid_col ${RSID_COL} --pos_col ${POS_COL} \
    	--pval_col ${PVAL_COL} --chr_col ${CHR_COL} ${TOP_SNPF} ${SUMMARY_FILE} \
    	${OUTDIR}/${OUTPREFIX}_pval_expanded_snps.txt

    ## now LD prune this
    echo "Performing LD pruning"
    time python data_preprocessing/ld_prune_snp_set.py ${LD_THRESH} \
    	${OUTDIR}/${OUTPREFIX}_pval_expanded_snps.txt ${KG_DIR}/${KG_POP}/ \
    	${OUTDIR}/${OUTPREFIX}_pruning/

    ## now this can go through INFERNO
    echo "Analyzing pruned variants"
    time python ./src/expand_and_annotate_snps.py --loglevel full --kg_pop ${KG_POP} \
    	--ld_threshold ${LD_THRESH} --ld_check_area ${LD_AREA} --gene_bed_file ${GENE_BED_FILE} \
    	--kgxref_file ${KGXREF_FILE} --unstranded_partition_dir ${UNSTRANDED_PARTITION_DIR} \
    	--fantom5_dir ${FANTOM5_DIR} --enhancer_locus_window ${ENH_LOCUS_WINDOW} \
    	--fantom5_correlation_file ${F5_CORR_FILE} --gtex_dir ${GTEX_DIR} --factorbook_file \
    	${FACTORBOOK_FILE} --roadmap_chromhmm_dir ${ROADMAP_HMM_DIR} --homer_motif_bed_file \
    	${HOMER_BED_FILE} --homer_motif_pwm_file ${HOMER_PWM_FILE} --homer_motif_seq_file \
    	${HOMER_SEQ_FILE} --dashr_locus_file ${DASHR_LOCUS_F} --targetscan_dir ${TARGETSCAN_DIR} \
    	${KG_DIR} ${OUTDIR}/${OUTPREFIX}_pruning/pruned_set_pipeline_input.txt \
    	${OUTDIR} ${OUTPREFIX}

    ## run the summarization script
    echo "Summarizing annotation results"
    time ./analysis_scripts/count_annotation_overlaps.sh -l ${ENH_LOCUS_WINDOW} \
    	-f ${F5_CLASSES} -g ${GTEX_CLASSES} -r ${ROADMAP_CLASSES} \
    	${OUTDIR} ${LD_THRESH} ${LD_AREA} ${OUTPREFIX} ${OUTDIR}/summaries/
    
    ## next run the full analysis script
    module load R/3.2.3
    # echo "Running R analysis scripts"
    ## need to get the parameter file from this most recent run
    PARAMF=`ls -t ${OUTDIR}/parameters/*parameters* | head -1`
    ## skip subtitles
    time Rscript ./analysis_scripts/Rscript_run_full_analysis.R ./analysis_scripts/ ${PARAMF} "${OUTPREFIX}" TRUE
    
    ## submit a job to run bootstrapping
    echo "Submitting bootstrapping job"
    bsub -M 40000 -J ${OUTPREFIX}_enh_bootstrapping -o ${OUTDIR}/logs/${OUTPREFIX}_enh_bootstrapping.o%J \
    	-e ${OUTDIR}/logs/${OUTPREFIX}_enh_bootstrapping.e%J ./bsub_wrappers/enhancer_bootstrap_bsub_wrapper.sh \
    	./src/enhancer_only_sample_and_expand_matched_input_variants.R \
    	${NUM_SAMPLES} ${MAF_BIN_SIZE} ${DIST_ROUND} ${DIST_THRESHOLD} ${LD_PARTNER_THRESHOLD} \
    	${BG_SNP_INFOF} ${LD_SETS_DIR} ${REF_SUMMARY_DIR} "${PARAMF}"

    ## submit a job to run co-localization analysis
    echo "Submitting co-localization analysis job"
    bsub -M 40000 -J ${OUTPREFIX}_gtex_colocalization -o ${OUTDIR}/logs/${OUTPREFIX}_gtex_coloc.o%J \
    	-e ${OUTDIR}/logs/${OUTPREFIX}_gtex_coloc.e%J ./bsub_wrappers/gtex_coloc_bsub_wrapper.sh \
    	./src/gtex_gwas_colocalization_analysis.R ${OUTDIR}/gtex_gwas_colocalization_analysis/ \
    	"${PARAMF}" ${COLOC_H4_THRESH} ${COLOC_ABF_THRESH} \
    	${OUTDIR}/${OUTPREFIX}_pruning/pruned_set_pipeline_input.txt ${SUMMARY_FILE} \
    	${COLOC_GTEX_DIR} ${GTEX_SAMPLE_SIZEF} ${GTEX_CLASSES} ${GTEX_RSID_MATCH} \
    	${HG19_ENSEMBL_REF_FILE} "${RELEVANT_CLASSES}" ${RSID_COL} ${POS_COL} ${PVAL_COL} \
    	${CHR_COL} ${ALLELE1_COL} ${ALLELE2_COL} ${MAF_COL} ${CASE_PROP} ${SAMPLE_SIZE}
    
    ## once that's done, do lncRNA correlation analysis
    echo "Submitting lncRNA correlation analysis job"
    bsub -M 40000 -J ${OUTPREFIX}_lncRNA_correlation -o ${OUTDIR}/logs/${OUTPREFIX}_gtex_lncRNA_corr.o%J \
    	-w "done(${OUTPREFIX}_gtex_colocalization)" -e \
    	${OUTDIR}/logs/${OUTPREFIX}_gtex_lncRNA_corr.e%J ./bsub_wrappers/gtex_lncRNA_corr_bsub_wrapper.sh \
    	./src/lncRNA_gtex_correlation.R ${OUTDIR}/gtex_lncRNA_correlation_analysis/ \
    	${OUTDIR}/gtex_gwas_colocalization_analysis/tables/${OUTPREFIX}_gtex_coloc_summaries.txt \
    	${GTEX_EXPR_DIR} ${SAMPLE_INFO_FILE} ${GENCODE_LNCRNA_FILE} ${F5_CLASSES} ${GTEX_CLASSES} \
    	${ROADMAP_CLASSES} ${COLOC_H4_THRESH} ${COR_THRESH}    
else
    echo "This script has two usages: 4 arguments for direct INFERNO analysis, or 15 arguments for p-value expansion, LD pruning, and colocalization"
    echo "Usage: $0 <INFERNO-formatted top SNP file> <Config file> <Output directory> <Output file prefix>"
    echo "Usage: $0 <RSID_COL> <POSITION_COL> <PVAL_COL> <CHR_COL> <ALLELE1_COL> <ALLELE2_COL> <MAF_COL> <PVAL MULTIPLIER> <case/control proportion> <GWAS sample size> <SUMMARY_FILE> <INFERNO-formatted top SNP file> <Config file> <Output directory> <Output file prefix>"
fi
    
