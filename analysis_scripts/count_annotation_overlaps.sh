#!/bin/bash

# set -x

## count_annotation_overlaps.sh
## alex amlie-wolf 05/04/2016
## part of the analysis pipeline, generates tables of SNPs and counts for each annotation,
## taking tissue classes into account.

if [[ "$#" -lt 5 ]]; then
    echo "Usage: $0 [-m midpoint_window] [-l locus_window] [-f FANTOM5_CLASSES] [-g GTEX_CLASSES] [-r ROADMAP_CLASSES] <DATA DIRECTORY> <LD THRESHOLD> <DISTANCE THRESHOLD> <ANALYSIS PREFIX> <OUTPUT DIRECTORY>"
else
    ## parse the optional arguments (enhancer overlap approach and class files)
    while getopts "m:l:f:g:r:" key; do
	case "$key" in
	    m)
	    	MIDPOINT_WINDOW="${OPTARG}";;
	    l)
		LOCUS_WINDOW="${OPTARG}";;
	    f)
		F5_FLAG=true; FANTOM5_CLASSES="${OPTARG}";;
	    g)
		GT_FLAG=true; GTEX_CLASSES="${OPTARG}";;
	    r)
		RM_FLAG=true; ROADMAP_CLASSES="${OPTARG}";;
	esac
    done

    ## now get the correct arguments:
    DATADIR=${@:$OPTIND:1}
    LD_THRESH=${@:$OPTIND+1:1}
    DIST_THRESH=${@:$OPTIND+2:1}
    PREFIX=${@:$OPTIND+3:1}
    OUTDIR=${@:$OPTIND+4:1}

    ## define a convenience variable for the parameter suffix of each output file
    PARAM_STRING=${LD_THRESH}_ld_${DIST_THRESH}_dist
    
    mkdir -p ${OUTDIR}
    ## write a convenience file just saying when this data was generated..
    date > ${OUTDIR}/command_run_date.txt
    
    ## create the overlap count file
    echo -e "tissue_class\tcount\tannotation" > ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt
    ## also create one split by tag region and tissue class
    echo -e "tag_region\ttissue_class\tcount\tannotation" > ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt
    ## finally, make one with annotation counts by tag region (across all tissues)
    ## this will also have the number of LD SNPs in each region
    echo -e "tag_region\tannotation\tcount" > ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
    
    ## we also need to get all the unique tissue classes that were provided (start by creating
    ## empty file). this file name doesn't need parameters, because it's always the same
    >${OUTDIR}/all_classes.txt
    
    if ${F5_FLAG}; then
	tail -n +2 ${FANTOM5_CLASSES} | cut -f1 | sort -u >> ${OUTDIR}/all_classes.txt
    fi

    if ${GT_FLAG}; then
	tail -n +2 ${GTEX_CLASSES} | cut -f1 | sort -u >> ${OUTDIR}/all_classes.txt
    fi

    if ${RM_FLAG}; then
	tail -n +2 ${ROADMAP_CLASSES} | cut -f1 | sort -u >> ${OUTDIR}/all_classes.txt
    fi

    sort -u ${OUTDIR}/all_classes.txt > ${OUTDIR}/all_classes_sorted.txt
    mv ${OUTDIR}/all_classes_sorted.txt ${OUTDIR}/all_classes.txt

    ## also make one that has all combinations of classes and tag regions
    tail -n +2 ${DATADIR}/ld_expansion/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}.txt | cut -f10 | sort -u | \
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]; next}
                                    {for(tiss in c) printf "%s\t%s\n", $1, tiss}' \
	     ${OUTDIR}/all_classes.txt - | sort -t$'\t' -k1,1 -k2,2 > ${OUTDIR}/all_classes_tag_regions.txt

    ## finally, make one that has all the tag regions
    tail -n +2 ${DATADIR}/ld_expansion/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}.txt | cut -f10 | sort -u  > ${OUTDIR}/all_tag_regions.txt    
    
    ## count the total number of SNPs
    UNIQ_SNPS=`tail -n +2 ${DATADIR}/ld_expansion/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}.txt | cut -f2 | sort -u | wc -l`
    echo "${UNIQ_SNPS} unique SNPs found by LD expansion" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt

    ## count the number in each tag region (don't need to refer to the tag region because they
    ## will all be in this file)
    tail -n +2 ${DATADIR}/ld_expansion/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}.txt | cut -f2,10 | sort -u | cut -f2 | sort | uniq -c | \
	awk 'BEGIN{OFS="\t"} {print $2, "num_ld_snps", $1}' >> \
	${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
    
    ## number of SNPs with enhancer overlap:
    if [[ ! -z "${MIDPOINT_WINDOW}" && ${F5_FLAG} ]]; then
	if [ ! -e "${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt" ]; then
	    tail -n +2 ${DATADIR}/fantom5_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_${MIDPOINT_WINDOW}bp_around_midpoint_enh_overlaps.txt | \
		cut -f1-3,10,13 | sort -k1,1V -k2,2 | uniq | \
		awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {sub("_expressed_enhancers.bed", "", $2); c[$2]=$1; next}
                             {print $1,$2,$3,$4,$5,c[$5]}' \
	            ${FANTOM5_CLASSES} - > ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt
	fi
	## now count the number of unique overlapping SNPs in each tissue category
	cut -f2,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt \
	    | sort -u | cut -f2 | sort | \
	    ## use awk to count the number of overlaps
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
                         END {for(tiss in c) {printf "%s\t%d\tmidpoint_enh\n", tiss, c[tiss]}}' \
	        ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt
	    
        ## also do this split by tag region and tissue
	cut -f2,4,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt \
	    | sort -u | cut -f2,3 | sort | \
	    ## use awk to count the number of overlaps
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
                         END {for(tiss_tag in c) {printf "%s\t%s\tmidpoint_enh\n", tiss_tag, c[tiss_tag]}}' \
	        ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

        ## also do this split by tag region only
	cut -f2,4 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | uniq -c \
	    | awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                   END {for(tag in c) {printf "%s\tmidpoint_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
	    
	ENH_SNPS=`cut -f1-2 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u | wc -l | cut -d' ' -f1`
	ENH_PROP=`echo "scale=8; ${ENH_SNPS} / ${UNIQ_SNPS}" | bc`
	echo "Number of SNPs overlapping enhancer by midpoint: ${ENH_SNPS} (${ENH_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt
    fi

    if [[ ! -z "${LOCUS_WINDOW}" && ${F5_FLAG} ]]; then
	if [ ! -e "${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt" ]; then
	    tail -n +2 ${DATADIR}/fantom5_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_${LOCUS_WINDOW}bp_around_orig_locus_enh_overlaps.txt | \
		cut -f1-3,10,13 | sort -k1,1V -k2,2 | uniq | \
		awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {sub("_expressed_enhancers.bed", "", $2); c[$2]=$1; next}
                             {print $1,$2,$3,$4,$5,c[$5]}' \
		    ${FANTOM5_CLASSES} - > ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt
	fi
	## now count the number of unique overlapping SNPs in each tissue category
	cut -f2,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | \
	    ## use awk to count the number of overlaps
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {printf "%s\t%d\tlocus_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

        ## also do this split by tag region and tissue class
	cut -f2,4,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt \
	    | sort -u | cut -f2,3 | sort | \
	    ## use awk to count the number of overlaps
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
                         END {for(tiss_tag in c) {printf "%s\t%d\tlocus_enh\n", tiss_tag, c[tiss_tag]}}' \
	        ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

        ## also do this split by tag region only
	cut -f2,4 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | uniq -c \
	    | awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                   END {for(tag in c) {printf "%s\tlocus_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
	    
	ENH_SNPS=`cut -f1-2 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u | wc -l | cut -d' ' -f1`
	ENH_PROP=`echo "scale=8; ${ENH_SNPS} / ${UNIQ_SNPS}" | bc`
	echo "Number of SNPs overlapping extended enhancer locus: ${ENH_SNPS} (${ENH_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt
    fi

    ## number of SNPs with eQTL overlap:
    if [[ ! -e "${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt" && ${GT_FLAG} ]]; then
	tail -n +2 ${DATADIR}/gtex_eqtl_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_eqtl_overlaps.txt | \
	    cut -f1-3,10,13 | sort -k1,1V -k2,2 | uniq | \
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {sub("_Analysis.snpgenes", "", $2); c[$2]=$1; next}
                         {print $1,$2,$3,$4,$5,c[$5]}' ${GTEX_CLASSES} - > ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt
    fi
    ## count the number of unique eQTL SNPs in each tissue category
    cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt \
	| sort -u | cut -f2 | sort | \
	## use awk to count the number of overlaps
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
		     END {for(tiss in c) {printf "%s\t%d\tgtex_eqtl\n", tiss, c[tiss]}}' \
	  ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

    ## count the number of unique eQTL SNPs in each tissue category and tag region
    cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt \
	| sort -u | cut -f2,3 | sort | \
	## use awk to count the number of overlaps
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
		     END {for(tiss_tag in c) {printf "%s\t%d\tgtex_eqtl\n", tiss_tag, c[tiss_tag]}}' \
	  ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

    ## also do this split by tag region only
    cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | uniq -c \
	| awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
	       END {for(tag in c) {printf "%s\tgtex_eqtl\t%s\n", tag, c[tag]}}' \
	  ${OUTDIR}/all_tag_regions.txt - >> \
	  ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
	
    EQTL_SNPS=`cut -f1-2 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u | wc -l | cut -d' ' -f1`
    EQTL_PROP=`echo "scale=8; ${EQTL_SNPS}/${UNIQ_SNPS}" | bc`
    echo "Number of SNPs that are eQTLs by direct overlap: ${EQTL_SNPS} (${EQTL_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt

    ## number of SNPs with roadmap enhancer state
    if [[ ! -e "${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt" && ${RM_FLAG} ]]; then
	tail -n +2 ${DATADIR}/roadmap_chromhmm_states/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_roadmap_chromHMM_states.txt | \
		    awk -F$'\t' 'BEGIN{OFS=FS; j=1}
		    	 NR==FNR {sub("E0*", "", $3); c[$3]=$1; id_map[j]=$3; j+=1; next}
		    	 {for(i=13; i <= NF; i++)
		    	   {if($i=="6_EnhG" || $i=="7_Enh" || $i=="12_EnhBiv")
		    	     printf "%s\t%s\t%s\t%s\tE%03d\t%s\n", $1, $2, $3, $10, id_map[(i-12)], c[id_map[(i-12)]]}}' \
				 <(tail -n +2 ${ROADMAP_CLASSES} | sort -t$'\t' -k3,3V) - > ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt
    fi
    ## count the number of unique roadmap enhancer SNPs in each tissue category
    cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | \
	## use awk to count the number of overlaps
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
		     END {for(tiss in c) {printf "%s\t%d\troadmap_hmm_enh\n", tiss, c[tiss]}}' \
	  ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

    ## count the number of unique roadmap enhancer SNPs in each tissue category and tag region
    cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u | cut -f2,3 | sort | \
	## use awk to count the number of overlaps
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
		     END {for(tiss_tag in c) {printf "%s\t%d\troadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	  ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

    ## also do this split by tag region only
    cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u | cut -f2 | sort | uniq -c \
	| awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
	       END {for(tag in c) {printf "%s\troadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	  ${OUTDIR}/all_tag_regions.txt - >> \
	  ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt
	
    HMM_SNPS=`cut -f1-2 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u | wc -l | cut -d' ' -f1`
    HMM_PROP=`echo "scale=8; ${HMM_SNPS}/${UNIQ_SNPS}" | bc`
    echo "Number of SNPs that overlap a roadmap chromHMM enhancer state: ${HMM_SNPS} (${HMM_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt

    ## number of SNPs with factorbook overlap:
    ## currently, no further analysis is done of these SNPs..
    if [ -e "${DATADIR}/factorbook_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_tfbs_overlaps.txt" ]; then
	if [ ! -e "${OUTDIR}/factorbook_tfbs_snps_${PARAM_STRING}.txt" ]; then
	    tail -n +2 ${DATADIR}/factorbook_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_tfbs_overlaps.txt | cut -f1-2 | sort -k1,1V -k2,2 -u > ${OUTDIR}/factorbook_tfbs_snps_${PARAM_STRING}.txt
	fi
	TFBS_SNPS=`wc -l ${OUTDIR}/factorbook_tfbs_snps_${PARAM_STRING}.txt | cut -d' ' -f1`
	TFBS_PROP=`echo "scale=8; ${TFBS_SNPS} / ${UNIQ_SNPS}" | bc`
	echo "Number of SNPs overlapping FactorBook TFBSs: ${TFBS_SNPS} (${TFBS_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt
    fi
    
    ## number of SNPs with HOMER overlap
    ## first try with the disruption file
    if [ -e "${DATADIR}/homer_motif_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_motif_overlap_and_disruption.txt" ]; then
	if [ ! -e "${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt" ]; then
	    tail -n +2 ${DATADIR}/homer_motif_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_motif_overlap_and_disruption.txt \
		| cut -f1-2 | sort -k1,1V -k2,2 -u > ${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt
	fi
	TFBS_SNPS=`wc -l ${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt | cut -d' ' -f1`
	TFBS_PROP=`echo "scale=8; ${TFBS_SNPS} / ${UNIQ_SNPS}" | bc`
	echo "Number of SNPs overlapping HOMER-predicted TFBSs: ${TFBS_SNPS} (${TFBS_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt
    ## if not, try the overlap file
    elif [ -e "${DATADIR}/homer_motif_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_motif_overlap.txt" ]; then
	if [ ! -e "${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt" ]; then
	    tail -n +2 ${DATADIR}/homer_motif_overlap/${PREFIX}_${LD_THRESH}_ld_cutoff_snps_within_${DIST_THRESH}_motif_overlap.txt \
		| cut -f1-2 | sort -k1,1V -k2,2 -u > ${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt
	fi
	TFBS_SNPS=`wc -l ${OUTDIR}/homer_tfbs_snps_${PARAM_STRING}.txt | cut -d' ' -f1`
	TFBS_PROP=`echo "scale=8; ${TFBS_SNPS} / ${UNIQ_SNPS}" | bc`
	echo "Number of SNPs overlapping HOMER-predicted TFBSs: ${TFBS_SNPS} (${TFBS_PROP} of total)" | tee ${OUTDIR}/${PREFIX}_annotation_overlap_summary.txt
    fi
    
    ## now look at combination overlap counts by category
    ## fantom5 enhancer and eQTL overlap
    if [[ ! -z "${MIDPOINT_WINDOW}" && ${F5_FLAG} && ${GT_FLAG} ]]; then
	comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tmid_enh+gtex_eqtl\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt
	## also split by tag region and tissue class
	comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                  printf "%s\t%d\tmid_enh+gtex_eqtl\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt
	    
	## finally, do this split by tag region only
	comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c |\
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tmid_enh+gtex_eqtl\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    if [[ ! -z "${LOCUS_WINDOW}" && ${F5_FLAG} && ${GT_FLAG} ]]; then
	comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                  printf "%s\t%d\tlocus_enh+gtex_eqtl\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

	## also split by tag region and tissue class
	comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tlocus_enh+gtex_eqtl\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## finally, do this split by tag region only
	comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c |\
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tlocus_enh+gtex_eqtl\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    ## fantom5 enhancer and roadmap overlap
    if [[ ! -z "${MIDPOINT_WINDOW}" && ${F5_FLAG} && ${RM_FLAG} ]]; then
	comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tmid_enh+roadmap_hmm_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

	## also split by tag region and tissue class
	comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tmid_enh+roadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## finally, do this split by tag region only
	comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c |\
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tmid_enh+roadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    if [[ ! -z "${LOCUS_WINDOW}" && ${F5_FLAG} && ${RM_FLAG} ]]; then
	comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tlocus_enh+roadmap_hmm_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

	## also split by tag region and tissue class
	comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tlocus_enh+roadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## finally, do this split by tag region only
	comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c | \
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tlocus_enh+roadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    ## eQTL and roadmap overlap
    if [[ ${GT_FLAG} && ${RM_FLAG} ]]; then
	comm -1 -2 <(cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tgtex_eqtl+roadmap_hmm_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

	## also split by tag region and tissue class
	comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tgtex_eqtl+roadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## finally, do this split by tag region only
	comm -1 -2 <(cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u) \
	    <(cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c | \
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tgtex_eqtl+roadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    ## fantom5, eqtl, and roadmap overlap
    if [[ ! -z "${MIDPOINT_WINDOW}" && ${F5_FLAG} && ${GT_FLAG} && ${RM_FLAG} ]]; then
	comm -1 -2 <(comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tmid_enh+gtex_eqtl+roadmap_hmm_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt
	    
	## split by tag region and tissue class
	comm -1 -2 <(comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tmid_enh+gtex_eqtl+roadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## split by tag region only
	comm -1 -2 <(comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_mid_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c | \
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tmid_enh+gtex_eqtl+roadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    if [[ ! -z "${LOCUS_WINDOW}" && ${F5_FLAG} && ${GT_FLAG} && ${RM_FLAG} ]]; then
	comm -1 -2 <(comm -1 -2 <(cut -f2,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1]=0; next} {c[$1]+=1}
			 END {for(tiss in c) {
                                 printf "%s\t%d\tlocus_enh+gtex_eqtl+roadmap_hmm_enh\n", tiss, c[tiss]}}' \
	      ${OUTDIR}/all_classes.txt - >> ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

	## split by tag region and tissue class
	comm -1 -2 <(comm -1 -2 <(cut -f2,4,6 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,4,6 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,4,6 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2,3 | \
	    ## count the overlap number
	    awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {c[$1"\t"$2]=0; next} {c[$1"\t"$2]+=1}
			 END {for(tiss_tag in c) {
                                 printf "%s\t%d\tlocus_enh+gtex_eqtl+roadmap_hmm_enh\n", tiss_tag, c[tiss_tag]}}' \
	      ${OUTDIR}/all_classes_tag_regions.txt - >> ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

	## split by tag region only
	comm -1 -2 <(comm -1 -2 <(cut -f2,4 ${OUTDIR}/enh_locus_snps_${PARAM_STRING}.txt | sort -u) \
	                        <(cut -f2,4 ${OUTDIR}/eqtl_snps_${PARAM_STRING}.txt | sort -u)) \
	    <(cut -f2,4 ${OUTDIR}/roadmap_hmm_snps_${PARAM_STRING}.txt | sort -u) | cut -f2 | sort | uniq -c | \
	    ## count the overlap number
	    awk 'BEGIN{OFS="\t"} NR==FNR {c[$1]=0; next} {c[$2]=$1}
                 END {for(tag in c) {printf "%s\tlocus_enh+gtex_eqtl+roadmap_hmm_enh\t%s\n", tag, c[tag]}}' \
	      ${OUTDIR}/all_tag_regions.txt - >> \
	      ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt	    
    fi

    ## for convenience, sort the annotation count files
    head -1 ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt > ${OUTDIR}/tag_region_counts_sorted.txt
    tail -n +2 ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt | sort -t$'\t' -k1,1 -k2,2 -k4,4 >> ${OUTDIR}/tag_region_counts_sorted.txt
    mv ${OUTDIR}/tag_region_counts_sorted.txt ${OUTDIR}/tag_region_tissue_class_annotation_counts_${PARAM_STRING}.txt

    head -1 ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt > ${OUTDIR}/tissue_counts_sorted.txt
    tail -n +2 ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt | sort -t$'\t' -k1,1 -k3,3 >> ${OUTDIR}/tissue_counts_sorted.txt
    mv ${OUTDIR}/tissue_counts_sorted.txt ${OUTDIR}/tissue_class_annotation_counts_${PARAM_STRING}.txt

    head -1 ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt > ${OUTDIR}/tag_region_annot_counts_sorted.txt
    tail -n +2 ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt | sort -t$'\t' -k1,1 -k2,2 >> ${OUTDIR}/tag_region_annot_counts_sorted.txt
    mv ${OUTDIR}/tag_region_annot_counts_sorted.txt ${OUTDIR}/tag_region_annotation_counts_${PARAM_STRING}.txt    
fi
