#!/bin/bash

## gtex_download_and_sort_full_v6p_data.sh
## alex amlie-wolf 04/26/17
## command to get the full GTEx v6p eQTL data, extract it, and sort it

if [ $# == 1 ]; then
    CWD=$(pwd)
    cd $1
    # wget http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_all-associations.tar 
    # tar -xvf GTEx_Analysis_v6p_all-associations.tar
    # rm GTEx_Analysis_v6p_all-associations.tar
    cd GTEx_Analysis_v6p_all-associations/
    mkdir -p sorted.temp/

    for f in *all_snpgene_pairs.txt.gz; do
	FNAME=`basename ${f%.txt.gz}`
	TISS=`basename ${f%_Analysis*}`
	mkdir -p sorted.temp/${TISS}
	# ## just to remove any previous data
	# rm sorted.temp/${TISS}/*
#	echo "Splitting file $FNAME"
#	zcat $f | head -1 > sorted.temp/$FNAME
	## make a separate file for each chromosome
#	zcat $f | tail -n +2 | awk -v OUTDIR="sorted.temp/${TISS}/" -F$'\t' 'BEGIN{OFS=FS} {split($2, SNP_INFO, "_"); print SNP_INFO[1], SNP_INFO[2], $0 >> OUTDIR"/chr"SNP_INFO[1]".txt"}'
	echo "Sorting chromosomal files for $TISS"
	for chrf in sorted.temp/${TISS}/chr*.txt; do	    
	    THIS_CHR=`basename ${chrf} | cut -d'.' -f1`
	    echo -n "${THIS_CHR} "
	    sort -k1,1V -k2,2n ${chrf} | cut -f3-8 | gzip - > sorted.temp/${TISS}/${FNAME}.${THIS_CHR}.txt.gz
	    rm ${chrf}
	done
	echo ""
    done
    cd ${CWD}
fi
