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
    mkdir -p sorted/

    for f in *all_snpgene_pairs.txt.gz; do
	FNAME=`basename ${f%.gz}`
	echo "Sorting file ${FNAME%.gz}"
	zcat $f | head -1 > sorted/$FNAME
	zcat $f | tail -n +2 | awk -F$'\t' 'BEGIN{OFS=FS} {split($2, SNP_INFO, "_"); print SNP_INFO[1], SNP_INFO[2], $0}' | sort -k1,1n -k2,2n | cut -f3-8 >> sorted/$FNAME
	## now compress the target
	gzip sorted/${FNAME}
    done
    cd ${CWD}
fi
