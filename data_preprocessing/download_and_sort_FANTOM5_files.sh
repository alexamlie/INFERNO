#!/bin/sh

## download_and_sort_FANTOM5_files.sh
## alex amlie-wolf 
## downloads and sorts the FANTOM5 data files, for use in the enhancer SNP pipeline

if [ $# == 1 ]; then
    FANTOM5_DIR=$1
    mkdir -p ${FANTOM5_DIR}
    cd ${FANTOM5_DIR}
    
    #make directories
    mkdir FANTOM5 FANTOM5/Enhancers 
    mkdir FANTOM5/Enhancers/facet_expressed_enhancers FANTOM5/Enhancers/facet_diffexp_enhancers
    mkdir FANTOM5/Enhancers/enhancer_expression

    #download enhancer lists
    wget http://enhancer.binf.ku.dk/presets/Ubiquitous_enhancers_S9.bed -P FANTOM5/Enhancers 
    wget http://enhancer.binf.ku.dk/presets/Ubiquitous_enhancers_S10.bed -P FANTOM5/Enhancers 
    wget http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed -P FANTOM5/Enhancers 
    wget http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed -P FANTOM5/Enhancers 
    wget http://enhancer.binf.ku.dk/presets/robust_enhancers.bed -P FANTOM5/Enhancers 

    #download tarball containing tissue expressed/differentially expressed enhancers
    wget http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz
    tar xvzf facet_expressed_enhancers.tgz -C FANTOM5/Enhancers/facet_expressed_enhancers 
    rm facet_expressed_enhancers.tgz

    wget http://enhancer.binf.ku.dk/presets/facet_differentially_expressed_0.05.tgz
    tar xvzf facet_differentially_expressed_0.05.tgz -C FANTOM5/Enhancers/facet_diffexp_enhancers
    rm facet_differentially_expressed_0.05.tgz

    mv FANTOM5/Enhancers/facet_diffexp_enhancers/0_05/* FANTOM5/Enhancers/facet_diffexp_enhancers
    rm -r FANTOM5/Enhancers/facet_diffexp_enhancers/0_05

    #download enhancer expression tables
    wget http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancer_usage.csv.gz
    gunzip hg19_permissive_enhancer_usage.csv.gz 
    mv hg19_permissive_enhancer_usage.csv FANTOM5/Enhancers/enhancer_expression

    wget http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancers_expression_rle_tpm.csv.gz
    gunzip hg19_permissive_enhancers_expression_rle_tpm.csv.gz
    mv hg19_permissive_enhancers_expression_rle_tpm.csv FANTOM5/Enhancers/enhancer_expression

    #download the latest version of FANTOM5 data
    # wget -r --no-parent --reject "index.html*" -nH --cut-dirs=2 http://fantom.gsc.riken.jp/5/datafiles/latest/ -P FANTOM5/
    
    mkdir -p ${FANTOM5_DIR}/sorted/

    for BEDF in ${FANTOM5_DIR}/*.bed; do
	echo ${BEDF}
	FNAME=`basename ${BEDF}`
	sort -k1,1V -k2,2n ${BEDF} > ${FANTOM5_DIR}/sorted/${FNAME}
    done
	
else
    echo "Usage: $0 <FANTOM5 dir>"
fi
