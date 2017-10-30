#!/bin/sh

## all_chr_gtex_expression_filter.sh
## alex amlie-wolf 09/06/2016
## pulls out all the GTEx expression data. requires that the script filter_gtex_rnaseq_genes.py
## is in the same directory as this one!

# DATADIR=/home/alexaml/data/ad_enhancer_analysis/epha1-as1_correlation_analysis/all_chr_gtex_expression/
# GTEX_DATA=/home/alexaml/data/GTEx/rnaseq/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz

if [ $# == 2 ]; then
   DATADIR=$1
   GTEX_DATA=$2

   ## loop through all the gene files and perform the analysis
   for F in ${DATADIR}/chr*_genes.txt; do
       THISCHR=`basename $F | cut -d'_' -f1`
       echo "Parsing ${THISCHR}"
       python filter_gtex_rnaseq_genes.py ${GTEX_DATA} ${F} > ${DATADIR}/${THISCHR}_gtex_expression.txt
   done
fi
