#!/bin/sh

## calculate_roadmap_enh_coverage.sh
## alex amlie-wolf 05/31/16
## just recording the commands I use to check the coverage of the roadmap data

cd ~/data/roadmap/chromHMM/sorted/
mkdir -p ../enhancer_coverage_analysis/

## first have to get all the intervals in one file
## rm -rf ../enhancer_coverage_analysis/all_enhs.bed
for EBED in E*.bed; do
    echo ${EBED}
    ## use awk idioms: it sees we have a condition, and default action on true condition is to print the line
    awk -F$'\t' '$4=="6_EnhG" || $4=="7_Enh" || $4=="12_EnhBiv"' ${EBED} >> ../enhancer_coverage_analysis/all_enhs.bed
done

## now we have to sort this file
cd ../enhancer_coverage_analysis/
time sort -k1,1 -k2,2n all_enhs.bed > all_enhs.sorted.bed

## load bedtools
module load bedtools2
## collapse into non-unique intervals
time bedtools merge -i all_enhs.sorted.bed > all_enhs.merged.bed
## finally, calculate the genome coverage
time bedtools genomecov -i all_enhs.merged.bed -g ~/data/refgenomes/hg19/hg19.chrom.sizes > all_enhs_genomecov.txt

cd -
