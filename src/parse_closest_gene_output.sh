#!/bin/bash

## parse_closest_gene_output.sh
## alex amlie-wolf
## used to parse the output of the bedtools wrapper script

RAW_OUTPUT_F=$1

sed -e "s/:/\t/; s/:/\t/" ${RAW_OUTPUT_F} | sort -k1,1V -k2,2n -k4,4 -k6,6 -k5,5 | \
    awk -F$'\t' 'BEGIN{OFS=FS} {$4=$4":"$5":"$6; print $0}' | \
    cut -f5,6 --complement > ${RAW_OUTPUT_F}.temp
mv ${RAW_OUTPUT_F}.temp ${RAW_OUTPUT_F}
