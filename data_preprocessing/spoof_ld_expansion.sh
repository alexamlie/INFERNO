#!/bin/sh

## alex amlie-wolf
## recording the command used to 'spoof' the LD expansion script

INFILE=$1
OUTFILE=$2

echo -e "chr\trsID\tpos\tref\talt\tMAF\ttag_rsID\ttag_pos\ttag_MAF\ttag_name\tR2\tDprime" > ${OUTFILE}

tail -n +2 ${INFILE} | awk -v OFS=$'\t' '{if($2==".") {$2=$1"-"$3} print $1, $2, $3, $7, $8, $5, $2, $3, $5, "NA", "1.0", "1.0"}' - >> ${OUTFILE}
