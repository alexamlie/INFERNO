mkdir -p bed_files

for EQTL_F in sorted/*snpgenes; do
    FNAME=`basename ${EQTL_F}`
    TISS=${FNAME%_Analysis.snpgenes}
    tail -n +2 ${EQTL_F} | awk -F$'\t' -v TISS=${TISS} 'BEGIN{OFS=FS} {printf "chr%s\t%d\t%d\t%s_%s\n", $14, $15, $15+1, TISS, $27}' > bed_files/${TISS}.bed
done
