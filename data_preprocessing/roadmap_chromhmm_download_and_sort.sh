## command to download the roadmap ChromHMM files and sort them

wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
tar -xzvf all.mnemonics.bedFiles.tgz

mkdir -p sorted

for BEDF in *.bed.gz; do
    echo ${BEDF}
    zcat ${BEDF} | sort -k1,1V -k2,2n > sorted/${BEDF%.gz}
done
