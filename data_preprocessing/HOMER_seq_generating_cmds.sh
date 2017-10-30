module load bedtools2
## we need to subtract one from the start to get the full correct sequence
time bedtools getfasta -fi ~/data/refgenomes/hg19/hg19.fa -bed <(awk -F$'\t' 'BEGIN{OFS=FS} {$2=$2-1; print $0}' homer.sorted.KnownMotifs.hg19.bed) -s -tab -fo homer.sorted.KnownMotifs.hg19.sequence.txt
