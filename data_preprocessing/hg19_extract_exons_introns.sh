#INFILE=hg38_ucsc_refGene 
#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
#0	NM_001276351	chr1	-	67092175	67134971	6709300467127240	8	67092175,67095234,67096251,67115351,67125751,67127165,67131141,67134929,	67093604,67095421,67096321,67115464,67125909,67127257,67131227,67134971,	0	C1orf141	cmpl	cmpl	0,2,1,2,0,0,-1,-1,

INFILE=hg19_ucsc_refGene_with_kgXref
#hg38.refGene.bin	hg38.refGene.name	hg38.refGene.chrom	hg38.refGene.strand	hg38.refGene.txStart	hg38.refGene.txEnd	hg38.refGene.cdsStart	hg38.refGene.cdsEnd	hg38.refGene.exonCount	hg38.refGene.exonStartshg38.refGene.exonEnds	hg38.refGene.score	hg38.refGene.name2	hg38.refGene.cdsStartStat	hg38.refGene.cdsEndStat	hg38.refGene.exonFrames	hg38.kgXref.kgID	hg38.kgXref.mRNA	hg38.kgXref.spID	hg38.kgXref.spDisplayID	hg38.kgXref.geneSymbol	hg38.kgXref.refseq	hg38.kgXref.protAcc	hg38.kgXref.description	hg38.kgXref.rfamAcc	hg38.kgXref.tRnaName
#0	NM_001276351	chr1	-	67092175	67134971	6709300467127240	8	67092175,67095234,67096251,67115351,67125751,67127165,67131141,67134929,	67093604,67095421,67096321,67115464,67125909,67127257,67131227,67134971,	0	C1orf141	cmpl	cmpl	0,2,1,2,0,0,-1,-1,	uc057hjb.1	NM_001276351	Q5JVX7	CA141_HUMAN	C1orf141	NM_001276351	NM_001276351	Homo sapiens chromosome 1 open reading frame 141 (C1orf141), transcript variant 1, mRNA. (from RefSeq NM_001276351)

# extract CDS
# cdsStart != cdsEnd
awk 'BEGIN{OFS="\t"}{if (NR==1) next; cdsStart=$7; cdsEnd=$8; txStart=$5; txEnd=$6; if (cdsStart!=cdsEnd) print;}' ${INFILE} > ${INFILE}.cds

awk 'BEGIN{OFS="\t"}
       {
       chr = $3;
       strand = $4;
       parent = $2;
       ucscID = $17;
       uniprotID = $19;
       geneSymbol = $21;
       refseqID = $22;
       gsub(/,$/,"",ucscID);
       gsub(/,$/,"",uniprotID);
       gsub(/,$/,"",geneSymbol);
       gsub(/,$/,"",refseqID);
       parent=($2 ";ucscID=" ucscID ";uniprotID=" uniprotID ";refseqID=" refseqID ";geneSymbol=" geneSymbol);
       nExons=$9+0;
       cdsStart=$7+1; # 1-based
       cdsEnd=$8+0; # 1-based
       txStart = $5;
       txEnd = $6;
       split($10, exonStarts, ",");
       split($11, exonEnds, ",");
       nIntrons = nExons - 1;
       #print parent, nExons, nIntrons, cdsStart, cdsEnd
       for (i=1; i<=nIntrons; ++i)
       {
         intronStarts[i] = exonEnds[i]+1; # 0-based
         intronEnds[i] = exonStarts[i+1]-1; # 0-based
         #print "intron", i, intronStarts[i], intronEnds[i]
       }

       if (txStart!=exonStarts[1] )
       {
         print "ERROR:", parent, txStart, "!=", exonStarts[1]
         exit;
       }
       if (txEnd!=exonEnds[nExons])
       {
         print "ERROR:", parent, txEnd, "!=", exonEnds[nExons]
         exit;
       }
 

       for (i=1; i<=nIntrons; ++i)
       {
         intStart = intronStarts[i]+0; # 0-based
         intEnd = intronEnds[i]+1; # 1-based
         if (intEnd<=intStart) {continue}
         if (sqrt((intEnd-intStart)*(intEnd-intStart)) <= 1) {continue}
         caseType="0";
         if (intEnd < cdsStart)
         {
            print chr,"ucsc_refGene","UTR5_intron",intStart,intEnd,".",strand,".","Parent="parent 
            caseType = "Left";
         }
         else if (intStart<cdsStart && intEnd < cdsEnd)
         {
            caseType = "cdsStartCross";
            print chr,"ucsc_refGene","UTR5_intron", intStart, cdsStart-1,".",strand,".","Parent="parent 
            print chr,"ucsc_refGene","mRNA_intron", cdsStart, intEnd,".",strand,".","Parent="parent
         }
         else if (intStart>=cdsStart && intStart < cdsEnd && intEnd > cdsEnd) {
            caseType = "cdsEndCross";
            print chr,"ucsc_refGene","mRNA_intron",intStart,cdsEnd,".",strand,".","Parent="parent 
            print chr,"ucsc_refGene","UTR3_intron",cdsEnd+1,intEnd,".",strand,".","Parent="parent 
         } else if (intStart > cdsEnd) {
            caseType="Right";
            print chr,"ucsc_refGene","UTR3_intron",intStart,intEnd,".",strand,".","Parent="parent 
         } else {
            caseType="Inside";
            print chr,"ucsc_refGene","mRNA_intron",intStart,intEnd,".",strand,".","Parent="parent 
         }
         #print "Intron ", i, parent, caseType;
       }       

      for (i=1; i<=nExons; ++i)
       { exStart = exonStarts[i]+1; # since UCSC start is 0-based
         exEnd = exonEnds[i]+0; # since UCSC end is 1-based
         exType="mRNA_exon";
         if (exEnd <= exStart) {continue}
         if (exEnd < cdsStart) {
            caseType = "Left";
            print chr,"ucsc_refGene", "UTR5_exon", exStart, exEnd,".",strand,".","Parent="parent
         } else if (exStart < cdsStart && exEnd < cdsEnd) {
                caseType = "cdsStartCross";
                # split into UTR and coding part
                print chr,"ucsc_refGene","UTR5_exon", exStart, cdsStart-1,".",strand,".","Parent="parent
                print chr,"ucsc_refGene", "mRNA_exon", cdsStart, exEnd,".",strand,".","Parent="parent
         } else if (exStart >= cdsStart && exStart < cdsEnd && exEnd > cdsEnd) {
                  caseType = "cdsEndCross";
                  print chr,"ucsc_refGene","mRNA_exon", exStart, cdsEnd,".",strand,".","Parent="parent
                  print chr,"ucsc_refGene","UTR3_exon", cdsEnd+1, exEnd,".",strand,".","Parent="parent
         } else if (exStart < cdsStart && exEnd > cdsEnd) {
                  caseType = "cdsStartAndEndCross"
                  print chr,"ucsc_refGene","UTR5_exon", exStart, cdsStart-1,".",strand,".","Parent="parent
                  print chr,"ucsc_refGene","mRNA_exon", cdsStart, cdsEnd,".",strand,".","Parent="parent
                  print chr,"ucsc_refGene","UTR3_exon", cdsEnd+1, exEnd,".",strand,".","Parent="parent
         } else if (exStart > cdsEnd) {
                  caseType = "Right"
                  print chr,"ucsc_refGene","UTR3_exon", exStart, exEnd,".",strand,".","Parent="parent
         } else { caseType = "Inside"; print chr,"ucsc_refGene","mRNA_exon", exStart, exEnd,".",strand,".","Parent="parent } 
         #print "Exon", i, parent, caseType;
       }
     }' ${INFILE}.cds | sort -k1,1 -k4,4n -k5,5n -k7,7 > ${INFILE}.exons_introns_utrs
    awk 'BEGIN{OFS="\t"}{h=($1 "\t" $4 "\t" $5 "\t" $7); if (!a[h]++) print}' ${INFILE}.exons_introns_utrs > ${INFILE}.exons_introns_utrs.unique_LOC

#| awk '$5-$4>10' # output only exons/introns with length > 10  
