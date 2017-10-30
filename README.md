#INFERNO
## INFERring the molecular mechanisms of NOncoding genetic variants

### Method description:
The majority of variants identified by genome-wide association studies (GWAS) reside in the
noncoding genome, where they affect regulatory elements including transcriptional enhancers. We
propose INFERNO (INFERring the molecular mechanisms of NOncoding genetic variants), a novel
method which integrates hundreds of diverse functional genomics data sources with GWAS summary
statistics to identify putatively causal noncoding variants underlying association
signals. INFERNO comprehensively characterizes the relevant tissue contexts, target genes, and
downstream biological processes affected by functional variants. 

### Instructions for use:
To use INFERNO, pull the source code from [the bitbucket
repository](https://bitbucket.org/alexamlie/INFERNO/). The full processed annotation datasets
used for the tool are [available for
download](http://tesla.pcbi.upenn.edu/~alexaml/INFERNO/full_INFERNO_annotations.tar.gz). There
is also a [web server](http://inferno.lisanwanglab.org) that runs a subset of the INFERNO
analyses. The full INFERNO pipeline currently assumes you are running it on an interactive node
of a bsub-based cluster system. To extract the annotation data and set up the configuration
file, run these steps:

```bash
## start from the directory containing full_INFERNO_annotations.tar.gz
$ tar -xzvf full_INFERNO_annotations.tar.gz
$ cd full_INFERNO_annotations/
$ ./update_config_file.sh
```

Then, there are two ways to run the pipeline. The more general approach is to make sure that
the following scripts are in your $PATH:

* Python v2.7.9 or a higher version of Python 2.7
* bedtools v2.25.0 or greater version of bedtools v2
* R 3.2.3
* plink v1.90b2i or greater

With this approach, you can run the full INFERNO pipeline using the INFERNO.py script:

```bash
$ python INFERNO.py -h
usage: INFERNO.py [-h] [--skip_ld_expansion] [--rsid_column RSID_COLUMN] [--pos_column POS_COLUMN] [--pval_column PVAL_COLUMN] [--chr_column CHR_COLUMN]
                  [--allele1_column ALLELE1_COLUMN] [--allele2_column ALLELE2_COLUMN] [--maf_column MAF_COLUMN] [--beta_column BETA_COLUMN] [--run_pval_expansion]
                  [--sig_mult SIG_MULT] [--case_prop CASE_PROP] [--sample_size SAMPLE_SIZE] [--summary_file SUMMARY_FILE] [--run_enhancer_sampling] [--run_gtex_coloc]
                  [--run_lncrna_correlation]
                  top_snpf cfg_file outdir outprefix

Driver script for the INFERNO pipeline

positional arguments:
  top_snpf              The tab separated file of the tag SNPs you want to analyze. Should be formatted with four columns: chromosome, rsID, region naming information, and
                        position (in that order). IMPORTANT: Note that SNPs without dbSNP rsIDs should use 'chr-pos' naming format, not 'chr:pos', which is incompatible with
                        this pipeline!
  cfg_file              The configuration file containing paths to all the required functional annotation files. Should be formatted as a bash configuration file i.e.
                        VARIABLE=DEFINTIION on each line.
  outdir                The directory to write all the results to.
  outprefix             The desired prefix for all the output file names

optional arguments:
  -h, --help            show this help message and exit
  --skip_ld_expansion   Give this flag to skip any LD (and p-value-based) expansion and just run analysis directly on the provided list of input variants.
  --rsid_column RSID_COLUMN
                        The summary statistics column number containing the rsIDs
  --pos_column POS_COLUMN
                        The summary statistics column number containing the positions
  --pval_column PVAL_COLUMN
                        The summary statistics column number containing the p-values
  --chr_column CHR_COLUMN
                        The summary statistics column number containing the chromosomes
  --allele1_column ALLELE1_COLUMN
                        The summary statistics column number containing allele 1, which should correspond to the major allele.
  --allele2_column ALLELE2_COLUMN
                        The summary statistics column number containing allele 2, which should correspond to the minor allele.
  --maf_column MAF_COLUMN
                        The summary statistics column number containing the minor allele frequency. Note that if this is ever greater than 0.5 and a beta column is provided,
                        the effect direction will be flipped to be defined relative to the minor allele.
  --beta_column BETA_COLUMN
                        The summary statistics column number containing the beta estimate (used for p-value expansion with consistent directions)
  --run_pval_expansion  If you want to do expansion by p-values when you have summary statistics, provide this flag. Otherwise, the top SNP file will be directly expanded.
  --sig_mult SIG_MULT   The multiplier range for significance of the p-value expanded variant set (e.g. a value of 10 means one order of magnitude)
  --case_prop CASE_PROP
                        The proportion of cases in the GWAS, for colocalization.
  --sample_size SAMPLE_SIZE
                        The total number of samples in the GWAS, for colocalization.
  --summary_file SUMMARY_FILE
                        The path to the full summary statistics file, required for p-value expansion, colocalization analysis, and lncRNA target analysis.
  --run_enhancer_sampling
                        If you want to run the enhancer bootstrapping analysis, provide this flag.
  --run_gtex_coloc      If you want to run COLOC analysis of your summary statistics against GTEx eQTLs from 44 tissues (requires summary statistics)
  --run_lncrna_correlation
                        If you want to analyze expression correlation of any lncRNAs identified by COLOC analysis (--run_gtex_coloc flag) against all other RNAseq-based GTEx
                        genes to find lncRNA targets.
```

For example, the full pipeline call for the schizophrenia analysis in the method manuscript
would look like this using the parsed configuration file from the full annotation file, where
the /path/to/XXX/ directories would be changed to wherever the files are on your server and XXX
refers to the relevant path:

```bash
$ cd /path/to/INFERNO_code/
$ python ./INFERNO.py --rsid_column 2 --pos_column 5 --pval_column 9 --chr_column 1 --allele1_column 3 \
  --allele2_column 4 --maf_column 15 --sig_mult 10 --case_prop 0.2464882 --sample_size 150064 \
  --summary_file /path/to/SCZ_data/scz2.snp.results.1kg_annotations.txt --run_pval_expansion \
  --run_enhancer_sampling --run_gtex_coloc --run_lncrna_correlation \
  /path/to/SCZ_data/SCZ2_128_top_variants_INFERNO_input.no_chrX.tsv \
  /path/to/annotations/full_INFERNO_annotations/INFERNO_annotation_config.cfg \
  /path/to/output_folder/ SCZ2_128_top_variants
```

The other way to run the pipeline is based on your cluster submission system using the 'module'
approach to loading the relevant packages and dependencies. In this case, make sure that the
system contains the following modules:

* python/2.7.9
* bedtools2
* R/3.2.3
* plink/1.90Beta

Then, you can use the INFERNO.sh script to run the pipeline without requiring that the relevant
dependencies are in your $PATH:

```bash
$ cd /path/to/INFERNO_code/
$ ./INFERNO.sh --rsid_column 2 --pos_column 5 --pval_column 9 --chr_column 1 --allele1_column 3 \
  --allele2_column 4 --maf_column 15 --sig_mult 10 --case_prop 0.2464882 --sample_size 150064 \
  --summary_file /path/to/SCZ_data/scz2.snp.results.1kg_annotations.txt --run_pval_expansion \
  --run_enhancer_sampling --run_gtex_coloc --run_lncrna_correlation \
  /path/to/SCZ_data/SCZ2_128_top_variants_INFERNO_input.no_chrX.tsv \
  /path/to/annotations/full_INFERNO_annotations/INFERNO_annotation_config.cfg \
  /path/to/output_folder/ SCZ2_128_top_variants
```

If you wanted to submit the same job in bsub, it would look like this:

```bash
$ cd /path/to/INFERNO_code/
$ bsub -J INFERNO_PGC_SCZ2_ld_pruned -o /path/to/output_folder/SCZ_analysis.o%J -e \
  /path/to/output_folder/SCZ_analysis.e%J \
  ./INFERNO.sh --rsid_column 2 --pos_column 5 --pval_column 9 --chr_column 1 --allele1_column 3 \
  --allele2_column 4 --maf_column 15 --sig_mult 10 --case_prop 0.2464882 --sample_size 150064 \
  --summary_file /path/to/SCZ_data/scz2.snp.results.1kg_annotations.txt --run_pval_expansion \
  --run_enhancer_sampling --run_gtex_coloc --run_lncrna_correlation \
  /path/to/SCZ_data/SCZ2_128_top_variants_INFERNO_input.no_chrX.tsv \
  /path/to/annotations/full_INFERNO_annotations/INFERNO_annotation_config.cfg \
  /path/to/output_folder/ SCZ2_128_top_variants
```

If you want to run INFERNO as a cluster job but your cluster system doesn't use modules, you
can go into the INFERNO.sh script and comment out lines 28-31 and then submit the INFERNO.sh
script to your cluster system. 

### INFERNO output organization:
After the pipeline runs, several folders and files containing output tables and figures are
generated, where outprefix is the last argument to the INFERNO script ("SCZ2_128_top_variants"
in the above example):

* logs/, which contains the main annotation overlap log as well as any bsub output logs
* P-value and LD expansion outputs:
    * outprefix_pval_expanded_snps.txt: the full list of variants expanded by p-value, if this analysis is performed
    * outprefix_pruning/: the results of the LD pruning analysis on the p-value expanded variant set
    * ld_expansion/: the result of the LD expansion on either the LD pruned p-value expanded variants, or the input variants

* From the main INFERNO annotation overlaps (note that not all of these folders will be generated if you are missing some annotation data or skip some steps):
    * parameters/: this folder contains the parameter and specific command run for the annotation overlaps
    * closest_gene/: contains the closest genes to each variant
    * unstranded_genomic_partition/: contains the results of the genomic partition analysis (i.e. annotation of variants overlapping promoters, exons, etc)
    * fantom5_overlap/: variant overlap results with FANTOM5 enhancers
    * closest_fantom5_enh/: contains the closest FANTOM5 enhancers to each variant
    * correlation_enh_targets/: contains the correlation-based target genes of FANTOM5 enhancers with overlapping genetic variants
    * gtex_eqtl_overlap/: the results of the direct eQTL overlap analysis (**NOTE: do not use these results if you have summary statistics and can perform co-localization analysis**)
    * factorbook_overlap/: overlap with FactorBook TFBSs (these are essentially a subset of the HOMER TFBS annotations)
    * homer_motif_overlap/: overlap with the more comprehensive set of HOMER motifs, including PWM calculations
    * roadmap_chromhmm_states/: annotation of variants with Roadmap ChromHMM states across all tissues and cell types
    * dashr_ncrna_loci_overlap/: contains overlaps with [DASHR](http://dashr2.lisanwanglab.org/) small noncoding RNA loci
    * targetscan_miRNA_overlap/: contains overlaps with predicted miRNA binding sites from TargetScan

* summaries/: this folder contains various summaries of annotation overlaps in each tissue class; these files are used for further analysis including enrichment sampling

* analysis_results/: this folder contains all the results from the R analysis script which performs summarization and visualization. This contains the following subfolders, each of which is further split into plots/ and tables/ folders:
    * ld_stats/: visualizations of the LD expansion statistics for this analysis
    * unstranded_genomic_partition/: visualizations of the genomic partition results
    * fantom5_overlap/: visualizations of various characteristics of the FANTOM5 enhancer overlap results
    * closest_fantom5_enhs/: analysis of closest FANTOM5 enhancers
    * gtex_eqtl_overlap/: results from direct eQTL overlap with GTEx data. Again, **do not use these results if you perform co-localization analysis**
    * factorbook_overlap/: analysis of FactorBook TFBS overlaps
    * homer_motif_overlap/: visualizations and analysis of HOMER TFBS results
    * targetscan_miRNA_overlap/: analysis of miRNA binding site overlaps from TargetScan
    * fantom5_roadmap_overlap/: integrative analysis of FANTOM5 and Roadmap enhancer overlaps
    * fantom5_eqtl_chromHMM_overlap/: integrative analysis of FANTOM5 and Roadmap enhancer overlaps with direct GTEx eQTL overlaps

* background_enh_sampling_match_and_expand/: this folder contains the results of the enhancer sampling analysis
    * plots/: contains many visualizations of the sampling results. One important
      characteristic is that the plots with '\_collapsed\_' in the name refer to the sampling
      analysis where multiple overlapping variants in one LD block only contribute one overlap
      count, and file names without it refer to the analysis where multiple variants always are
      counted individually, even if they are in the same LD block. The main visualization used
      for this analysis is
      outprefix_collapsed_annotation_by_tissue_bh_pval_heatmap_with_text.pdf, which contains
      the collapsed analysis with multiple testing corrected p-values. **INFERNO will report
      raw enrichment p-values, but these should not be used for measuring significance without
      multiple testing correction!**
        * input_sampling/: this contains diagnostic plots for the sampling against the input variants
	* ld_collapsed_split_tag_regions/: this folder contains plots for the sampling enrichment analysis using collapsed LD blocks for individual tag regions, with subfolders for each specific region
	* split_tag_regions/: this folder contains plots for the sampling analysis without collapsing by LD blocks for individual tag regions
    * samples/: this (very large) folder contains tables of the specific variants sampled
      against the input as well as expanded variants, and also contains files describing the
      specific LD blocks and tag regions for each sampled variant. These files are mostly just
      index tables for the R script to use and are not very informative to directly look at
    * tables/: this contains easy to parse tables of the outputs of the various steps of the
      sampling analysis, including raw and corrected p-values for both the cross-tag-region
      analysis and the split tag region analyses as well as the LD-collapsed and non-collapsed
      analyses

* gtex_gwas_colocalization_analysis/: this folder contains the results of the eQTL
  co-localization analysis, if that step was run.
    * plots/: this contains several diagnostic and summary plots for the analysis results. The
      main summary plots for the co-localization analysis are
      outprefix_ABF_and_motif_summary_barplot_0.5_prob_thresh.pdf and
      outprefix_ABF_and_motif_nonzero_summary_barplot_0.5_prob_thresh_no_title.pdf, which
      display summaries of the annotation overlaps compared with co-localization results across
      tag regions and different variant prioritization approaches.
        * locuszoom_plots/: if the LOCUSZOOM variable for the path to your locuszoom executable
          is defined in the config file, the colocalization analysis script will use it to
          generate locusZoom plots of the GWAS and GTEx eQTL signals at each strongly
          co-localized locus
    * tables/: the top level of this folder contains tables of various characteristics of the
      co-localization results. These are all tab separated files and can be directly opened in
      Excel for ease of data exploration. The main file containing annotation overlaps and
      co-localization signals for all the ABF-expanded variants is
      outprefix_gtex_coloc_enh_motif_eqtl_info.0.5_thresh_expanded.txt. **Note: Excel will
      often complain if the file name is too long, and this happens quite often with these
      files as they are deep in a nested folder hierarchy. To get around this, just copy the
      full file and name it something like outprefix_COLOC_summary.txt, and then you should be
      able to open it in Excel.**
        * gtex_coloc/: this (very large) folder contains the full co-localization results
          across tag regions (individual subfolders) and GTEx tissues. Typically you wouldn't
          need to look at the specific files here because they are summarized by the analysis
          script.

* gtex_lncRNA_correlation_analysis/: if the correlation analysis to find lncRNA target genes was performed, this folder contains the results.
    * full_correlation_tables/: this contains the full correlation results across all genes in
      the genome for each lncRNA identified in the co-localization analysis. These are the raw
      files and are summarized in a more easily digestible format in the other subfolders from
      this analysis.
    * plots/: this contains summary and diagnostic plots across all lncRNAs as well as
      subfolders for each specific lncRNA containing correlation scatterplots
    * tables/: this contains the lists of target genes identified by correlation threshold for
      each specific lncRNA as well as across all lncRNAs
      (all_lncRNA_genes_0.5_correlation_threshold.txt), as well as lists split by the GTEx
      tissue of the original lncRNA signal (in tissue_specific_gene_lists) and split by the
      tissue category of that original GTEx signal (class_specific_gene_lists). For all of
      these gene lists, if you want to do pathway analysis using
      [WebGestalt](http://webgestalt.org), use the 'genesymbol' option.

### Data sources and pre-processing steps:
#### 1,000 genomes data:
Files downloaded from http://csg.sph.umich.edu/abecasis/mach/download/1000G.2012-03-14.html for
1,000 Genomes data split into 4 populations. Then, several scripts, all available in the
data_preprocessing/ folder of the INFERNO code, are used to process these datasets. All the
steps past the sorting command are for the enhancer sampling analysis. For example, the EUR
population analysis:

```bash
$ cd /path/to/EUR_data/
$ tar -xzvf phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz.tgz
## sort the files
$ mkdir sorted_files/
$ for i in {1..22}; do \
  echo "chr${i}"; \
  zcat chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz | grep "#" > sorted_files/chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf; \
  zcat chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz | grep -v "#" | sort -k2,2n >> sorted_files/chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf; \
done
## calculate minor allele frequencies
$ cd /path/to/INFERNO_code/data_preprocessing/
$ ./calculate_maf.sh /path/to/EUR_data/ chr /path/to/EUR_data/MAF_info
## calculate distance to nearest TSS
$ ./calculate_tss_distance.sh /path/to/EUR_data/ chr /path/to/hg19_ref/hg19_refseq_tss.bed /path/to/bedtools_bin/ /path/to/EUR_data/dist_to_tss/
## calculate pairwise LD (assumes bsub-based cluster system)
$ ./bsub_calculate_pairwise_ld.sh /path/to/EUR_data/ chr /path/to/EUR_data/pairwise_ld/
## can also just run this sequentially, which will take a while:
$ ./calculate_pairwise_ld.sh /path/to/EUR_data/ chr /path/to/EUR_data/pairwise_ld/
## use these LD pairs to precompute LD blocks (in this case, using 0.7 as the R^2 threshold):
$ python compute_ld_sets.py 0.7 /path/to/EUR_data/pairwise_ld/ /path/to/EUR_data/precomputed_ld_sets/
## finally, summarize all of these quantities for each variant, for sampling purposes:
$ python summarize_ld_maf_dist_results.py --ld_threshold 0.7 /path/to/EUR_data/pairwise_ld/ /path/to/EUR_data/MAF_info/ \
  /path/to/EUR_data/dist_to_tss/ /path/to/EUR_data/snp_maf_tss_ld_summary/
```

#### Unstranded genomic partition:
Tables of UCSC genes were downloaded from the UCSC Table Browser, and only chr1-22, X and Y are
used in INFERNO. The 5’ UTR exons and introns, 3’ UTR exons and introns, and exons and introns
were extracted from the knownGene annotation for each protein-coding gene, and all overlapping
exons were merged together. Promoter annotations were defined as 1,000bp upstream of the first
exon in the transcript, either coding or in the UTR.  Variants were then assigned to mutually
exclusive genomic element annotations using the hierarchy: 5’ UTR exon > 5’ UTR intron > 3’ UTR
exon > 3’ UTR intron > promoter > mRNA exon > mRNA intron > repeat. A variant not overlapping
with any class of elements above was classified as intergenic.

#### FANTOM5 enhancers:
Facet-level enhancer expression BED files were
[downloaded](http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz) from
http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz, extracted, and sorted using
the Unix sort tool with arguments ‘-k1,1V -k2,2n’. Specific commands are found in the file
download_and_sort_FANTOM5_files.sh in the data_preprocessing/ directory.

#### Roadmap ChromHMM:
Roadmap 15-state ChromHMM BED files for the 5 core marks (H3K4me3, H3K4me1, H3K36me3, H3K27me3,
H3K9me3) were
[downloaded](http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz)
from
http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz,
extracted, and sorted similarly to the FANTOM5 annotations; commands are available in the file
roadmap_chromhmm_download_and_sort.sh.

#### HOMER TFBSs:
The set of HOMER TFBS annotations was
[downloaded](http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg19.bed.gz) from
http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg19.bed.gz and PWM annotations were
[downloaded](http://homer.ucsd.edu/homer/custom.motifs) from
http://homer.ucsd.edu/homer/custom.motifs. The BED file was extracted and sorted using the same
approach as the FANTOM5 and Roadmap data, and the getfasta tool from the bedtools suite was
used to generate sequences for each BED interval in order to calculate the ∆PWM score, using
the hg19 fasta file downloaded from the UCSC table browser. The specific command used is
available in HOMER_seq_generating_cmds.sh in the data_preprocessing/ folder.

#### GTEx data:
For direct eQTL overlap analysis, files were downloaded and sorted as follows:

```bash
$ cd /path/to/eQTL_dir/
$ wget http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_V6p_eQTLs.tar.gz
$ tar -xzvf GTEx_Analysis_V6p_eQTLs.tar.gz
$ mkdir sorted
$ for f in *.snpgenes; do \
    FNAME=`basename $f`; \
    echo "Sorting file ${FNAME}"; \
    sort -k14,14n -k15,15n $f > sorted/$FNAME; \
    rm $f; \
done
```

For colocalization analysis, the set of all GTEx eQTL tests (not just significant ones) was
[downloaded](http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_all-associations.tar)
from
http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_all-associations.tar,
extracted, and sorted. Specific commands are available in the
gtex_download_and_sort_full_v6p_data.sh in the data preprocessing/ folder.

For the lncRNA correlation analysis, RNA-seq read per kilobase per million (RPKM) values across
all 44 tissues were
[downloaded](http://gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz)
and extracted from
http://gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz.

To identify lncRNA eQTL targets, the GENCODE lncRNA annotations were
[downloaded](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz)
and extracted from
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz.

Finally, the sample sizes for each GTEx tissue were also
[downloaded](http://www.gtexportal.org/home/tissueSummaryPage#sampleCountsPerTissue) as a csv
from http://www.gtexportal.org/home/tissueSummaryPage#sampleCountsPerTissue.

