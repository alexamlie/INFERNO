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
#### Obtaining source code and annotations:
To use INFERNO, pull the source code from [the bitbucket
repository](https://bitbucket.org/wanglab-upenn/INFERNO/). The full processed annotation datasets
used for the tool are [available for
download](https://lisanwanglab.org/FILER/full_INFERNO_annotations.tar.gz). There
is also a [web server](http://inferno.lisanwanglab.org) that runs a subset of the INFERNO
analyses. The full INFERNO pipeline can be run on a bsub-based cluster system (default) or through direct execution, using the '--cluster_system' flag. To extract the annotation data and set up the configuration
file, run these steps:

```bash
## start from the directory containing full_INFERNO_annotations.tar.gz
$ tar -xzvf full_INFERNO_annotations.tar.gz
$ cd full_INFERNO_annotations/
$ ./update_config_file.sh
```

#### Options in configuration file:
To use the INFERNO pipeline, some arguments are given to the command line tool (see below), and some are set up in the configuration file, mostly for setting up where various data files are:

Config file variable | Value
-------------------- | -----
**Parameters for annotation overlaps** |
KG_POP		     | Desired 1,000 genomes population to use (EUR, AMR, ASN, AFR), EUR by default
LD_THRESH	     | Threshold for R^2 values of the LD expansion (Default = 0.7)
LD_AREA		     | Distance around each tag variant to check (Default = 500000)
KG_DIR		     | The folder containing the population-specific 1,000 genomes data (vcf files, precomputed LD sets, and summary files)
GENE_BED_FILE	     | Bed file containing exons of protein-coding genes
KGXREF_FILE	     | Reference file to match exon IDs to gene names
UNSTRANDED_PARTITION_DIR	 | Directory containing parsed partition information
FANTOM5_DIR			 | Directory containing sorted FANTOM5 facet-level enhancer bed files
ENHANCER_LOCUS_WINDOW		 | Base-pair window around FANTOM5 enhancer loci to look for overlaps (Default = 1000)
FANTOM5_CORRELATION_FILE	 | File containing the correlation-based TSS targets of FANTOM5 enhancers
ROADMAP_CHROMHMM_DIR		 | Directory containing the sorted bed files of ChromHMM states across Roadmap samples
BEDTOOLS_BIN_DIR		 | Specific path to bedtools, or can be left undefined if you have bedtools in your path
GTEX_DIR			 | Directory containing the sorted files of significant eQTL signals (for direct overlap)
FACTORBOOK_FILE			 | If desired, the file containing the FactorBook TFBS annotations
HOMER_MOTIF_BED_FILE		 | The bed file containing the positions of HOMER-identified TFBSs
HOMER_MOTIF_PWM_FILE		 | The custom.motifs file containing the PWMs for each TF analyzed by HOMER
HOMER_MOTIF_SEQ_FILE		 | The processed file containing the sequences of each HOMER TFBS, for delta PWM calculation
DASHR_LOCUS_FILE		 | The path to the BED file containing DASHR small noncoding RNA loci
TARGETSCAN_DIR			 | The path to the directory containing TargetScan miRNA binding site predictions
F5_CLASSES			 | The text file containing the FANTOM5 tissue class assignments
GTEX_CLASSES			 | The text file containing the GTEx tissue class assignments
ROADMAP_CLASSES			 | The text file containing the Roadmap tissue class assignments
**Parameters for enhancer enrichment analysis** |
NUM_SAMPLES	     | The number of control variant sets to sample (Default = 10,000)
MAF_BIN_SIZE	     | The size of the bins to group variants by minor allele frequency (Default = 0.01)
DIST_ROUND	     | The rounding constant / bin size for distance to the nearest TSS (Default = 1000)
DIST_THRESHOLD	     | An upper threshold on the distance to the nearest TSS, after which all variants are grouped together as 'high' (Default = "Inf" -> no limit)
LD_PARTNER_THRESHOLD	     | An upper threshold on the number of LD partners, after which all variants are grouped together as 'high' (Default = "Inf" -> no limit)
BG_SNP_INFOF		     | The file containing the sampling characteristics for all variants (see below for generation)
LD_SETS_DIR		     | The directory containing the precomputed sets of LD blocks for all variants
REF_SUMMARY_DIR		     | The directory containing the summary files of annotation overlaps for all 1,000 genomes variants
**Parameters for co-localization analysis** |
COLOC_H4_THRESH	 | The threshold on P(H_4) to define a strong colocalization signal (Default = 0.5)
COLOC_ABF_THRESH | The threshold on the amount of cumulative ABF density that should be accounted for by expansion (Default = 0.5)
LOCUSZOOM_PATH	 | If desired, the direct path to the locusZoom executable to generate GWAS and eQTL locusZoom plots for strong colocalization signals
COLOC_GTEX_DIR 	 | The directory containing the full GTEx dataset (not just the significant eQTLs, but all associations)
GTEX_SAMPLE_SIZEF  | The csv containing the GTEx sample sizes for each eQTL dataset, for COLOC
GTEX_RSID_MATCH	   | The file that is used to cross-reference GTEx IDs with rsIDs
HG19_ENSEMBL_REF_FILE	 | The file that is used to cross-reference Ensembl gene IDs with gene names
RELEVANT_CLASSES	 | The set of tissue categories that you especially care about, if any. There must be no spaces in this, or the script breaks! (Default = "'Blood','Brain','Connective Tissue'"). Follow the formatting of the default argument to make sure it gets parsed correctly
**Parameters for lncRNA correlation analysis** |
GTEX_EXPR_DIR	 | The directory containing the parsed RNAseq data (see below)
SAMPLE_INFO_FILE | The file containing the GTEx sample attributes to match IDs with tissues
GENCODE_LNCRNA_FILE    | The file containing the GENCODE lncRNA annotations, used to detect lncRNA eQTL targets
SPEARMAN_THRESH	       | The absolute value threshold on Spearman correlation to define strong lncRNA targets (Default = 0.5)
PEARSON_THRESH	       | The absolute value threshold on Pearson correlation to define strong lncRNA targets (Default = 0.5)
NUM_PCS		       | The number of principal components to regress out for the tissue-specific analysis
**Parameters for metaXcan analysis** |
METAXCAN_DIR |	 The code directory containing the software for MetaXcan (i.e. ~/code/MetaXcan/software/)
GTEX_V7_DBDIR	 |   The directory containing the GTEx v7 databases for PredictDB
**Parameter for LD score regression**	      |
LDSC_CODE_DIR	| The code directory for LD score regression (i.e. ~/code/ldsc)
MUNGE_SNPLIST	| The HapMap3 list of SNPs from MetaXcan (w_hm3.snplist)
LDSC_BASELINE_DIR | The directory containing the extracted baseline annotations from metaXcan
LDSC_WEIGHTS_DIR  | The directory containing the HapMap3 weights for metaXcan
LDSC_FRQ_DIR 	  | The directory containing the 1,000 Genomes frequency information


#### Dependencies and requirements
To run the full pipeline, around 40Gb of memory is required, as the enhancer sampling,
co-localization, and lncRNA correlation analyses are computationally intensive and run in R,
which loads objects into memory. A typical run for the full pipeline will generate around
10-20Gb of result data. The [full annotation
data](http://inferno.lisanwanglab.org/full_INFERNO_annotations.tar.gz) is around
130Gb compressed and expands to around 320Gb of annotation data.

To run the pipeline, the following main programs are required:

* Python v2.7.9 or a higher version of Python 2.7
* bedtools v2.25.0 or greater version of bedtools v2
* R 3.2.3
* plink v1.90b2i or greater

Additionally, the following packages are required:

**For Python:**

* argparse
* subprocess
* datetime
* os
* time
* sys
* glob
* gzip
* re
* subprocess
* copy
* commands
* pickle
* math

**For R:**

* data.table
* ggplot2
* gplots
* gtools
* plyr
* psych
* coloc
* RColorBrewer
* reshape2
* scales

Here is the output of a sessionInfo() call for R package versions used to develop the pipeline:

```R
> sessionInfo()
R version 3.2.3 (2015-12-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] psych_1.7.3.21     coloc_2.3-1        BMA_3.18.6         rrcov_1.4-3       
 [5] inline_0.3.14      robustbase_0.92-7  leaps_3.0          survival_2.41-3   
 [9] MASS_7.3-45        colorspace_1.3-2   gtools_3.5.0       RColorBrewer_1.1-2
[13] gplots_3.0.1       reshape2_1.4.2     scales_0.4.1       plyr_1.8.4        
[17] ggplot2_2.2.1      data.table_1.10.4 

loaded via a namespace (and not attached):
 [1] pcaPP_1.9-61       Rcpp_0.12.10       compiler_3.2.3     DEoptimR_1.0-8    
 [5] bitops_1.0-6       tools_3.2.3        nlme_3.1-131       tibble_1.3.4      
 [9] gtable_0.2.0       lattice_0.20-35    rlang_0.1.2        Matrix_1.2-11     
[13] parallel_3.2.3     mvtnorm_1.0-6      stringr_1.2.0      cluster_2.0.6     
[17] caTools_1.17.1     stats4_3.2.3       grid_3.2.3         foreign_0.8-67    
[21] gdata_2.17.0       magrittr_1.5       splines_3.2.3      mnormt_1.5-5      
[25] KernSmooth_2.23-15 stringi_1.1.5      lazyeval_0.2.0     munsell_0.4.3     
```

#### Running the tool
Then, there are two ways to run the pipeline. The more general approach is to make sure that
python, bedtools, R, and plink are in your path. With this approach, you can run the full
INFERNO pipeline using the INFERNO.py script:

```bash
$ python INFERNO.py -h
usage: INFERNO.py [-h] [--skip_ld_expansion] [--skip_annotation] [--rsid_column RSID_COLUMN] [--pos_column POS_COLUMN] [--pval_column PVAL_COLUMN] [--chr_column CHR_COLUMN]
                  [--allele1_column ALLELE1_COLUMN] [--allele2_column ALLELE2_COLUMN] [--maf_column MAF_COLUMN] [--beta_column BETA_COLUMN] [--run_pval_expansion]
		  [--consistent_direction] [--sig_mult SIG_MULT] [--case_prop CASE_PROP] [--sample_size SAMPLE_SIZE] [--summary_file SUMMARY_FILE] [--run_enhancer_sampling]
		  [--run_gtex_coloc] [--run_lncrna_correlation] [--run_metaXcan] [--summary_has_header] [--run_LDSC] [--run_pathway_analysis] [--cluster_system {bsub,shell}]
		  top_snpf cfg_file outdir outprefix

Driver script for the INFERNO pipeline

positional arguments:
  top_snpf              The tab separated file of the tag SNPs you want to analyze. Should be formatted with four columns: chromosome, rsID, region naming information, and position (in that order). IMPORTANT:
                        Note that SNPs without dbSNP rsIDs should use 'chr-pos' naming format, not 'chr:pos', which is incompatible with this pipeline!
  cfg_file              The configuration file containing paths to all the required functional annotation files. Should be formatted as a bash configuration file i.e. VARIABLE=DEFINTIION on each line.
  outdir                The directory to write all the results to.
  outprefix             The desired prefix for all the output file names

optional arguments:
  -h, --help            show this help message and exit
  --skip_ld_expansion   Give this flag to skip any LD (and p-value-based) expansion and just run analysis directly on the provided list of input variants.
  --skip_annotation     Give this flag to skip all the annotation, enhancer sampling, and co-localization steps. Used if you just want to run MetaXcan or LD score regression with summary stats. You still need
                        to give a top SNP file, but it can be a spoof file since it won't get used for anything.
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
                        The summary statistics column number containing the minor allele frequency. Note that if this is ever greater than 0.5 and a beta column is provided, the effect direction will be
                        flipped to be defined relative to the minor allele.
  --beta_column BETA_COLUMN
                        The summary statistics column number containing the beta estimate (used for p-value expansion with consistent directions). Providing this means that the p-value expansion will consider
                        effect direction given the --consistent_direction flag. This is also required for MetaXcan analysis.
  --run_pval_expansion  If you want to do expansion by p-values when you have summary statistics, provide this flag. Otherwise, the top SNP file will be directly expanded.
  --consistent_direction
                        If you want to do expansion by p-values and also consider effect direction when you have summary statistics, provide this flag. --beta_column is also required for this.
  --sig_mult SIG_MULT   The multiplier range for significance of the p-value expanded variant set (e.g. a value of 10 means one order of magnitude)
  --case_prop CASE_PROP
                        The proportion of cases in the GWAS, for colocalization. If the GWAS is quantitative, set this to 'QUANT' or 'NA'.
  --sample_size SAMPLE_SIZE
                        The total number of samples in the GWAS, for colocalization.
  --summary_file SUMMARY_FILE
                        The path to the full summary statistics file, required for p-value expansion, colocalization analysis, and lncRNA target analysis.
  --run_enhancer_sampling
                        If you want to run the enhancer bootstrapping analysis, provide this flag.
  --run_gtex_coloc      If you want to run COLOC analysis of your summary statistics against GTEx eQTLs from 44 tissues (requires summary statistics)
  --run_lncrna_correlation
                        If you want to analyze expression correlation of any lncRNAs identified by COLOC analysis (--run_gtex_coloc flag) against all other RNAseq-based GTEx genes to find lncRNA targets.
  --run_metaXcan        If you want to run MetaXcan analysis of your summary statistics against GTEx v7 eQTLs from 48 tissues (requires summary statistics)
  --summary_has_header  If you want to run metaXcan or LD score regression, this flag says whether your summary statistics file has a header or not. If not, one will be appended for use with these tools
  --run_LDSC            If you want to run partitioned LD score regression analysis on your summary statistics against the 53 default LD score regression annotations.
  --run_pathway_analysis
                        If you want to run pathway analysis (only applies if you also run lncRNA correlation). Note that the WebGestaltR package must be installed for this, and your environment must have
                        internet access for that package to work, so this will only run if you're directly running the pipeline (i.e. --cluster_system='shell').
  --cluster_system {bsub,shell}
                        If running enhancer sampling, GTEx co-localization, and/or lncRNA correlation, this flag describes how those computationally intensive jobs will be run. The bsub option submits them as
                        separate bsub jobs, while the shell option just runs them sequentially from the same shell as the other INFERNO.py analyses
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

* metaXcan_GTEx_v7/: The results of the MetaXcan analysis on GTEx v7 data containing the
  MetaXcan csv file outputs.
    * metaxcan_analysis/: plots and tables summarizing the metaXcan results.

* LD_score_regression/: The results of the LD score regression output, which follows the
  default organization provided by the tool. 

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

Then, to generate annotation overlaps and summaries for all of these variants, all variants from 1,000 genomes were run through the pipeline:

```bash
$ cd /path/to/INFERNO/src/
$ python expand_and_annotate_snps.py --loglevel save --kg_pop EUR --skip_ld_expansion --unstranded_partition_dir \
   ~/data/refgenomes/hg19/unstranded_partitions/utr_annotations/final_files/ \
   --fantom5_dir ~/data/FANTOM5/Enhancers/facet_expressed_enhancers/sorted/ --enhancer_locus_window 1000 --skip_closest_enh --skip_enh_summary --fantom5_correlation_file \
   ~/data/FANTOM5/Enhancers/enhancer_tss_associations.bed --gtex_dir ~/data/GTEx/single_cell_sig_eqtls_v6/sorted/ --factorbook_file \
   ~/data/factorbook/wgEncodeRegTfbsClusteredWithCellsV3.sorted.bed \
   --roadmap_chromhmm_dir ~/data/roadmap/chromHMM/sorted/ --homer_motif_bed_file ~/data/HOMER_motifs/homer.sorted.KnownMotifs.hg19.bed --homer_motif_pwm_file \
   ~/data/HOMER_motifs/custom.motifs \
   --homer_motif_seq_file ~/data/HOMER_motifs/homer.sorted.KnownMotifs.hg19.sequence.txt ~/data/1000_genomes/phase1_release_v3/EUR/sorted_files/ \
   ~/data/enhancer_snp_pipeline/input_data/no_indels_1kg_EUR_phase1_v3_snps.txt ~/data/enhancer_snp_pipeline/output/all_1kg_EUR_phase1_v3_snps_11_14_16/ all_1kg_EUR
$ cd ../analysis_scripts/
$ time ./large_scale_count_annotation_overlaps.sh -l 1000 -f ~/data/FANTOM5/Enhancers/fantom5_classes.txt -g ~/data/GTEx/gtex_classes.txt -r ~/data/roadmap/roadmap_classes.txt \
  ~/data/enhancer_snp_pipeline/output/all_1kg_EUR_phase1_v3_snps_11_14_16/ 1.0 0 all_1kg_EUR ~/data/enhancer_snp_pipeline/output/all_1kg_EUR_phase1_v3_snps_11_14_16/summaries/
$ cd ~/data/enhancer_snp_pipeline/output/all_1kg_EUR_phase1_v3_snps_11_14_16/summaries/
$ time cut -f1-3,6 enh_locus_snps_1.0_ld_0_dist.txt | sort -u > uniq_class_enh_locus_snps.txt; \
  time cut -f1-3,6 eqtl_snps_1.0_ld_0_dist.txt | sort -u > uniq_class_eqtl_snps.txt; \
  time cut -f1-3,6 roadmap_hmm_snps_1.0_ld_0_dist.txt | sort -u > uniq_class_roadmap_hmm_snps.txt
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
from
http://gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz.

Then, a few pre-processing scripts were used to parse these for use in the correlation script:

```bash
$ cd /path/to/INFERNO/data_preprocessing/
$ python find_ensembl_genes_by_chr.py /path/to/hg19_reference/ensembl_hg19_genes.txt /path/to/RNAseq_output_dir/
$ ./all_chr_gtex_expression_filter.sh /path/to/RNAseq_output_dir/ /path/to/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz
```

To identify lncRNA eQTL targets, the GENCODE lncRNA annotations were
[downloaded](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz)
and extracted from
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz.

Finally, the sample sizes for each GTEx tissue were also
[downloaded](http://www.gtexportal.org/home/tissueSummaryPage#sampleCountsPerTissue) as a csv
from http://www.gtexportal.org/home/tissueSummaryPage#sampleCountsPerTissue.

#### MetaXcan implementation:
INFERNO implements [MetaXcan](https://github.com/hakyimlab/MetaXcan) analysis using the
MetaMany script and GTEx v7 databases from the [PredictDB site](http://predictdb.org/). The
following steps were taken to download and process these annotations:

```bash
$ cd /path/to/code/
$ git clone https://github.com/hakyimlab/MetaXcan
$ cd /path/to/PredictDB/
$ wget https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29.tar.gz
$ wget https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29_README.txt
$ mkdir GTEx_v7_dbs
$ time tar -C GTEx_v7_dbs/ -xzvf GTEx-V7_HapMap-2017-11-29.tar.gz
```

Then, the parameters in the config file are set as follows:

METAXCAN_DIR="/path/to/code/MetaXcan/software/"

GTEX_V7_DBDIR="/path/to/PredictDB/GTEx_v7_dbs/"

Once this is set up, INFERNO takes care of running the MetaXcan analysis using the
"--run_metaXcan" flag.

#### LD score regression implementation:
INFERNO implements partitioned heritability [LD score
regression](https://github.com/bulik/ldsc) analysis using the baseline model LD scores,
regression weights, and allele frequencies, as described in the [partitioned heritability
wiki](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability). The following steps were
taken to download and process these annotations:

```bash
$ cd /path/to/code/
$ git clone https://github.com/bulik/ldsc.git
$ module load anaconda/3 ## on a module-based server; otherwise, just have anaconda in your $PATH
$ cd ldsc/
$ conda env create --file environment.yml
$ cd /path/to/ldsc_annotations/
$ wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_baseline_ldscores.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_frq.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_cell_type_groups.tgz
$ wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
$ bzip2 -d w_hm3.snplist.bz2
$ time for f in *tgz; do tar -xzf $f; done
```

Then, the parameters in the config file are set as follows:

LDSC_CODE_DIR="/path/to/code/ldsc/"

MUNGE_SNPLIST="/path/to/ldsc_annotations/w_hm3.snplist"

LDSC_BASELINE_DIR="/path/to/ldsc_annotations/baseline/"

LDSC_WEIGHTS_DIR="/path/to/ldsc_annotations/weights_hm3_no_hla/"

LDSC_FRQ_DIR="/path/to/ldsc_annotations/1000G_frq/"

Once these are set up, INFERNO will run LD score regression using the "--run_LDSC" flag.