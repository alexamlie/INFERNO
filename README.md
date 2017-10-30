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
	
### Other instructions:
To perform pathway analysis using WebGestalt on any list of target genes, use the 'genesymbol'
option.

### Data sources and pre-processing steps:
