## gtex_gwas_colocalization_analysis.R
## alex amlie-wolf 04/20/17
## using the coloc package to check for colocalization between any given GWAS and GTEx eQTL signals
## this script requires a lot of memory for parsing the full GTEx eQTL data in step 4 (even
## though i use system calls, it still has to read stuff into memory) step 4 takes a really
## long time to run, while steps 5 and 6 are meant to be run (separately) after that
## it also runs in around 3 days, typically..

library(coloc)
library(ggplot2)
library(reshape2)
library(plyr)
sessionInfo()

## for locuszoom analysis, set the path to the binary:
locuszoom <- "/appl/locuszoom-1.3/bin/locuszoom"

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------
## 0. Table of Contents
## 1. Function Definitions
## 2. Read in data, set data directories
## 3. Run COLOC analysis across all top GWAS loci and GTEx data
## 4. High-level analysis of GTEx COLOC results with top GWAS hits
## 5. Compare GTEx COLOC hits with top GWAS enhancer overlaps and eQTL effect direction
## 6. Compare GTEx COLOC hits with motif disruptions
## 7. Compare GTEx COLOC hits with direct eQTL overlaps

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------

## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {    
make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
    if(type=='pdf') {
        pdf(file=paste0(filename, ".pdf"), width=10*width_ratio, height=10*height_ratio, pointsize=12, onefile=FALSE)
    } else if(type=='png') {
        ## use type='cairo' for when X11 doesn't work
        png(filename=paste0(filename, ".png"), width=10*width_ratio, height=10*height_ratio, res=300, units='in', type='cairo')
    } else {
        cat('filetype not supported\n')
    }
}

## a convenience function from INFERNO to check if a parameter is in the reference vector. only
## returns true if the parameter is not set to None or False (from Python), i.e. the only time
## we want a true return value is if the parameter is set for real
check_param <- function(param_vec, param) {
    ret_val <- FALSE
    if (param %in% names(param_vec)) {
        if(param_vec[[param]] != "None" & param_vec[[param]] != "False") {
            ret_val <- TRUE
        }
    }
    ret_val
}

## -----------------------------------------------------------------------------
## 2. Read in data, set data directories
## -----------------------------------------------------------------------------
## set a variable for printing messages or not
PRINT_MSGS <- FALSE

outdir <- "/home/alexaml/data/ad_enhancer_analysis/gtex_gwas_colocalization_analysis/"
dir.create(paste0(outdir, "/plots/"), F, T)
dir.create(paste0(outdir, "/tables/"), F, T)

## read in GWAS data:
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==4) {
    ## set up the INFERNO directory for enhancer and direct eQTL overlaps
    inferno_param_f <- args[1]
    ## define the P(H_4) threshold for finding the strongest hits
    coloc_thresh <- args[2]
    ## get the top GWAS SNPs (INFERNO input_)
}


## read in the parameters and make a named vector for reference
parameter_tab <- read.table(inferno_param_f, header=F, sep="\t", quote="", as.is=T, col.names=c("param", "value"))
param_ref <- parameter_tab$value
names(param_ref) <- parameter_tab$param
if(check_param(param_ref, 'skip_ld_expansion')) {
    r2_thresh <- "1.0"
    dist_thresh <- "0"
} else {
    r2_thresh <- param_ref[['ld_threshold']]
    dist_thresh <- param_ref[['ld_check_area']]
}
## define a prefix for the output files that use the INFERNO data
outprefix <- param_ref[['outprefix']]

## read in the SNPs that were analyzed by INFERNO, to cross reference
ld_stats_file <- paste0(param_ref[['outdir']], '/ld_expansion/', param_ref[['outprefix']],
                        "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, ".txt")
ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)


## read in the list of the top IGAP SNPs
top_igap_snpf <- "/home/alexaml/data/IGAP_data/IGAP_phase1_top_hits_no_HLA_DSG2.txt"
top_igap_snps <- read.table(top_igap_snpf, header=F, sep="\t", as.is=T,
                            col.names=c("chr", "rsID", "region", "pos"))

## set up the IGAP directory (we don't read all the data in at once, just parse it for the
## specific locus we're looking at)
igap_summary_dir <- "/home/alexaml/data/IGAP_data/"

## same for the GTEx directory, since we read in specific tissues
gtex_dir <- "/home/alexaml/data/GTEx/GTEx_v6p/GTEx_Analysis_v6p_all-associations/"

## read in the sample sizes
gtex_sample_sizes <- read.csv("/home/alexaml/data/GTEx/GTEx_v6p/GTEx_sample_sizes.csv")
colnames(gtex_sample_sizes) <- c("Tissue", "eqtl_samp_num", "rnaseq_samp_num", "egene_num")
gtex_sample_sizes$tissmatch <- gsub("- ", "", gsub("\\(|\\)", "", gtex_sample_sizes$Tissue))

## also read in the GTEx tissue classes
gtex_class_file <- '/home/alexaml/data/GTEx/gtex_classes.txt'
gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
## to match up with the GTEx samples
gtex_category_df$coloc_match <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data, fixed=T)

## define the file that is used to match GTEx IDs with rsIDs
gtex_rsid_match <- "/home/alexaml/data/GTEx/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt"

## also set up the eGWAS data directory
egwas_dir <- "/home/alexaml/data/mayo_egwas/"

## define the tissue classes that we care about here
relevant_classes <- c("Blood", "Brain", "Connective Tissue")

## finally, read in the reference to get gene names
hg19_ensembl_ref_file <- "/home/alexaml/data/refgenomes/hg19/hg19_ensembl_to_name.txt"
ensembl_xref <- read.table(hg19_ensembl_ref_file, header=T, sep="\t", quote="", comment.char="", as.is=T, col.names=c("tx_id", "gene_id", "gene_name"))

## -----------------------------------------------------------------------------
## 3. Run COLOC analysis across all top GWAS loci and GTEx data
## -----------------------------------------------------------------------------
{
coloc_start_time <- proc.time()
## save the GTEx eQTL column names
gtex_eqtl_header <- c("gene_id", "variant_id", "tss_distance", "pval_nominal", "slope", "slope_se")

## go through all the chromosomes since we read in the IGAP data by chromosome
for(this_chr in unique(top_igap_snps$chr)) {
    cat("Parsing", this_chr, "\n")

    ## read in the IGAP data
    this_igap_data <- read.table(paste0(igap_summary_dir, "/rs_tbl_sorted_", this_chr, ".txt"), header=F, sep="\t", quote="", as.is=T)
    ## set the header
    colnames(this_igap_data) <- c("rsID", "pos", "Marker", "Allele1", "Allele2", "Freq1",
                                  "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr",
                                  "pval", "Direction", "HetISq", "HetChiSq", "df", "HetPVal")

    ## read in the GTEx rsID matches for this chromosome    
    this_chr_gtex_rsid_match <- read.table(pipe(paste0("awk -F$'\t' '$1==\"Chr\" || $1==",
                                                       gsub("chr", "", this_chr),
                                                       "' ", gtex_rsid_match)),
                                                header=T, sep="\t", quote="", as.is=T)
    
    ## now analyze each tag region
    for(this_tag in unique(top_igap_snps$region[top_igap_snps$chr==this_chr])) {
        cat("Analyzing", this_tag, "region\n")
        dir.create(paste0(outdir, '/tables/coloc_analysis/gtex_coloc/top_regions/', gsub("/", "_", this_tag, fixed=T)), F, T)

        this_tag_data <- top_igap_snps[top_igap_snps$region==this_tag,]
        
        this_tag_igap_data <- this_igap_data[this_igap_data$rsID==this_tag_data$rsID,]
        
        if(nrow(this_tag_igap_data) != 1) {
            cat("Could not find single match for tag snp", this_tag_data$rsID, "in IGAP data\n")
            next
        }

        this_region_igap_snps <- this_igap_data[abs(this_igap_data$pos - this_tag_igap_data$pos) <= 500000,]
        cat(nrow(this_region_igap_snps), "SNPs found in IGAP around", this_tag_data$rsID, "\n")

        ## add a column for the variant IDs in the GTEx format
        ## ({chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37)
        ## the IGAP alleles are defined as the ones that were first identified, so we have to try both
        ## orderings to get as many overlaps as we possibly can..
        this_region_igap_snps$gtex_id1 <-
            ## allele 1 is minor
            paste(gsub("chr", "", this_chr), this_region_igap_snps$pos,
                  toupper(this_region_igap_snps$Allele2),
                  toupper(this_region_igap_snps$Allele1), "b37", sep="_")
        
        this_region_igap_snps$gtex_id2 <-
            ## allele 2 is minor
            paste(gsub("chr", "", this_chr), this_region_igap_snps$pos,
                  toupper(this_region_igap_snps$Allele1),
                  toupper(this_region_igap_snps$Allele2), "b37", sep="_") 

        ## to match up with GTEx ID, we need to try to manually figure out if our tag SNP is
        ## oriented correctly or not. note that this uses IGAP frequencies and it might not
        ## actually match up with the GTEx annotations (i deal with this in the gtex_file loop)
        if(this_tag_igap_data$Freq1 > 0.50) {
            ## allele 1 is major (so we need to flip)
            this_tag_gtex_id <-
                paste(gsub("chr", "", this_chr), this_tag_igap_data$pos,
                      toupper(this_tag_igap_data$Allele1),
                      toupper(this_tag_igap_data$Allele2), "b37", sep="_") 
        } else {
            ## allele 1 is minor
            this_tag_gtex_id <-
                paste(gsub("chr", "", this_chr), this_tag_igap_data$pos,
                      toupper(this_tag_igap_data$Allele2),
                      toupper(this_tag_igap_data$Allele1), "b37", sep="_") 
        }

        ## now we have to loop through each GTEx tissue and define the gene sets to perform colocalization analysis on
        for(gtex_file in list.files(gtex_dir, "*.txt.gz")) {
            gtex_tiss_match <- gsub("_", " ", gsub("_Analysis", "", strsplit(gtex_file, ".", fixed=T)[[1]][1]))
            cat("Analyzing tissue", gtex_tiss_match, "\n")

            ## try to find all the genes tested with this tag SNP
            all_genes <- system2("zcat", paste0(gtex_dir, gtex_file, " | awk -F$'\t' '$2==\"", this_tag_gtex_id, "\"' | cut -f1 | sort -u"), stdout=TRUE)
                        
            ## check that we actually found some genes, and flip the alleles if not
            if(length(all_genes)==0) {
                ## we only need to change the ID we use to search here because the direction of
                ## effect is taken care of when we search the IDs of the full IGAP dataset
                cat("IGAP minor allele for tag SNP", this_tag_igap_data$rsID, "does not match up with GTEx ID (no genes found using IGAP minor allele); flipping direction and re-testing\n")
                
                ## if we thought allele 1 was major, it's actually minor
                if(this_tag_igap_data$Freq1 > 0.50) {
                    this_tag_gtex_id <-
                        paste(gsub("chr", "", this_chr), this_tag_igap_data$pos,
                              toupper(this_tag_igap_data$Allele2),
                              toupper(this_tag_igap_data$Allele1), "b37", sep="_") 
                } else {
                    ## otherwise, we thought allele 1 was minor but it's actually major
                    this_tag_gtex_id <-
                        paste(gsub("chr", "", this_chr), this_tag_igap_data$pos,
                              toupper(this_tag_igap_data$Allele1),
                              toupper(this_tag_igap_data$Allele2), "b37", sep="_") 
                }
                                
                all_genes <- system2("zcat", paste0(gtex_dir, gtex_file, " | awk -F$'\t' '$2==\"", this_tag_gtex_id, "\"' | cut -f1 | sort -u"), stdout=TRUE)

                ## if this is still 0, just continue
                if(length(all_genes)==0) {
                    cat("No genes found for", this_tag_igap_data$rsID, this_tag, "in tissue",
                        strsplit(gtex_file, ".", fixed=T)[[1]][1], "!\n")
                    next
                }
            }

            ## write the gene list to a file for awk purposes
            gene_outf <- paste0(outdir, '/tables/coloc_analysis/gtex_coloc/top_regions/', gsub("/", "_", this_tag, fixed=T), "/", gsub("_Analysis", "", strsplit(gtex_file, ".", fixed=T)[[1]][1]), "_gene_ids.txt")
            write.table(all_genes, gene_outf, quote=F, sep="\t", row.names=F, col.names=F)
            
            ## read in all the data for all these genes at once
            all_eqtl_data <- data.frame(do.call(rbind, strsplit(
                system2("zcat", paste0(gtex_dir, gtex_file,
                                       " | awk -F$'\t' 'NR==FNR {gene_ids[$1]; next} {if($1 in gene_ids) {print $0}}' ",
                                       gene_outf, " - "), stdout=TRUE),
                split="\t")), stringsAsFactors=F)
            colnames(all_eqtl_data) <- gtex_eqtl_header

            cat("Found", length(unique(all_eqtl_data$gene_id)), "genes tested in tissue",
                strsplit(gtex_file, ".", fixed=T)[[1]][1], "taking up",
                format(object.size(all_eqtl_data), units="auto"), "with", nrow(all_eqtl_data),
                "total eQTL pairs tested\n")

            ## only make the directory if we found any genes to test
            this_tiss_outdir <- paste0(outdir, '/tables/coloc_analysis/gtex_coloc/top_regions/', gsub("/", "_", this_tag, fixed=T), "/", gsub(" ", "_", gtex_tiss_match), "/")
            dir.create(this_tiss_outdir, F, T)
            
            ## finally, loop through all these genes and perform colocalization analysis on each
            for(this_gene_id in all_genes) {
                gene_name <- unique(ensembl_xref$gene_name[grep(strsplit(this_gene_id, ".", fixed=T)[[1]][1], ensembl_xref$gene_id)])
                if(length(gene_name) > 1) {
                    cat("Multiple gene names found for", this_gene_id, ":", gene_name, ".. using the first one found\n")
                    gene_name <- gene_name[1]
                }

                this_eqtl_data <- all_eqtl_data[all_eqtl_data$gene_id==this_gene_id,]

                cat("Found", nrow(this_eqtl_data), "SNPs tested for eQTL with", gene_name, this_gene_id, "\n")

                ## reset some of the columns to have the correct type for the coloc package
                this_eqtl_data$pval_nominal <- as.numeric(this_eqtl_data$pval_nominal)
                this_eqtl_data$slope <- as.numeric(this_eqtl_data$slope)
                this_eqtl_data$slope_se <- as.numeric(this_eqtl_data$slope_se)
                
                ## parse these datasets to only include overlapping SNPs
                gtex_id1_match <- this_region_igap_snps$gtex_id1 %in% this_eqtl_data$variant_id
                gtex_id2_match <- this_region_igap_snps$gtex_id2 %in% this_eqtl_data$variant_id
                this_region_igap_snps.parsed <- this_region_igap_snps[gtex_id1_match | gtex_id2_match,]
                
                ## set a new column for the actual matching
                ## we use this matching a lot so save it 
                gtex_id1_parsed_match <- this_region_igap_snps.parsed$gtex_id1 %in% this_eqtl_data$variant_id
                this_region_igap_snps.parsed$gtex_id <- ifelse(gtex_id1_parsed_match, this_region_igap_snps.parsed$gtex_id1, this_region_igap_snps.parsed$gtex_id2)
                ## we also need to set the frequency to always be the minor allele frequency
                ## gtex_id1 corresponds to allele 1 being minor
                this_region_igap_snps.parsed$MAF <- ifelse(gtex_id1_parsed_match, this_region_igap_snps.parsed$Freq1, 1-this_region_igap_snps.parsed$Freq1)
                
                ## similarly, we need to flip the effect directions when we have flipped alleles (the IGAP
                ## effect is calculated relative to allele 1, and the GTEx effects are calculated relative to
                ## the alternate allele, so this needs to be consistent)
                this_region_igap_snps.parsed$correct_beta <- ifelse(gtex_id1_parsed_match, this_region_igap_snps.parsed$Effect, -this_region_igap_snps.parsed$Effect)
                
                ## also parse the eQTL data
                this_eqtl_data.parsed <- this_eqtl_data[this_eqtl_data$variant_id %in% this_region_igap_snps.parsed$gtex_id,]
                
                ## check that the number of SNPs are equal
                cat("Analyzing", nrow(this_region_igap_snps.parsed), "SNPs from IGAP and", nrow(this_eqtl_data.parsed), "matching SNPs tested for eQTL\n")
                if(nrow(this_region_igap_snps.parsed)!=nrow(this_eqtl_data.parsed)) {
                    cat("Unequal number of SNPs in region", this_tag, "and gene", gene_name, this_gene_id, "\n")
                    next
                }
                
                ## check that everything is in the same order
                if(sum(this_eqtl_data.parsed$variant_id != this_region_igap_snps.parsed$gtex_id) != 0) {
                    cat("Wrong order in region", this_tag, "and gene", gene_name, this_gene_id, "\n")
                    next
                }

                this_gtex_samp_size <- gtex_sample_sizes$eqtl_samp_num[gtex_sample_sizes$tissmatch==gtex_tiss_match]
                
                ## need to use the squared standard error for the variance argument
                this_beta_coloc_results <- coloc.abf(
                    dataset1=list(beta=this_region_igap_snps.parsed$correct_beta,
                        varbeta=this_region_igap_snps.parsed$StdErr^2,
                        type="cc", s=0.3140209, N=54162),
                    dataset2=list(beta=this_eqtl_data.parsed$slope,
                        varbeta=this_eqtl_data.parsed$slope_se^2,
                        type="quant", N=this_gtex_samp_size),
                    MAF=this_region_igap_snps.parsed$MAF)

                write.table(this_beta_coloc_results[['summary']],
                            file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_summary.txt"),
                            quote=F, sep="\t", col.names=F)

                ## before writing out the full results, need to add a column so that we can
                ## actually match up the SNPs
                ## define the indexing vector
                coloc_snp_idx <- as.numeric(unlist(lapply(strsplit(this_beta_coloc_results[['results']]$snp, ".", fixed=T), "[[", 2)))
                
                this_beta_coloc_results[['results']]$rsID <- as.character(this_region_igap_snps.parsed$rsID[coloc_snp_idx])
                
                write.table(this_beta_coloc_results[['results']],
                            file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_full_results.txt"),
                            quote=F, sep="\t", col.names=T, row.names=F)

                ## if this comparison meets the threshold, make locuszoom plots:
                if (this_beta_coloc_results[['summary']]['PP.H4.abf'] >= coloc_thresh) {
                    ## annotate the GTEx data with the matching rsIDs
                    this_eqtl_data.parsed$dbsnp135_rsid <- this_chr_gtex_rsid_match$RS_ID_dbSNP135_original_VCF[match(this_eqtl_data.parsed$variant_id, this_chr_gtex_rsid_match$VariantID)]
                    ## now write out that info
                    write.table(this_eqtl_data.parsed[,c("dbsnp135_rsid", "pval_nominal")],
                                paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_gtex_pvals.txt"),
                                quote=F, sep="\t", col.names=T, row.names=F)

                    ## now make a folder within the plots for this tag region
                    this_tag_lz_out <- paste0(outdir, '/plots/locuszoom_plots/', gsub("/", "_", this_tag, fixed=T), "_region/")
                    dir.create(this_tag_lz_out, F, T)

                    ## do the IGAP analysis, only if we didn't already (this only needs to be run once, even with different GTEx genes and tissues)
                    if(!file.exists(paste0(this_tag_lz_out, gsub("/", "_", this_tag, fixed=T), "_region_igap_locuszoom_", this_tag_data$rsID, ".pdf"))) {
                        write.table(this_region_igap_snps.parsed[,c("rsID", "pval")],
                                    paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_igap_pvals.txt"),
                                    quote=F, sep="\t", col.names=T, row.names=F)
                        
                        ## make an IGAP locuszoom plot around the tag SNP
                        system2(locuszoom, paste("--metal",
                                                 paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_igap_pvals.txt"),
                                                 "--markercol rsID --pvalcol pval --refsnp",
                                                 this_tag_data$rsID, "--flank 500kb --pop EUR --build hg19 --source 1000G_Nov2010 --plotonly --prefix",
                                                 paste0(this_tag_lz_out, gsub("/", "_", this_tag, fixed=T), "_region_igap_locuszoom"),
                                                 "--no-date legend='right'"))
                    }
                    
                    ## same thing for the eQTL data
                    ## first find the top ABF SNP so that we can label that
                    top_abf_snp <- this_beta_coloc_results[['results']]$rsID[which.max(this_beta_coloc_results[['results']]$SNP.PP.H4)]
                    max_abf <- max(this_beta_coloc_results[['results']]$SNP.PP.H4)

                    ## if it's the same, just label the ABF of the tag variant
                    if(top_abf_snp == this_tag_data$rsID) {
                        write.table(data.frame(snp=this_tag_data$rsID, string=paste0("GWAS tag (", round(max_abf, digits=2), ")"), color="purple"),
                                    file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_lz_labels.txt"),
                                    quote=F, sep="\t", row.names=F)
                        
                    } else {
                        ## find the tag variant ABF
                        tag_abf <- this_beta_coloc_results[['results']]$SNP.PP.H4[this_beta_coloc_results[['results']]$rsID==this_tag_data$rsID]
                        
                        ## make the label file
                        write.table(data.frame(snp=c(this_tag_data$rsID, top_abf_snp),
                                               string=c(paste0("GWAS tag (", round(tag_abf, digits=2), ")"), paste0("Top ABF (", round(max_abf, digits=2), ")")), color=c("purple", "blue")),
                                    file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_lz_labels.txt"),
                                    quote=F, sep="\t", row.names=F)                        
                    }

                    ## now run the command
                    system2(locuszoom, paste("--metal",
                                             paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_gtex_pvals.txt"),
                                             "--markercol dbsnp135_rsid --pvalcol pval_nominal --refsnp",
                                             this_tag_data$rsID, "--flank 500kb --pop EUR --build hg19 --source 1000G_Nov2010 --plotonly --prefix",
                                             paste0(this_tag_lz_out, gsub(" ", "_", gtex_tiss_match), "_", gene_name, "_gtex_locuszoom"),
                                             paste0("title='", gene_name, " eQTL in ", gtex_tiss_match,"'"),
                                             "--denote-markers-file",
                                             paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_lz_labels.txt"),                                                 
                                             "--no-date legend='right'"))
                    
                }

                ## don't save this one, just keep it for comparison
                cat("P-value-based colocalization, for comparison\n")
                this_pval_coloc_results <- coloc.abf(
                    dataset1=list(pvalues=this_region_igap_snps.parsed$pval,
                        type="cc", s=0.3140209, N=54162),
                    dataset2=list(pvalues=this_eqtl_data.parsed$pval_nominal,
                        type="quant", N=this_gtex_samp_size),
                    MAF=this_region_igap_snps.parsed$MAF)
                
                rm(this_eqtl_data)
            }
        }
        
    }
    
}
cat("Colocalization analysis took", (proc.time() - coloc_start_time)[["elapsed"]],
    "seconds\n")
}

## -----------------------------------------------------------------------------
## 4. High-level analysis of GTEx COLOC results with top GWAS hits
## -----------------------------------------------------------------------------
{
cur_tag <- ""
summ_start <- proc.time()

## i want to read in all the summary data into one unified data frame
all_summ_files <- list.files(paste0(outdir, '/tables/coloc_analysis/gtex_coloc/top_regions/'), pattern="*.summary.txt",
                             full.names=T, recursive=T)

all_summary_data <- data.frame(stringsAsFactors = F)
for(f in all_summ_files) {
    ## split the path to get the tag region, tissue, and target gene
    this_info <- strsplit(gsub(outdir, "", gsub("//", "/", f)), "/")[[1]]    
    
    this_tag <- this_info[5]
    if(this_tag != cur_tag) {
        cur_tag <- this_tag
        cat("Analyzing", cur_tag, "region\n")
    }
    this_tiss <- this_info[6]

    this_eqtl_info <- strsplit(this_info[7], "_")[[1]]
    this_gene_name <- this_eqtl_info[1]
    this_gene_id <- this_eqtl_info[2]

    ## now actually read in the summary data
    this_data <- read.table(f, header=F, sep="\t", quote="", as.is=T, nrows=6)

    all_summary_data <- rbind(all_summary_data, data.frame(tag_region=this_tag, tissue=this_tiss,
                                                           eqtl_gene_name=this_gene_name,
                                                           eqtl_gene_id=this_gene_id,
                                                           nsnps=this_data[1,2],
                                                           PP.H0.abf=this_data[2,2],
                                                           PP.H1.abf=this_data[3,2],
                                                           PP.H2.abf=this_data[4,2],
                                                           PP.H3.abf=this_data[5,2],
                                                           PP.H4.abf=this_data[6,2]))
}
cat("Reading in all summaries took", (proc.time() - summ_start)[["elapsed"]], "seconds\n")
}

## add a tissue class category
all_summary_data$gtex_tissue_class <- gtex_category_df$Class[match(all_summary_data$tissue, gtex_category_df$coloc_match)]
## and reorder
all_summary_data <- all_summary_data[,c(1,2,11,3:10)]

## write this out in summarized form
write.table(all_summary_data, paste0(outdir, '/tables/coloc_analysis/IGAP_top_hits_gtex_coloc_summaries.txt'), quote=F, sep="\t", row.names=F)

## ## to read this in directly:
## all_summary_data <- read.table(paste0(outdir, '/tables/coloc_analysis/IGAP_top_hits_gtex_coloc_summaries.txt'), header=T, sep="\t", quote="")

## make histograms for each of the different hypotheses by tag region
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H0_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H0.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H0 (no causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H0 (no causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H1_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H1.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H1 (GWAS causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H1 (GWAS causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H2_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H2.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H2 (eQTL causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H2 (eQTL causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H3_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H3.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H3 (unshared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H3 (unshared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H4_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H4_highest_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0.5, 1), breaks=seq(0.5, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

## also do this split by tissue
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H0_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H0.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H0 (no causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H0 (no causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H1_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H1.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H1 (GWAS causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H1 (GWAS causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H2_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H2.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H2 (eQTL causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H2 (eQTL causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H3_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H3.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H3 (unshared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H3 (unshared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H4_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_H4_highest_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) +
    scale_x_continuous(limits=c(0.5, 1), breaks=seq(0.5, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

## now make a single histogram looking at the distributions of all hypotheses together
## need to melt the data to do this
melted_summary_data <- melt(all_summary_data, id.vars=1:6, variable.name="hypothesis", value.name="probability")

make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_all_hypotheses_prob_density'))
ggplot(melted_summary_data, aes(x=probability, fill=hypothesis, color=hypothesis)) +
    geom_density(alpha=0.5) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    scale_fill_brewer(palette="Set1") +
    scale_colour_brewer(palette="Set1") + 
    ggtitle("Combined density plot of all hypotheses, GTEx") + 
    theme_bw() + xlab("Probability of hypothesis") + ylab("Density") +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15),
          title=element_text(size=25), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5))
dev.off()

## -----------------------------------------------------------------------------
## 5. Compare GTEx COLOC hits with top GWAS enhancer overlaps and eQTL effect direction
## -----------------------------------------------------------------------------
## use the colocalization threshold parameter
top_coloc_hits <- all_summary_data[all_summary_data$PP.H4.abf > coloc_thresh,]
nrow(top_coloc_hits)

## define the probability threshold to grab the set of SNPs up to
coloc_prob_thresh <- 0.5
## coloc_prob_thresh <- 0.9
## coloc_prob_thresh <- 0.99

{
enh_start_time <- proc.time()

## read in the FANTOM5 and roadmap enhancer hits
fantom5_overlap_f <- paste0(param_ref[['outdir']], '/analysis_results/fantom5_overlap/tables/', param_ref[['outprefix']], "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, "_locus_", param_ref[['enhancer_locus_window']], "bp_window_uniq_fantom5_overlap_snps.txt")
fantom5_enh_overlaps <- read.table(fantom5_overlap_f, header=T, sep="\t", quote="", as.is=T)

roadmap_overlap_f <- paste0(param_ref[['outdir']], '/analysis_results/roadmap_chromHMM/tables/', param_ref[['outprefix']], "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, "_uniq_chromHMM_enh_snps.txt")
roadmap_enh_overlaps <- read.table(roadmap_overlap_f, header=T, sep="\t", quote="", as.is=T)

## now for each of the most colocalized hits, figure out if the top SNP is in one of these
## datasets.  we want to store all this information: which SNP, which tag region, which eQTL
## dataset and target gene, etc
top_coloc_enh_overlaps <- data.frame(stringsAsFactors = F)
## also make one to store more variants, accounting for some amount (coloc_prob_thresh) of the
## individual SNP probability
top_coloc_enh_overlaps.expanded <- data.frame(stringsAsFactors = F)

for(i in seq(nrow(top_coloc_hits))) {
    ## this starts out as just a row of the coloc data but is supplemented
    this_comparison <- top_coloc_hits[i,]
    ## save the original row, for the expanded SNP set analysis
    this_orig_comparison <- this_comparison

    ## to read in the full data, we need to reconstruct the path
    this_data_file <- paste0(outdir, '/tables/coloc_analysis/gtex_coloc/top_regions/', this_comparison$tag_region, '/',
                             this_comparison$tissue, '/', this_comparison$eqtl_gene_name, '_',
                             this_comparison$eqtl_gene_id, '_full_results.txt')

    this_coloc_data <- read.table(this_data_file, header=T, sep="\t", quote="", as.is=T)
    
    ## get the best SNP for the summary data frame. start by ordering all the probabilities
    indiv_snp_prob_order <- order(this_coloc_data$SNP.PP.H4, decreasing=T)

    top_causal_snp_data <- this_coloc_data[indiv_snp_prob_order[1],]
    top_causal_snp_rsid <- top_causal_snp_data$rsID
    eqtl_z_score <- top_causal_snp_data$z.df2
    eqtl_variance <- top_causal_snp_data$V.df2
    
    ## store this SNP, its ABF, and eQTL statistics including Z score, variance, and beta
    this_comparison <- cbind(this_comparison, top_coloc_snp=top_causal_snp_rsid,
                             max_coloc_prob=max(this_coloc_data$SNP.PP.H4),
                             eqtl_z_score=eqtl_z_score,
                             eqtl_variance=eqtl_variance,
                             eqtl_beta=eqtl_z_score * sqrt(eqtl_variance))

    ## also get the set of SNPs responsible for at least 50% of the probability
    ## to do this, calculate the cumulative sum of the ordered probabilities
    snp_cumsum_values <- cumsum(this_coloc_data$SNP.PP.H4[indiv_snp_prob_order])
    ## get the SNP set, enough to get 0.5 probability
    snp_set_idx <- indiv_snp_prob_order[1:(sum(snp_cumsum_values <= coloc_prob_thresh)+1)]
    
    top_causal_snp_set <- this_coloc_data$rsID[snp_set_idx]
    expanded_snp_data <- this_coloc_data[snp_set_idx,] 
    
    ## -----------------------------
    ## first do the analysis against the single top SNP
    top_causal_f5_hits <- fantom5_enh_overlaps[fantom5_enh_overlaps$rsID==top_causal_snp_rsid,]
    top_causal_hmm_hits <- roadmap_enh_overlaps[roadmap_enh_overlaps$rsID==top_causal_snp_rsid,]

    if(nrow(top_causal_f5_hits) > 0) {
        if(PRINT_MSGS) {
            cat("Shared causal SNP", top_causal_snp_rsid, "for", as.character(this_comparison$eqtl_gene_name), "eQTL in GTEx",
                as.character(this_comparison$tissue), "in", as.character(this_comparison$tag_region),
                "region found to overlap FANTOM5 enhancer\n")
        }
        ## store information about the fantom5 overlaps
        all_f5_tissues <- unique(top_causal_f5_hits$enh_source)
        all_f5_classes <- unique(top_causal_f5_hits$enh_class)

        ## look for class overlap
        f5_consistent_classes <- "None"
        if(any(all_f5_classes==this_comparison$gtex_tissue_class)) {
            if(PRINT_MSGS) {            
                cat("GTEx tissue category", this_comparison$gtex_tissue_class, "is consistent with FANTOM5 overlap\n")
            }
        }
        
        this_comparison <- cbind(this_comparison, num_f5_tissues=length(all_f5_tissues), f5_tissues=paste(all_f5_tissues, collapse=","), num_f5_tissue_classes=length(all_f5_classes), f5_tissue_classes=paste(all_f5_classes, collapse=","))
    } else {
        ## fill in with negative data
        this_comparison <- cbind(this_comparison, num_f5_tissues=0, f5_tissues="None", num_f5_tissue_classes=0, f5_tissue_classes="None")        
    }

    if(nrow(top_causal_hmm_hits)==1) {
        if(PRINT_MSGS) {
            cat("Shared causal SNP", top_causal_snp_rsid, "for", as.character(this_comparison$eqtl_gene_name), "eQTL in",
                as.character(this_comparison$tissue), "in", as.character(this_comparison$tag_region),
                "region found to overlap Roadmap enhancer\n")
        }

        ## look for class overlap
        all_hmm_classes <- strsplit(top_causal_hmm_hits$hmm_classes, ",")[[1]]
        if(any(all_hmm_classes==this_comparison$gtex_tissue_class)) {
            if(PRINT_MSGS) {
                cat("GTEx tissue category", this_comparison$gtex_tissue_class, "is consistent with Roadmap overlap\n")
            }
        }

        ## store information about the roadmap overlaps (this one is easier because it's
        ## already computed)        
        this_comparison <- cbind(this_comparison, num_hmm_tissues=top_causal_hmm_hits$num_hmm_tissues, hmm_tissues=top_causal_hmm_hits$hmm_tissues, num_hmm_enh_states=top_causal_hmm_hits$num_hmm_enh_states, hmm_enh_states=top_causal_hmm_hits$hmm_enh_states, num_hmm_classes=top_causal_hmm_hits$num_hmm_classes, hmm_classes=top_causal_hmm_hits$hmm_classes)
    } else if(nrow(top_causal_hmm_hits) > 1) {
        ## check this just for debugging
        cat("Too many Roadmap entries for SNP", top_causal_snp_rsid, "\n")
    } else {
        this_comparison <- cbind(this_comparison, num_hmm_tissues=0, hmm_tissues="None", num_hmm_enh_states=0, hmm_enh_states="None", num_hmm_classes=0, hmm_classes="None")        
    }

    top_coloc_enh_overlaps <- rbind(top_coloc_enh_overlaps, this_comparison)

    ## -----------------------------    
    ## now go through the full set of overlapping SNPs and check for enhancer hits    
    ## store all of these SNPs with their ABFs in a new data structure (which gets modified
    ## later and then incorporated row by row into top_coloc_enh_overlaps.expanded)
    this_expanded_comparison <- cbind(this_orig_comparison[rep(1, length(top_causal_snp_set)),],
                                      num_prob_expanded_snps=length(top_causal_snp_set),
                                      high_coloc_snp=top_causal_snp_set,
                                      coloc_snp_prob=this_coloc_data$SNP.PP.H4[snp_set_idx],
                                      eqtl_z_score=expanded_snp_data$z.df2,
                                      eqtl_variance=expanded_snp_data$V.df2,
                                      eqtl_beta=expanded_snp_data$z.df2*sqrt(expanded_snp_data$V.df2))

    for(j in seq(length(top_causal_snp_set))) {
        this_rsid <- top_causal_snp_set[j]
        ## grab the row from the first expanded data frame
        this_snp_comparison <- this_expanded_comparison[j,]

        this_snp_causal_f5_hits <- fantom5_enh_overlaps[fantom5_enh_overlaps$rsID==this_rsid,]
        this_snp_causal_hmm_hits <- roadmap_enh_overlaps[roadmap_enh_overlaps$rsID==this_rsid,]

        if(nrow(this_snp_causal_f5_hits) > 0) {
            if(PRINT_MSGS) {
                cat("Shared causal SNP", this_rsid, "for", as.character(this_comparison$eqtl_gene_name), "eQTL in GTEx",
                    as.character(this_comparison$tissue), "in", as.character(this_comparison$tag_region),
                    "region found to overlap FANTOM5 enhancer\n")
            }
            
            ## store information about the fantom5 overlaps
            all_f5_tissues <- unique(this_snp_causal_f5_hits$enh_source)
            all_f5_classes <- unique(this_snp_causal_f5_hits$enh_class)

            this_snp_comparison <- cbind(this_snp_comparison, num_f5_tissues=length(all_f5_tissues), f5_tissues=paste(all_f5_tissues, collapse=","), num_f5_tissue_classes=length(all_f5_classes), f5_tissue_classes=paste(all_f5_classes, collapse=","))
        } else {
            ## fill in with negative data
            this_snp_comparison <- cbind(this_snp_comparison, num_f5_tissues=0, f5_tissues="None", num_f5_tissue_classes=0, f5_tissue_classes="None")        
        }

        if(nrow(this_snp_causal_hmm_hits)==1) {
            if(PRINT_MSGS) {
                cat("Shared causal SNP", this_rsid, "for", as.character(this_comparison$eqtl_gene_name), "eQTL in",
                    as.character(this_comparison$tissue), "in", as.character(this_comparison$tag_region),
                    "region found to overlap Roadmap enhancer\n")
            }
            
            ## store information about the roadmap overlaps (this one is easier because it's
            ## already computed)        
            this_snp_comparison <- cbind(this_snp_comparison, num_hmm_tissues=this_snp_causal_hmm_hits$num_hmm_tissues, hmm_tissues=this_snp_causal_hmm_hits$hmm_tissues, num_hmm_enh_states=this_snp_causal_hmm_hits$num_hmm_enh_states, hmm_enh_states=this_snp_causal_hmm_hits$hmm_enh_states, num_hmm_classes=this_snp_causal_hmm_hits$num_hmm_classes, hmm_classes=this_snp_causal_hmm_hits$hmm_classes)
        } else if(nrow(this_snp_causal_hmm_hits) > 1) {
            ## check this just for debugging
            cat("Too many Roadmap entries for SNP", this_snp_causal_snp_rsid, "\n")
        } else {
            this_snp_comparison <- cbind(this_snp_comparison, num_hmm_tissues=0, hmm_tissues="None", num_hmm_enh_states=0, hmm_enh_states="None", num_hmm_classes=0, hmm_classes="None")        
        }

        top_coloc_enh_overlaps.expanded <- rbind(top_coloc_enh_overlaps.expanded, this_snp_comparison)        
    }
}

## add a column for whether there is enhancer overlap or not
top_coloc_enh_overlaps$any_enh_overlap <- ifelse(top_coloc_enh_overlaps$num_f5_tissues > 0 | top_coloc_enh_overlaps$num_hmm_tissues > 0, "yes", "no")
## also add a column for whether the GTEx tissue category is found in roadmap or fantom5
top_coloc_enh_overlaps$f5_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps$gtex_tissue_class, x=top_coloc_enh_overlaps$f5_tissue_classes), "yes", "no")
top_coloc_enh_overlaps$hmm_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps$gtex_tissue_class, x=top_coloc_enh_overlaps$hmm_classes), "yes", "no")
## finally, add a column to describe the three levels: enhancer overlap without tissue match, with tissue match, and no enhancer overlap
top_coloc_enh_overlaps$enh_levels <- with(top_coloc_enh_overlaps, ifelse(any_enh_overlap=="no", "no_overlap", ifelse(f5_category_match=="yes" | hmm_category_match=="yes", "tiss_match", "tiss_non_match")))

write.table(top_coloc_enh_overlaps, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_overlaps.txt'), quote=F, sep="\t", row.names=F)

## do the same thing for the expanded data
## add a column for whether there is enhancer overlap or not
top_coloc_enh_overlaps.expanded$any_enh_overlap <- ifelse(top_coloc_enh_overlaps.expanded$num_f5_tissues > 0 | top_coloc_enh_overlaps.expanded$num_hmm_tissues > 0, "yes", "no")
## also add a column for whether the GTEx tissue category is found in roadmap or fantom5
top_coloc_enh_overlaps.expanded$f5_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps.expanded$gtex_tissue_class, x=top_coloc_enh_overlaps.expanded$f5_tissue_classes), "yes", "no")
top_coloc_enh_overlaps.expanded$hmm_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps.expanded$gtex_tissue_class, x=top_coloc_enh_overlaps.expanded$hmm_classes), "yes", "no")
## use these new columns to describe another convenient column for summarizing the tissue matches
top_coloc_enh_overlaps.expanded$enh_levels <- with(top_coloc_enh_overlaps.expanded, ifelse(any_enh_overlap=="no", "no_overlap", ifelse(f5_category_match=="yes" | hmm_category_match=="yes", "tiss_match", "tiss_non_match")))

write.table(top_coloc_enh_overlaps.expanded, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_overlaps.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

cat("Enhancer overlap analysis took", (proc.time() - enh_start_time)[["elapsed"]], "seconds\n")
}

## check how many match up
cat(sum(unique(top_coloc_enh_overlaps$top_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps$top_coloc_snp)),
    "unique max colocalized SNPs match up with INFERNO SNPs\n")

cat(sum(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp)),
    "unique high colocalized SNPs in", coloc_prob_thresh, "expanded set match up with INFERNO SNPs\n")

## -----------------------------
## do some visualizations of the single max coloc SNP overlaps
## start with histograms of the maximum colocalization probabilities
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(top_coloc_enh_overlaps, aes(x=max_coloc_prob, fill=tag_region)) +
    geom_histogram(binwidth=0.01) +
    scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of maximum SNP causal probability across tag regions, GTEx") + 
    theme_bw() + xlab("Highest probability of SNP being causal") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

## now look at the distributions split by annotation overlap
## set the factor levels to get the plots in the right order
top_coloc_enh_overlaps$enh_levels <- factor(top_coloc_enh_overlaps$enh_levels, levels=c("no_overlap", "tiss_non_match", "tiss_match"), ordered=T)

## look at the individual SNP probability distributions by enhancer overlap and tag region
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_overlap_and_tag_region'), height_ratio=1.5, width_ratio=1.5)
ggplot(top_coloc_enh_overlaps, aes(x=enh_levels, y=max_coloc_prob, color=enh_levels)) +
    geom_boxplot() +
    facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
    ggtitle("Distributions of maximum SNP causal probability, GTEx") + 
    theme_bw() + xlab("Type of enhancer overlap") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20))
dev.off()

## same thing, but collapsed across tag regions
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_overlap'))
ggplot(top_coloc_enh_overlaps, aes(x=enh_levels, y=max_coloc_prob, color=enh_levels)) +
    geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
    ggtitle("Distributions of maximum SNP causal probability across tag regions, GTEx") + 
    theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20))
dev.off()

## next, summarize the counts of highly colocalized SNPs in each region, whether they overlap an enhancer, and whether the category matches or not
tag_region_enh_count_summary <- ddply(top_coloc_enh_overlaps, .(tag_region), summarize,
                                      num_coloc_hits = length(tag_region),
                                      enh_overlap_hits = sum(any_enh_overlap=="yes"),
                                      enh_non_overlap_hits = sum(any_enh_overlap=="no"),
                                      enh_tiss_match = sum(f5_category_match=="yes" | hmm_category_match=="yes"),
                                      enh_tiss_non_match = sum(any_enh_overlap=="yes" & f5_category_match=="no" & hmm_category_match=="no"))

## add in 0 counts for the regions we did not observe
all_regions <- unique(all_summary_data$tag_region)
tag_region_enh_count_summary <- rbind(tag_region_enh_count_summary,
                                      data.frame(tag_region=all_regions[!(all_regions %in% tag_region_enh_count_summary$tag_region)], num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))

## now melt this down into the three categories we care about
melted_enh_count_summ <- melt(tag_region_enh_count_summary, id.vars="tag_region", measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

melted_enh_count_summ$variable <- factor(melted_enh_count_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_overlaps_barplot'))
ggplot(melted_enh_count_summ, aes(x=tag_region, y=count, fill=variable)) +
    scale_fill_manual(
        labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
        values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle("Enhancer analysis of colocalized SNPs") + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=3, title="")) + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5))
dev.off()

## eQTL effect direction analysis for the unexpanded set
## number of positive and negative hits
table(ifelse(top_coloc_enh_overlaps$eqtl_beta > 0, "positive", "negative"))
## negative positive 
##       85       69 

## now plot beta distributions across tag regions
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_beta_dists_across_tag_regions'))
print(ggplot(top_coloc_enh_overlaps, aes(x=tag_region, y=eqtl_beta, color=eqtl_beta, fill=tag_region)) + 
xlab("Tag region") + ylab("eQTL beta value") + theme_bw() +
ggtitle("COLOC eQTL beta distributions across tag regions for top ABF SNPs") + 
geom_point() + geom_violin(alpha=0.5) + geom_hline(color="red", yintercept=0, linetype=3) + 
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot Z-score distributions across tag regions
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_zscore_dists_across_tag_regions'))
print(ggplot(top_coloc_enh_overlaps, aes(x=tag_region, y=eqtl_z_score, color=eqtl_z_score, fill=tag_region)) + 
xlab("Tag region") + ylab("eQTL Z score") + theme_bw() +
ggtitle("COLOC eQTL Z score distributions across tag regions for top ABF SNPs") + 
geom_point() + geom_violin(alpha=0.5) + geom_hline(color="red", yintercept=0, linetype=3) + 
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## variance distributions across tag regions
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_variance_dists_across_tag_regions'))
print(ggplot(top_coloc_enh_overlaps, aes(x=tag_region, y=eqtl_variance, color=eqtl_variance, fill=tag_region)) + 
xlab("Tag region") + ylab("eQTL variance") + theme_bw() +
ggtitle("COLOC eQTL variance distributions across tag regions for top ABF SNPs") + 
geom_point() + geom_violin(alpha=0.5) + 
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## -----------------------------
## next look at the expanded dataset
## first we want histograms of the number of SNPs for comparison by tag region
max_snp_set_size <- max(top_coloc_enh_overlaps.expanded$num_prob_expanded_snps)

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_num_snps_to_', coloc_prob_thresh, '_thresh_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
ggplot(unique(top_coloc_enh_overlaps.expanded[,c("tag_region", "tissue", "eqtl_gene_name", "num_prob_expanded_snps")]), aes(x=num_prob_expanded_snps, fill=tag_region)) +
    geom_histogram(binwidth=1) +
    scale_x_continuous(breaks=seq(0, max_snp_set_size, by=ifelse(max_snp_set_size > 50, 10, 2))) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle(paste("Histograms of numbers of SNPs to reach", coloc_prob_thresh, "probability")) + 
    theme_bw() + xlab("Number of SNPs to reach probability threshold") + ylab("Number of GTEx - GWAS comparisons") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5))
dev.off()

## set the factor levels to get the plots in the right order
top_coloc_enh_overlaps.expanded$enh_levels <- factor(top_coloc_enh_overlaps.expanded$enh_levels, levels=c("no_overlap", "tiss_non_match", "tiss_match"), ordered=T)

## look at the individual SNP probability distributions by enhancer overlap and tag region
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_overlap_and_tag_region_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
ggplot(top_coloc_enh_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
    geom_boxplot() +
    facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
    ggtitle(paste("Distributions of SNP ABFs, expanded to", coloc_prob_thresh, "probability, GTEx")) + 
    theme_bw() + xlab("Type of enhancer overlap") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20))
dev.off()

## same thing, but collapsed across tag regions
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_overlap_', coloc_prob_thresh, '_prob_thresh'))
ggplot(top_coloc_enh_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
    geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
    ggtitle(paste("Distributions of SNP ABFs across tag regions\nexpanded to", coloc_prob_thresh, "probability, GTEx")) + 
    theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20))
dev.off()

## next, summarize the counts of highly colocalized SNPs in each region, whether they overlap an enhancer, and whether the category matches or not
expanded_tag_region_enh_count_summary <- ddply(top_coloc_enh_overlaps.expanded, .(tag_region), summarize,
                                               num_uniq_snp_comparisons = length(tag_region),
                                               enh_overlap_hits = sum(any_enh_overlap=="yes"),
                                               enh_non_overlap_hits = sum(any_enh_overlap=="no"),
                                               enh_tiss_match = sum(f5_category_match=="yes" | hmm_category_match=="yes"),
                                               enh_tiss_non_match = sum(any_enh_overlap=="yes" & f5_category_match=="no" & hmm_category_match=="no"))

## add in 0 counts for the regions we did not observe
expanded_tag_region_enh_count_summary <- rbind(expanded_tag_region_enh_count_summary,
                                               data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_enh_count_summary$tag_region)], num_uniq_snp_comparisons=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))

## now melt this down into the three categories we care about
melted_expanded_enh_count_summ <- melt(expanded_tag_region_enh_count_summary, id.vars="tag_region", measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

melted_expanded_enh_count_summ$variable <- factor(melted_expanded_enh_count_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_overlaps_barplot_', coloc_prob_thresh, '_prob_thresh'))
ggplot(melted_expanded_enh_count_summ, aes(x=tag_region, y=count, fill=variable)) +
    scale_fill_manual(
        labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
        values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer overlaps of colocalized SNPs expanded to", coloc_prob_thresh, "probability, GTEx")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=3, title="")) + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5))
dev.off()

## eQTL effect direction analysis for the expanded set
## check whether the eQTLs in each comparison line up or not:
expanded_eqtl_info <- ddply(top_coloc_enh_overlaps.expanded, .(tag_region, tissue, gtex_tissue_class, eqtl_gene_name), summarize, positive_betas = sum(eqtl_beta > 0, na.rm=T), negative_betas = sum(eqtl_beta < 0, na.rm=T))

cat(sum(with(expanded_eqtl_info, positive_betas > 0 & negative_betas > 0)),
    "of", nrow(expanded_eqtl_info), "ABF-expanded sets have variants with different beta directions\n")
cat("Of the rest,", sum(with(expanded_eqtl_info, positive_betas > 0 & negative_betas==0)), "are all positive and",
    sum(with(expanded_eqtl_info, positive_betas==0 & negative_betas > 0)), "are all negative\n")

## now plot beta distributions across comparisons, for the full set of 154 comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_beta_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_beta, color=eqtl_beta)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL beta values") + theme_bw() +
ggtitle(paste("COLOC eQTL beta distributions, expanded to", coloc_prob_thresh, "probability")) + 
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot z-score distributions across comparisons, for the full set of 154 comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_zscore_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_z_score, color=eqtl_z_score)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL Z scores") + theme_bw() +
ggtitle(paste("COLOC eQTL Z score distributions, expanded to", coloc_prob_thresh, "probability")) + 
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot variance distributions across comparisons, for the full set of 154 comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_eqtl_variance_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_variance, color=eqtl_variance)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL variance") + theme_bw() +
ggtitle(paste("COLOC eQTL variance distributions, expanded to", coloc_prob_thresh, "probability")) + 
geom_point() + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## make the same plots for just the relevant tissue categories
## now plot beta distributions across comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_relevant_class_eqtl_beta_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_beta, color=eqtl_beta)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL beta values") + theme_bw() +
ggtitle(paste("COLOC eQTL beta distributions, expanded to", coloc_prob_thresh, "probability\n",
              "Only comparisons in", paste(relevant_classes, collapse=", "))) + 
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot z-score distributions across comparisons, for the full set of 154 comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_relevant_class_eqtl_zscore_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_z_score, color=eqtl_z_score)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL Z scores") + theme_bw() +
ggtitle(paste("COLOC eQTL Z score distributions, expanded to", coloc_prob_thresh, "probability\n",
              "Only comparisons in", paste(relevant_classes, collapse=", "))) +
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot variance distributions across comparisons
make_graphic(paste0(outdir, '/plots/IGAP_top_hits_gtex_coloc_relevant_class_eqtl_variance_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_variance, color=eqtl_variance)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL variance") + theme_bw() +
ggtitle(paste("COLOC eQTL variance distributions, expanded to", coloc_prob_thresh, "probability\n",
              "Only comparisons in", paste(relevant_classes, collapse=", "))) +
geom_point() + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## -----------------------------------------------------------------------------
## 6. Compare GTEx COLOC hits with motif disruptions
## -----------------------------------------------------------------------------
## require PWM disruption calculations
if(check_param(param_ref, "homer_motif_pwm_file") & check_param(param_ref, "homer_motif_seq_file")) {
    ## read in the overlap data
    homer_pwm_f <- paste0(param_ref[['outdir']], '/analysis_results/homer_motif_overlap/tables/', param_ref[['outprefix']], "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, "_homer_motif_overlaps.txt")

    homer_pwm_overlaps <- read.table(homer_pwm_f, header=T, sep="\t", quote="", as.is=T)

    ## -----------------------------------    
    ## start with the single top colocalized SNP data
    ## -----------------------------------
    ## we want to add this information in summarized form to the enhancer overlap table but
    ## also have a file with all the unique motif overlaps in separate rows, so make two data
    ## frames
    top_coloc_enh_motif_overlaps <- data.frame(stringsAsFactors = F)
    ## set up an empty DF for the motif overlaps only
    top_coloc_motif_overlaps <- data.frame(stringsAsFactors = F)
    ## create an indexing vector for the non-enhancer columns of the DF
    non_enh_cols <- colnames(top_coloc_enh_overlaps) %in% colnames(top_coloc_hits) |
        colnames(top_coloc_enh_overlaps) %in% c("top_coloc_snp", "max_coloc_prob")
    
    ## now go through and find pwm overlaps
    for(i in seq(nrow(top_coloc_enh_overlaps))) {
        this_data <- top_coloc_enh_overlaps[i,]

        ## find all the motif overlaps for the top SNP here
        this_motif_overlaps <- homer_pwm_overlaps[homer_pwm_overlaps$rsID==this_data$top_coloc_snp,]
        ## only include it in the output if there are any hits
        if(nrow(this_motif_overlaps) > 0) {
            ## first add all these separate files to the motif-only DF
            top_coloc_motif_overlaps <- rbind(top_coloc_motif_overlaps, cbind(
                this_data[rep(1, times=nrow(this_motif_overlaps)),non_enh_cols],
                this_motif_overlaps))
            
            ## now make a summary and add it to the DF with enhancer info as well
            this_data <- cbind(this_data,
                               num_motif_hits=nrow(this_motif_overlaps), 
                               motif_chr=unique(this_motif_overlaps$chr), 
                               motif_starts=paste(this_motif_overlaps$motif_start, collapse=";"),
                               motif_ends=paste(this_motif_overlaps$motif_end, collapse=";"),
                               tfs=paste(this_motif_overlaps$tf_name, collapse=";"),
                               log_odds_scores=paste(this_motif_overlaps$log_odds_score, collapse=";"),
                               strands=paste(this_motif_overlaps$strand, collapse=";"),
                               relative_snp_pos=paste(this_motif_overlaps$relative_snp_pos, collapse=";"),
                               ref_pwms=paste(this_motif_overlaps$ref_pwm, collapse=";"),
                               alt_pwms=paste(this_motif_overlaps$alt_pwm, collapse=";"),
                               delta_pwms=paste(this_motif_overlaps$delta_pwm, collapse=";"),
                               rounded_delta_pwms=paste(round(this_motif_overlaps$delta_pwm, digits=2), collapse=";"), 
                               ref_probs=paste(this_motif_overlaps$ref_prob, collapse=";"),
                               alt_probs=paste(this_motif_overlaps$alt_prob, collapse=";"),
                               motif_seqs=paste(this_motif_overlaps$motif_seq, collapse=";"))
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_overlaps <- rbind(top_coloc_enh_motif_overlaps, this_data)
        } else {
            ## if there are no overlaps, note that in the summary
            ## now make a summary and add it to the DF with enhancer info as well
            this_data <- cbind(this_data, num_motif_hits=0, 
                               motif_chr="NA", motif_starts="NA", motif_ends="NA", tfs="NA",
                               log_odds_scores="NA", strands="NA", relative_snp_pos="NA",
                               ref_pwms="NA", alt_pwms="NA", delta_pwms="NA",
                               rounded_delta_pwms="NA", ref_probs="NA", alt_probs="NA",
                               motif_seqs="NA")
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_overlaps <- rbind(top_coloc_enh_motif_overlaps, this_data)
        }
        
    }

    ## add a column to easily describe the now 6 different levels: the 3 enhancer levels X
    ## having any motif hits or not
    top_coloc_enh_motif_overlaps$motif_levels <- with(top_coloc_enh_motif_overlaps,
        ## if we have a motif hit, deal with that
        ifelse(num_motif_hits > 0, ifelse(enh_levels=="no_overlap", "motif_no_enh",
                                   ifelse(enh_levels=="tiss_non_match", "motif_non_tiss_match", "motif_tiss_match")),
        ## otherwise, no motif hit
        ifelse(enh_levels=="no_overlap", "no_motif_no_enh",
               ifelse(enh_levels=="tiss_non_match", "no_motif_non_tiss_match", "no_motif_tiss_match"))))
    
    ## write these out
    write.table(top_coloc_enh_motif_overlaps, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_motif_overlap_summary.txt'), quote=F, sep="\t", row.names=F)

    write.table(top_coloc_motif_overlaps, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_motif_overlaps.txt'), quote=F, sep="\t", row.names=F)
    
    ## now look at the individual SNP probability distributions
    ## set the factor levels to get the plots in the right order
    top_coloc_enh_motif_overlaps$motif_levels <- factor(top_coloc_enh_motif_overlaps$motif_levels, levels=c("no_motif_no_enh", "no_motif_non_tiss_match", "no_motif_tiss_match", "motif_no_enh", "motif_non_tiss_match", "motif_tiss_match"), ordered=T)
    ## make another factor so that i can facet on motif ovelrap
    top_coloc_enh_motif_overlaps$motif_overlap <- factor(ifelse(top_coloc_enh_motif_overlaps$num_motif_hits > 0, "TFBS", "no TFBS"))
    
    ## look at the individual SNP probability distributions by enhancer overlap and tag region
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_motif_overlap_and_tag_region'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_overlaps, aes(x=motif_levels, y=max_coloc_prob, color=motif_levels)) +
        geom_boxplot() +
        facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
        ggtitle("Distributions of maximum SNP causal probability, GTEx") + 
        theme_bw() + xlab("Type of motif and enhancer overlap") + ylab("Individual SNP ABFs") +
        theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, size=15),
              axis.text.y = element_text(size=15), strip.text=element_text(size=25),
              legend.text = element_text(size=15), legend.title=element_text(size=25),
              axis.title=element_text(size=25),
              plot.title=element_text(hjust=0.5, size=20)))
    dev.off()

    ## same thing, but collapsed across tag regions
    ## make two plots for this: one with motif facets, and one straight one
    ## this one is not faceted
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_motif_overlap'))
    print(ggplot(top_coloc_enh_motif_overlaps, aes(x=motif_levels, y=max_coloc_prob, color=motif_levels)) +
        geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
        ggtitle("Distributions of maximum SNP causal probability across tag regions, GTEx") + 
        theme_bw() + xlab("Type of motif and enhancer overlaps") + ylab("Individual SNP ABFs") +
        theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
              axis.text.y = element_text(size=15), strip.text=element_text(size=25),
              legend.text = element_text(size=15), legend.title=element_text(size=25),
              axis.title=element_text(size=25),
              plot.title=element_text(hjust=0.5, size=20)))
    dev.off()    
    
    ## this one has facets
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_motif_overlap_motif_facet'))
    print(ggplot(top_coloc_enh_motif_overlaps, aes(x=enh_levels, y=max_coloc_prob, color=enh_levels)) +
        geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
        ggtitle("Distributions of maximum SNP causal probability across tag regions, GTEx") +
        scale_color_manual(values=c("tiss_match"="#B3112E", "tiss_non_match"="darkred", "no_overlap"="gray20")) +
        scale_x_discrete(breaks=c("no_overlap", "tiss_non_match", "tiss_match"), labels=c("No enhancer overlap", "Inconsistent tissue class", "Consistent tissue class")) +
        facet_grid(motif_overlap ~ ., scales="free") +           
        theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
        theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
              axis.text.y = element_text(size=15), strip.text=element_text(size=25),
              legend.text = element_text(size=15), legend.title=element_text(size=25),
              axis.title=element_text(size=25),
              plot.title=element_text(hjust=0.5, size=20)))
    dev.off()    
    
    ## now make a summary barplot split by tag region and whether there is a motif hit or not
    tag_region_motif_enh_summary <- ddply(top_coloc_enh_motif_overlaps, .(tag_region), function(x) {
        motif_hits <- x[x$num_motif_hits > 0,]
        non_motif_hits <- x[x$num_motif_hits==0,]

        ## summarize the characteristics of the motif hits
        motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(motif_hits),
            enh_overlap_hits = sum(motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(motif_hits$any_enh_overlap=="no"),
            enh_tiss_match = sum(motif_hits$f5_category_match=="yes" |
                motif_hits$hmm_category_match=="yes"),
            enh_tiss_non_match = sum(motif_hits$any_enh_overlap=="yes" &
                motif_hits$f5_category_match=="no" & motif_hits$hmm_category_match=="no"))
        
        ## summarize the characteristics of the non-motif hits
        non_motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(non_motif_hits),
            enh_overlap_hits = sum(non_motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(non_motif_hits$any_enh_overlap=="no"),
            enh_tiss_match = sum(non_motif_hits$f5_category_match=="yes" |
                non_motif_hits$hmm_category_match=="yes"),
            enh_tiss_non_match = sum(non_motif_hits$any_enh_overlap=="yes" &
                non_motif_hits$f5_category_match=="no" & non_motif_hits$hmm_category_match=="no"))

        return(rbind(cbind(hit_motif="motif_overlap", motif_hit_summary),
                     cbind(hit_motif="no_motif_overlap", non_motif_hit_summary)))
    })

    ## add in 0 counts for the regions we did not observe
    tag_region_motif_enh_summary <- rbind(tag_region_motif_enh_summary,
                                          data.frame(tag_region=all_regions[!(all_regions %in% tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0),
                                          data.frame(tag_region=all_regions[!(all_regions %in% tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))
    
    ## write this summary table out
    write.table(tag_region_motif_enh_summary, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_motif_enh_tag_region_summary.txt'), quote=F, sep="\t", row.names=F)
    
    ## melt this down
    melted_enh_motif_summ <- melt(tag_region_motif_enh_summary, id.vars=c("tag_region", "hit_motif"), measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

    melted_enh_motif_summ$variable <- factor(melted_enh_motif_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_and_motif_overlaps_barplot'), height_ratio=1.5)
    print(ggplot(melted_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle("Enhancer and motif analysis of GTEx-colocalized SNPs") + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=3, title="")) +
    facet_wrap(~hit_motif, nrow=2, scales="free") +     
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
    dev.off()
    
    ## -----------------------------------    
    ## next do the expanded SNP set
    ## -----------------------------------
    ## we want to add this information in summarized form to the enhancer overlap table but
    ## also have a file with all the unique motif overlaps in separate rows, so make two data
    ## frames
    top_coloc_enh_motif_overlaps.expanded <- data.frame(stringsAsFactors = F)
    ## set up an empty DF for the motif overlaps only
    top_coloc_motif_overlaps.expanded <- data.frame(stringsAsFactors = F)
    ## create an indexing vector for the non-enhancer columns of the DF
    non_enh_cols.expanded <- colnames(top_coloc_enh_overlaps.expanded) %in% colnames(top_coloc_hits) |
        colnames(top_coloc_enh_overlaps.expanded) %in% c("num_prob_expanded_snps", "high_coloc_snp", "coloc_snp_prob")
    
    ## now go through and find pwm overlaps
    for(i in seq(nrow(top_coloc_enh_overlaps.expanded))) {
        this_data <- top_coloc_enh_overlaps.expanded[i,]

        ## find all the motif overlaps for this particular SNP
        this_motif_overlaps <- homer_pwm_overlaps[homer_pwm_overlaps$rsID==this_data$high_coloc_snp,]
        ## only include it in the output if there are any hits
        if(nrow(this_motif_overlaps) > 0) {
            ## first add all these separate files to the motif-only DF
            top_coloc_motif_overlaps.expanded <- rbind(top_coloc_motif_overlaps.expanded, cbind(
                this_data[rep(1, times=nrow(this_motif_overlaps)),non_enh_cols.expanded],
                this_motif_overlaps))
            
            ## now make a summary and add it to the DF with enhancer info as well
            this_data <- cbind(this_data,
                               num_motif_hits=nrow(this_motif_overlaps), 
                               motif_chr=unique(this_motif_overlaps$chr), 
                               motif_starts=paste(this_motif_overlaps$motif_start, collapse=";"),
                               motif_ends=paste(this_motif_overlaps$motif_end, collapse=";"),
                               tfs=paste(this_motif_overlaps$tf_name, collapse=";"),
                               log_odds_scores=paste(this_motif_overlaps$log_odds_score, collapse=";"),
                               strands=paste(this_motif_overlaps$strand, collapse=";"),
                               relative_snp_pos=paste(this_motif_overlaps$relative_snp_pos, collapse=";"),
                               ref_pwms=paste(this_motif_overlaps$ref_pwm, collapse=";"),
                               alt_pwms=paste(this_motif_overlaps$alt_pwm, collapse=";"),
                               delta_pwms=paste(this_motif_overlaps$delta_pwm, collapse=";"),
                               rounded_delta_pwms=paste(round(this_motif_overlaps$delta_pwm, digits=2), collapse=";"),                                
                               ref_probs=paste(this_motif_overlaps$ref_prob, collapse=";"),
                               alt_probs=paste(this_motif_overlaps$alt_prob, collapse=";"),
                               motif_seqs=paste(this_motif_overlaps$motif_seq, collapse=";"))
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_overlaps.expanded <- rbind(top_coloc_enh_motif_overlaps.expanded, this_data)
        } else {
            ## if there are no overlaps, note that in the summary
            ## now make a summary and add it to the DF with enhancer info as well
            this_data <- cbind(this_data, num_motif_hits=0, 
                               motif_chr="NA", motif_starts="NA", motif_ends="NA", tfs="NA",
                               log_odds_scores="NA", strands="NA", relative_snp_pos="NA",
                               ref_pwms="NA", alt_pwms="NA", delta_pwms="NA",
                               rounded_delta_pwms="NA", ref_probs="NA", alt_probs="NA",
                               motif_seqs="NA")
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_overlaps.expanded <- rbind(top_coloc_enh_motif_overlaps.expanded, this_data)
        }
        
    }

    ## add a column to easily describe the now 6 different levels: the 3 enhancer levels X
    ## having any motif hits or not
    top_coloc_enh_motif_overlaps.expanded$motif_levels <- with(top_coloc_enh_motif_overlaps.expanded,
        ## if we have a motif hit, deal with that
        ifelse(num_motif_hits > 0, ifelse(enh_levels=="no_overlap", "motif_no_enh",
                                   ifelse(enh_levels=="tiss_non_match", "motif_non_tiss_match", "motif_tiss_match")),
        ## otherwise, no motif hit
        ifelse(enh_levels=="no_overlap", "no_motif_no_enh",
               ifelse(enh_levels=="tiss_non_match", "no_motif_non_tiss_match", "no_motif_tiss_match"))))

    ## write these out
    write.table(top_coloc_enh_motif_overlaps.expanded, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_motif_overlap_summary.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    write.table(top_coloc_motif_overlaps.expanded, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_motif_overlaps.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## now look at the individual SNP probability distributions
    ## set the factor levels to get the plots in the right order
    top_coloc_enh_motif_overlaps.expanded$motif_levels <- factor(top_coloc_enh_motif_overlaps.expanded$motif_levels, levels=c("no_motif_no_enh", "no_motif_non_tiss_match", "no_motif_tiss_match", "motif_no_enh", "motif_non_tiss_match", "motif_tiss_match"), ordered=T)
    ## make another factor so that i can facet on motif ovelrap
    top_coloc_enh_motif_overlaps.expanded$motif_overlap <- factor(ifelse(top_coloc_enh_motif_overlaps.expanded$num_motif_hits > 0, "TFBS", "no TFBS"))

    ## look at the individual SNP probability distributions by enhancer overlap and tag region
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_motif_overlap_and_tag_region_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_overlaps.expanded, aes(x=motif_levels, y=coloc_snp_prob, color=motif_levels)) +
        geom_boxplot() +
        facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
        ggtitle(paste("Distributions of SNP ABFs, expanded to", coloc_prob_thresh, "probability, GTEx")) + 
        theme_bw() + xlab("Type of motif and enhancer overlap") + ylab("Individual SNP ABFs") +
        theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, size=15),
              axis.text.y = element_text(size=15), strip.text=element_text(size=25),
              legend.text = element_text(size=15), legend.title=element_text(size=25),
              axis.title=element_text(size=25),
              plot.title=element_text(hjust=0.5, size=20)))
    dev.off()

    ## same thing, but collapsed across tag regions
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_motif_overlap_', coloc_prob_thresh, '_prob_thresh'))
    print(ggplot(top_coloc_enh_motif_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
        geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
        ggtitle(paste("Distributions of SNP ABFs across tag regions\nexpanded to", coloc_prob_thresh, "probability, GTEx")) +
        scale_color_manual(values=c("tiss_match"="#B3112E", "tiss_non_match"="darkred", "no_overlap"="gray20")) +
        scale_x_discrete(breaks=c("no_overlap", "tiss_non_match", "tiss_match"), labels=c("No enhancer overlap", "Inconsistent tissue class", "Consistent tissue class")) + 
        facet_grid(motif_overlap ~ ., scales="free") + 
        theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
        scale_y_log10(breaks=c(0, 0.1, 0.25, 0.5, 0.75, 1.0)) + 
        theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=25),
              axis.text.y = element_text(size=20), strip.text=element_text(size=35),
              legend.text = element_text(size=25), legend.title=element_text(size=35),
              axis.title=element_text(size=35),
              plot.title=element_text(hjust=0.5, size=30)))
    dev.off()    
    
    ## now make a summary split by tag region and whether there is a motif hit or not
    expanded_tag_region_motif_enh_summary <- ddply(top_coloc_enh_motif_overlaps.expanded, .(tag_region), function(x) {
        motif_hits <- x[x$num_motif_hits > 0,]
        non_motif_hits <- x[x$num_motif_hits==0,]

        ## summarize the characteristics of the motif hits
        motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(motif_hits),
            enh_overlap_hits = sum(motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(motif_hits$any_enh_overlap=="no"),
            enh_tiss_match = sum(motif_hits$f5_category_match=="yes" |
                motif_hits$hmm_category_match=="yes"),
            enh_tiss_non_match = sum(motif_hits$any_enh_overlap=="yes" &
                motif_hits$f5_category_match=="no" & motif_hits$hmm_category_match=="no"))
        
        ## summarize the characteristics of the non-motif hits
        non_motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(non_motif_hits),
            enh_overlap_hits = sum(non_motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(non_motif_hits$any_enh_overlap=="no"),
            enh_tiss_match = sum(non_motif_hits$f5_category_match=="yes" |
                non_motif_hits$hmm_category_match=="yes"),
            enh_tiss_non_match = sum(non_motif_hits$any_enh_overlap=="yes" &
                non_motif_hits$f5_category_match=="no" & non_motif_hits$hmm_category_match=="no"))

        return(rbind(cbind(hit_motif="motif_overlap", motif_hit_summary),
                     cbind(hit_motif="no_motif_overlap", non_motif_hit_summary)))
    })

    ## add in 0 counts for the regions we did not observe
    expanded_tag_region_motif_enh_summary <- rbind(expanded_tag_region_motif_enh_summary,
                                          data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0),
                                          data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))
    
    ## write this summary table out
    write.table(expanded_tag_region_motif_enh_summary, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_motif_enh_tag_region_summary.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## melt this down
    melted_expanded_enh_motif_summ <- melt(expanded_tag_region_motif_enh_summary, id.vars=c("tag_region", "hit_motif"), measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

    melted_expanded_enh_motif_summ$variable <- factor(melted_expanded_enh_motif_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_and_motif_overlaps_barplot_', coloc_prob_thresh, '_prob_thresh'), width_ratio=2.0)
    print(ggplot(melted_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_prob_thresh, "probability")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=3, title="")) +
    facet_wrap(~hit_motif, nrow=1, scales="free") + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
    dev.off()

    ## make a more specific summary that checks whether the consistent tissue category is one of the ones we care about    
    expanded_relevant_class_tag_region_motif_enh_summary <- ddply(top_coloc_enh_motif_overlaps.expanded, .(tag_region), function(x) {
        motif_hits <- x[x$num_motif_hits > 0,]
        non_motif_hits <- x[x$num_motif_hits==0,]

        ## summarize the characteristics of the motif hits
        motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(motif_hits),
            enh_overlap_hits = sum(motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(motif_hits$any_enh_overlap=="no"),
            enh_relevant_tiss_match = sum((motif_hits$f5_category_match=="yes" |
                motif_hits$hmm_category_match=="yes") &
                motif_hits$gtex_tissue_class %in% relevant_classes),
            enh_irrelevant_tiss_match = sum((motif_hits$f5_category_match=="yes" |
                motif_hits$hmm_category_match=="yes") &
                !(motif_hits$gtex_tissue_class %in% relevant_classes)),
            enh_tiss_non_match = sum(motif_hits$any_enh_overlap=="yes" &
                motif_hits$f5_category_match=="no" & motif_hits$hmm_category_match=="no"))
        
        ## summarize the characteristics of the non-motif hits
        non_motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(non_motif_hits),
            enh_overlap_hits = sum(non_motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(non_motif_hits$any_enh_overlap=="no"),
            enh_relevant_tiss_match = sum((non_motif_hits$f5_category_match=="yes" |
                non_motif_hits$hmm_category_match=="yes") &
                non_motif_hits$gtex_tissue_class %in% relevant_classes),
            enh_irrelevant_tiss_match = sum((non_motif_hits$f5_category_match=="yes" |
                non_motif_hits$hmm_category_match=="yes") &
                !(non_motif_hits$gtex_tissue_class %in% relevant_classes)),
            enh_tiss_non_match = sum(non_motif_hits$any_enh_overlap=="yes" &
                non_motif_hits$f5_category_match=="no" & non_motif_hits$hmm_category_match=="no"))

        return(rbind(cbind(hit_motif="motif_overlap", motif_hit_summary),
                     cbind(hit_motif="no_motif_overlap", non_motif_hit_summary)))
    })

    ## add in 0 counts for the regions we did not observe
    expanded_relevant_class_tag_region_motif_enh_summary <- rbind(expanded_relevant_class_tag_region_motif_enh_summary,
                                          data.frame(tag_region=all_regions[!(all_regions %in% expanded_relevant_class_tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_relevant_tiss_match=0, enh_irrelevant_tiss_match=0, enh_tiss_non_match=0),
                                          data.frame(tag_region=all_regions[!(all_regions %in% expanded_relevant_class_tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_relevant_tiss_match=0, enh_irrelevant_tiss_match=0, enh_tiss_non_match=0))
    
    ## write this summary table out
    write.table(expanded_relevant_class_tag_region_motif_enh_summary, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_motif_enh_tag_region_relevant_class_summary.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## melt this down
    melted_relevant_class_expanded_enh_motif_summ <- melt(expanded_relevant_class_tag_region_motif_enh_summary, id.vars=c("tag_region", "hit_motif"), measure.vars=c("enh_tiss_non_match", "enh_relevant_tiss_match", "enh_irrelevant_tiss_match", "enh_non_overlap_hits"), value.name="count")

    melted_relevant_class_expanded_enh_motif_summ$variable <- factor(melted_relevant_class_expanded_enh_motif_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_irrelevant_tiss_match", "enh_relevant_tiss_match"), ordered=T)

    ## make two barplots: one vertical (for paper) and one horizontal (for talks)
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_enh_and_motif_overlaps_barplot_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5)
    print(ggplot(melted_relevant_class_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in other tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_relevant_tiss_match"="red", "enh_irrelevant_tiss_match"="firebrick3", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_prob_thresh, "probability")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=2, title="")) +
    facet_wrap(~factor(hit_motif, levels=c("motif_overlap", "no_motif_overlap"),
                       labels=c("Motif Overlap", "No Motif Overlap")),
                       nrow=2, scales="free") + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_enh_and_motif_overlaps_horizontal_barplot_', coloc_prob_thresh, '_prob_thresh'), width_ratio=1.5)
    print(ggplot(melted_relevant_class_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in other tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_relevant_tiss_match"="red", "enh_irrelevant_tiss_match"="firebrick3", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_prob_thresh, "probability")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=2, title="")) +
    facet_wrap(~factor(hit_motif, levels=c("motif_overlap", "no_motif_overlap"),
                       labels=c("Motif Overlap", "No Motif Overlap")),
                       ncol=2, scales="free") + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
    dev.off()

} else {
    cat("This analysis requires HOMER PWM disruption calculations\n")
}

## -----------------------------------------------------------------------------
## 7. Compare GTEx COLOC hits with direct eQTL overlaps
## -----------------------------------------------------------------------------
## note: i wrote this section before realizing that I could directly pull out the beta directions
## from the COLOC _full_results.txt files, so this section has fewer variants than the eQTL
## plots in section 6, but has p-value information
if(check_param(param_ref, 'gtex_dir')) {
    eqtl_overlap_file <- paste0(param_ref[['outdir']], '/analysis_results/gtex_eqtl_overlap/tables/',
                                param_ref[['outprefix']], "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh,
                                "_uniq_snp_eqtl_overlaps_no_tagsnp_info.txt")
    
    uniq_snp_eqtl_overlaps <- read.table(eqtl_overlap_file, header=T, sep="\t", quote="", as.is=T)

    ## -----------------------------------    
    ## start with the single top colocalized SNP data
    ## -----------------------------------
    ## we want to add information about the eQTL effect direction to the data frame
    top_coloc_enh_motif_eqtl_overlaps <- data.frame(stringsAsFactors = F)
    ## indexing vector
    non_enh_motif_cols <- colnames(top_coloc_enh_motif_overlaps) %in% colnames(top_coloc_hits) |
        colnames(top_coloc_enh_motif_overlaps) %in% c("top_coloc_snp", "max_coloc_prob")
    
    ## now go through and find eQTL overlaps
    for(i in seq(nrow(top_coloc_enh_motif_overlaps))) {
        this_data <- top_coloc_enh_motif_overlaps[i,]

        ## find the eQTL data for the top SNP here. note that we aren't looking for all the
        ## eQTLs that this SNP has, just the one we found to be colocalized
        ## so this should give us exactly one hit
        this_eqtl_overlap <- uniq_snp_eqtl_overlaps[with(uniq_snp_eqtl_overlaps,
                                                         rsID==this_data$top_coloc_snp & tissue==this_data$tissue & gene_name==this_data$eqtl_gene_name),]
        ## only include it in the output if there are any hits
        if(nrow(this_eqtl_overlap)==1) {
            ## now make a summary and add it to the DF that has the other info as well
            this_data <- cbind(this_data,
                               beta=this_eqtl_overlap$beta,
                               t_stat=this_eqtl_overlap$t_stat,
                               se=this_eqtl_overlap$se,
                               p_value=this_eqtl_overlap$p_value)
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_eqtl_overlaps <- rbind(top_coloc_enh_motif_eqtl_overlaps, this_data)
        } else if(nrow(this_eqtl_overlap)>1) {
            cat("There were multiple hits for SNP", as.character(this_data$top_coloc_snp), "tissue", as.character(this_data$tissue), "and gene", as.character(this_data$eqtl_gene_name), "\n")
        } else {
            ## if there are no overlaps, note that in the summary
            ## now make a summary and add it to the DF with the other info as well
            this_data <- cbind(this_data, beta=NA, t_stat=NA, se=NA, p_value=NA)
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_eqtl_overlaps <- rbind(top_coloc_enh_motif_eqtl_overlaps, this_data)
        }
        
    }

    cat(sum(!is.na(top_coloc_enh_motif_eqtl_overlaps$beta)), "out of", nrow(top_coloc_enh_motif_eqtl_overlaps), "colocalization signals were found in direct eQTL overlap dataset\n")
    
    ## write this out
    write.table(top_coloc_enh_motif_eqtl_overlaps, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_motif_eqtl_info.txt'), quote=F, sep="\t", row.names=F)
    
    ## -----------------------------------    
    ## next do the expanded SNP set
    ## -----------------------------------
    ## we want to add information about the eQTL effect direction to the data frame
    top_coloc_enh_motif_eqtl_overlaps.expanded <- data.frame(stringsAsFactors = F)
    ## indexing vector
    non_enh_motif_cols.expanded <- colnames(top_coloc_enh_motif_overlaps.expanded) %in% colnames(top_coloc_hits) |
        colnames(top_coloc_enh_motif_overlaps.expanded) %in% c("num_prob_expanded_snps", "high_coloc_snp", "coloc_snp_prob")
    
    ## now go through and find pwm overlaps
    for(i in seq(nrow(top_coloc_enh_motif_overlaps.expanded))) {
        this_data <- top_coloc_enh_motif_overlaps.expanded[i,]

        ## find the eQTL data for this particular SNP. note that we aren't looking for all the
        ## eQTLs that this SNP has, just the one we found to be colocalized
        ## so this should give us exactly one hit
        this_eqtl_overlap <- uniq_snp_eqtl_overlaps[with(uniq_snp_eqtl_overlaps,
                                                         rsID==this_data$high_coloc_snp & tissue==this_data$tissue & gene_name==this_data$eqtl_gene_name),]
        ## only include it in the output if there are any hits
        if(nrow(this_eqtl_overlap)==1) {
            ## now make a summary and add it to the DF that has the other info as well
            this_data <- cbind(this_data,
                               beta=this_eqtl_overlap$beta,
                               t_stat=this_eqtl_overlap$t_stat,
                               se=this_eqtl_overlap$se,
                               p_value=this_eqtl_overlap$p_value)
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_eqtl_overlaps.expanded <- rbind(top_coloc_enh_motif_eqtl_overlaps.expanded, this_data)
        } else if(nrow(this_eqtl_overlap)>1) {
            cat("There were multiple hits for SNP", as.character(this_data$top_coloc_snp), "tissue", as.character(this_data$tissue), "and gene", as.character(this_data$eqtl_gene_name), "\n")
        } else {
            ## if there are no overlaps, note that in the summary
            ## now make a summary and add it to the DF with the other info as well
            this_data <- cbind(this_data, beta=NA, t_stat=NA, se=NA, p_value=NA)
            
            ## now add this row to the summary DF
            top_coloc_enh_motif_eqtl_overlaps.expanded <- rbind(top_coloc_enh_motif_eqtl_overlaps.expanded, this_data)
        }
        
    }

    cat(sum(!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$beta)), "out of", nrow(top_coloc_enh_motif_eqtl_overlaps.expanded), "ABF-expanded colocalization signals were found in direct eQTL overlap dataset\n")
    ## also ask how many unique variants this was
    cat(length(unique(top_coloc_enh_motif_eqtl_overlaps.expanded$high_coloc_snp[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$beta)])), "out of",
        length(unique(top_coloc_enh_motif_eqtl_overlaps.expanded$high_coloc_snp)), "unique ABF-expanded SNPs had at least one significant direct eQTL overlap\n")
    
    ## i also want to know how many of the unique colocalization signals had any overlap
    expanded_coloc_eqtl_overlap_count <- ddply(top_coloc_enh_motif_eqtl_overlaps.expanded, .(tag_region, tissue, eqtl_gene_name), summarize, eqtl_overlap = any(!is.na(beta)))
    table(expanded_coloc_eqtl_overlap_count$eqtl_overlap)
    
    cat(sum(expanded_coloc_eqtl_overlap_count$eqtl_overlap), "out of", nrow(expanded_coloc_eqtl_overlap_count), "colocalization signals were found in direct eQTL overlap dataset after ABF expansion\n")    
    
    ## write this out
    write.table(top_coloc_enh_motif_eqtl_overlaps.expanded, paste0(outdir, '/tables/coloc_analysis/', outprefix, '_gtex_coloc_enh_motif_eqtl_info.', coloc_prob_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## check whether the eQTLs in each comparison line up or not:
    expanded_direct_eqtl_info <- ddply(top_coloc_enh_motif_eqtl_overlaps.expanded, .(tag_region, tissue, gtex_tissue_class, eqtl_gene_name), summarize, eqtl_overlaps = sum(!is.na(beta)), eqtl_overlap_prop = sum(!is.na(beta)) / length(beta), positive_betas = sum(beta > 0, na.rm=T), negative_betas = sum(beta < 0, na.rm=T))

    cat(sum(expanded_direct_eqtl_info$eqtl_overlaps!=0), "comparisons have eQTL overlap and", sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas > 0 & negative_betas > 0)),
        "of these have variants with different beta directions\n")
    cat("Of the rest,", sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas > 0 & negative_betas==0)), "are all positive and",
        sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas==0 & negative_betas > 0)), "are all negative\n")

    ## make a plot of the beta distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_beta_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$beta),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=beta, color=beta)) +
    xlab("Tissue; target gene; tag region") + ylab("Beta value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap beta distributions, expanded to", coloc_prob_thresh, "probability")) + 
    geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## make a plot of the t_stat distributions across comparisons
    ## set the limits first
    tstat_min <- round_any(min(top_coloc_enh_motif_eqtl_overlaps.expanded$t_stat, na.rm=T), accuracy=1, f=floor)
    tstat_max <- round_any(max(top_coloc_enh_motif_eqtl_overlaps.expanded$t_stat, na.rm=T), accuracy=1, f=ceiling)
    
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_t_stat_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$t_stat),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=t_stat, color=t_stat)) +
    xlab("Tissue; target gene; tag region") + ylab("T statistics value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap T-statistic distributions, expanded to", coloc_prob_thresh, "probability")) +
    scale_y_continuous(breaks=seq(tstat_min, tstat_max, by=1)) + 
    geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## make a plot of the se distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_se_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$se),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=se, color=se)) +
    xlab("Tissue; target gene; tag region") + ylab("Standard error value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap SE distributions, expanded to", coloc_prob_thresh, "probability")) + 
    geom_point() + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## make a plot of the p value distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_p_value_dists_', coloc_prob_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$p_value),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=p_value, color=p_value)) +
    xlab("Tissue; target gene; tag region") + ylab("P-value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap P value distributions, expanded to", coloc_prob_thresh, "probability")) + 
    geom_point() + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()    
    
}

