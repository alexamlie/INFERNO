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

## read in parameters:
args <- commandArgs(trailingOnly=TRUE)
cat(args, '\n')   
cat(length(args), '\n')

if (length(args)==21) {
    ## the main output directory
    outdir <- args[1]
    ## set up the INFERNO directory for enhancer, motif, and direct eQTL overlaps
    inferno_param_f <- args[2]
    ## define the P(H_4) threshold for finding the strongest hits
    coloc_h4_thresh <- as.numeric(args[3])
    ## also get the probability threshold for ABF expansion
    coloc_abf_thresh <- as.numeric(args[4])
    ## get the top GWAS SNPs (this is the file after LD pruning, to make sure the number of tag
    ## regions are consistent across analysis scripts)
    top_snpf <- args[5]
    ## the file containing summary statistics
    gwas_summary_file <- args[6]
    ## the GTEx directory containing the full eQTL statistics for each tissue
    ## this MUST be sorted by position and split into individual files per chromosome
    ## (gtex_download_and_sort_full_v6p_data.sh)
    gtex_dir <- args[7]
    ## sample sizes file for GTEx (csv)
    gtex_sample_sizef <- args[8]
    ## the GTEx tissue class file for INFERNO
    gtex_class_file <- args[9]
    ## define the file that is used to match GTEx IDs with rsIDs
    gtex_rsid_match <- args[10]
    ## reference for gene names
    hg19_ensembl_ref_file <- args[11]
    ## define the relevant tissue classes. must be a list with no spaces e.g. "("Blood","Brain")"
    ## i try to make this flexible to different types of quotes
    relevant_classes <- strsplit(gsub("\'|\"", "", args[12]), ",")[[1]]
    rsid_col <- as.numeric(args[13])
    pos_col <- as.numeric(args[14])
    pval_col <- as.numeric(args[15])
    chr_col <- as.numeric(args[16])
    allele1_col <- as.numeric(args[17])
    allele2_col <- as.numeric(args[18])
    maf_col <- as.numeric(args[19])
    case_prop <- as.numeric(args[20])
    sample_size <- as.numeric(args[21])
} else if (length(args)==22) {
    ## the main output directory
    outdir <- args[1]
    ## set up the INFERNO directory for enhancer, motif, and direct eQTL overlaps
    inferno_param_f <- args[2]
    ## define the P(H_4) threshold for finding the strongest hits
    coloc_h4_thresh <- as.numeric(args[3])
    ## also get the probability threshold for ABF expansion
    coloc_abf_thresh <- as.numeric(args[4])
    ## get the top GWAS SNPs (this is the file after LD pruning, to make sure the number of tag
    ## regions are consistent across analysis scripts)
    top_snpf <- args[5]
    ## the file containing summary statistics
    gwas_summary_file <- args[6]
    ## the GTEx directory containing the full eQTL statistics for each tissue
    ## this MUST be sorted by position and split into individual files per chromosome
    ## (gtex_download_and_sort_full_v6p_data.sh)
    gtex_dir <- args[7]
    ## sample sizes file for GTEx (csv)
    gtex_sample_sizef <- args[8]
    ## the GTEx tissue class file for INFERNO
    gtex_class_file <- args[9]
    ## define the file that is used to match GTEx IDs with rsIDs
    gtex_rsid_match <- args[10]
    ## reference for gene names
    hg19_ensembl_ref_file <- args[11]
    ## define the relevant tissue classes. must be a list with no spaces e.g. "("Blood","Brain")"
    ## i try to make this flexible to different types of quotes
    relevant_classes <- strsplit(gsub("\'|\"", "", args[12]), ",")[[1]]
    rsid_col <- as.numeric(args[13])
    pos_col <- as.numeric(args[14])
    pval_col <- as.numeric(args[15])
    chr_col <- as.numeric(args[16])
    allele1_col <- as.numeric(args[17])
    allele2_col <- as.numeric(args[18])
    maf_col <- as.numeric(args[19])
    case_prop <- as.numeric(args[20])
    sample_size <- as.numeric(args[21])
    ## for locuszoom analysis, set the path to the binary:
    locuszoom <- args[22]
} else {
    ## TODO: fill in and make sure to mention that GTEx files MUST be sorted
    stop("Requires 21 or 22 arguments")    
}

dir.create(paste0(outdir, "/plots/"), F, T)
dir.create(paste0(outdir, "/tables/"), F, T)

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

## set up the summary output
summary_file <- paste0(output_dir, '/tables/', param_ref[['outprefix']], '_colocalization_summary.txt')

## read in the SNPs that were analyzed by INFERNO, to cross reference
ld_stats_file <- paste0(param_ref[['outdir']], '/ld_expansion/', param_ref[['outprefix']],
                        "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, ".txt")
ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)

## read in the list of the top GWAS SNPs in INFERNO format
top_snps <- read.table(top_snpf, header=F, sep="\t", as.is=T,
                            col.names=c("chr", "rsID", "region", "pos"))

## based on the number of regions, set up the graphics function
if(length(unique(top_snps$region)) > 30) {
    default_gfx_width <- 20
    default_gfx_height <- 20
} else {
    default_gfx_width <- 10
    default_gfx_height <- 10
}
## also define a text size multiplier based on these defaults
text_mult <- default_gfx_width / 15

## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
    if(type=='pdf') {
        pdf(file=paste0(filename, ".pdf"), width=default_gfx_width*width_ratio, height=default_gfx_height*height_ratio, pointsize=12, onefile=FALSE)
    } else if(type=='png') {
        ## use type='cairo' for when X11 doesn't work
        png(filename=paste0(filename, ".png"), width=default_gfx_width*width_ratio, height=default_gfx_height*height_ratio, res=300, units='in', type='cairo')
    } else {
        cat('filetype not supported\n')
    }
}

## read in the sample sizes
gtex_sample_sizes <- read.csv(gtex_sample_sizef)
colnames(gtex_sample_sizes) <- c("Tissue", "eqtl_samp_num", "rnaseq_samp_num", "egene_num")
gtex_sample_sizes$tissmatch <- gsub("- ", "", gsub("\\(|\\)", "", gtex_sample_sizes$Tissue))

## also read in the GTEx tissue classes
gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
## to match up with the GTEx samples
gtex_category_df$coloc_match <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data, fixed=T)

## finally, read in the reference to get gene names
ensembl_xref <- read.table(hg19_ensembl_ref_file, header=T, sep="\t", quote="", comment.char="", as.is=T, col.names=c("tx_id", "gene_id", "gene_name"))

## ## read in the header of the summary stats file so that i know the column names for locuszoom
## summary_f_header <- strsplit(readLines(gwas_summary_file, n=1), "\t")[[1]]
## rsid_colname <- summary_f_header[rsid_col]
## pval_colname <- summary_f_header[pval_col]

## -----------------------------------------------------------------------------
## 3. Run COLOC analysis across all top GWAS loci and GTEx data
## -----------------------------------------------------------------------------
cat("Performing colocalization analysis\n")
{
coloc_start_time <- proc.time()
## save the GTEx eQTL column names
## this is especially important since the chromosome-split files never have headers
gtex_eqtl_header <- c("gene_id", "variant_id", "tss_distance", "pval_nominal", "slope", "slope_se")

## go through all the chromosomes since we read in the GTEx data by chromosome
for(this_chr in unique(top_snps$chr)) {
    cat("Parsing", this_chr, "\n")

    ## read in the summary data from this chromosome
    this_gwas_data <- read.table(pipe(paste0("awk -F$'\t' '$", chr_col, "==\"", this_chr, "\" || $", chr_col, "==\"", gsub("chr", "", this_chr), "\"' ", gwas_summary_file)), header=F, sep="\t", quote="", as.is=T)
    
    ## read in the GTEx rsID matches for this chromosome    
    this_chr_gtex_rsid_match <- read.table(pipe(paste0("awk -F$'\t' '$1==\"Chr\" || $1==",
                                                       gsub("chr", "", this_chr),
                                                       "' ", gtex_rsid_match)),
                                                header=T, sep="\t", quote="", as.is=T)
    
    ## now analyze each tag region
    for(this_tag in unique(top_snps$region[top_snps$chr==this_chr])) {
        cat("Analyzing", this_tag, "region\n")
        dir.create(paste0(outdir, '/tables/gtex_coloc/', gsub("/", "_", this_tag, fixed=T)), F, T)

        this_tag_data <- top_snps[top_snps$region==this_tag,]
        ## for multiple variants in a region, just use the first one
        if(nrow(this_tag_data) > 1) {
            cat(nrow(this_tag_data), "variants found in", this_tag, "region, using the first one\n")
            this_tag_data <- this_tag_data[1,]
        }
        
        this_tag_gwas_data <- this_gwas_data[this_gwas_data[,rsid_col]==this_tag_data$rsID,]
        
        if(nrow(this_tag_gwas_data) != 1) {
            cat("Could not find single match for tag snp", this_tag_data$rsID, "in GWAS data\n")
            next
        }

        this_region_gwas_snps <- this_gwas_data[abs(this_gwas_data[,pos_col] - this_tag_gwas_data[,pos_col]) <= 500000,]
        cat(nrow(this_region_gwas_snps), "SNPs found in GWAS around", this_tag_data$rsID, "\n")

        ## add a column for the variant IDs in the GTEx format
        ## ({chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37)
        ## GWAS alleles may not match up with GTEx, so we have to try both ways
        this_region_gwas_snps$gtex_id1 <-
            ## allele 1 is minor
            paste(gsub("chr", "", this_chr), this_region_gwas_snps[,pos_col],
                  toupper(this_region_gwas_snps[,allele2_col]),
                  toupper(this_region_gwas_snps[,allele1_col]), "b37", sep="_")
        
        this_region_gwas_snps$gtex_id2 <-
            ## allele 2 is minor
            paste(gsub("chr", "", this_chr), this_region_gwas_snps[,pos_col],
                  toupper(this_region_gwas_snps[,allele1_col]),
                  toupper(this_region_gwas_snps[,allele2_col]), "b37", sep="_") 

        ## define the tag variant GTEx ID as well
        ## start by trying to use allele 1 as the reference allele, but it gets checked
        this_tag_gtex_id <-
            paste(gsub("chr", "", this_chr), this_tag_gwas_data[,pos_col],
                  toupper(this_tag_gwas_data[,allele1_col]),
                  toupper(this_tag_gwas_data[,allele2_col]), "b37", sep="_")         
        
        ## now we have to loop through each GTEx tissue and define the gene sets to perform
        ## colocalization analysis on. first loop through the top-level tissues and then pull
        ## out the chromosome data for that
        for(gtex_tiss in list.files(gtex_dir)) {
            gtex_tiss_match <- gsub("_", " ", gtex_tiss)
            cat("Analyzing tissue", gtex_tiss_match, "\n")

            gtex_file <- paste0(gtex_dir, "/", gtex_tiss, "/", gtex_tiss,
                                "_Analysis.v6p.all_snpgene_pairs.", this_chr, ".txt.gz")

            ## try to find all the genes tested with this tag SNP
            ## first set up a fast awk call
            ## note that this relies on the data being sorted by position
            awk_call <- paste0("awk -F$'\t' 'BEGIN{SNP_OBS=\"\"} {if($2==\"",
                               this_tag_gtex_id, "\") {SNP_OBS=\"Yes\"; print $1} else if(SNP_OBS) {exit;}}'")

            all_genes <- system2("zcat", paste0(gtex_file, " | ", awk_call, " | sort -u"), stdout=TRUE)
                        
            ## check that we actually found some genes, and flip the alleles if not
            if(length(all_genes)==0) {
                ## we only need to change the ID we use to search here because the direction of
                ## effect is taken care of when we search the IDs of the full GWAS dataset
                cat("GWAS minor allele for tag SNP", this_tag_gwas_data[,rsid_col], "does not match up with GTEx ID (no genes found using GWAS minor allele); flipping direction and re-testing\n")
                
                ## if we thought allele 1 was major, it's actually minor
                this_tag_gtex_id <-
                    paste(gsub("chr", "", this_chr), this_tag_gwas_data[,pos_col],
                          toupper(this_tag_gwas_data[,allele2_col]),
                          toupper(this_tag_gwas_data[,allele1_col]), "b37", sep="_") 

                awk_call <- paste0("awk -F$'\t' 'BEGIN{SNP_OBS=\"\"} {if($2==\"",
                                   this_tag_gtex_id, "\") {SNP_OBS=\"Yes\"; print $1} else if(SNP_OBS) {exit;}}'")
                
                all_genes <- system2("zcat", paste0(gtex_file, " | ", awk_call, " | sort -u"), stdout=TRUE)

                ## if this is still 0, just continue
                if(length(all_genes)==0) {
                    cat("No genes found for", this_tag_gwas_data[,rsid_col], this_tag, "in tissue",
                        gtex_tiss_match, "!\n")
                    next
                }
            }

            ## write the gene list to a file for awk purposes
            gene_outf <- paste0(outdir, '/tables/gtex_coloc/', gsub("/", "_", this_tag, fixed=T), "/", gtex_tiss, "_gene_ids.txt")
            write.table(all_genes, gene_outf, quote=F, sep="\t", row.names=F, col.names=F)
            
            ## read in all the data for all these genes at once
            all_eqtl_awk_call <- paste0("awk -F$'\t' 'NR==FNR {gene_ids[$1]; next} {if($1 in gene_ids) {print $0}}' ")
            
            all_eqtl_data <- data.frame(do.call(rbind, strsplit(
                system2("zcat", paste0(gtex_file,
                                       " | ", all_eqtl_awk_call,
                                       gene_outf, " - "), stdout=TRUE),
                split="\t")), stringsAsFactors=F)
            colnames(all_eqtl_data) <- gtex_eqtl_header

            cat("Found", length(unique(all_eqtl_data$gene_id)), "genes tested in tissue",
                gtex_tiss_match, "taking up",
                format(object.size(all_eqtl_data), units="auto"), "with", nrow(all_eqtl_data),
                "total eQTL pairs tested\n")

            ## only make the directory if we found any genes to test
            this_tiss_outdir <- paste0(outdir, '/tables/gtex_coloc/', gsub("/", "_", this_tag, fixed=T), "/", gtex_tiss, "/")
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
                gtex_id1_match <- this_region_gwas_snps$gtex_id1 %in% this_eqtl_data$variant_id
                gtex_id2_match <- this_region_gwas_snps$gtex_id2 %in% this_eqtl_data$variant_id
                this_region_gwas_snps.parsed <- this_region_gwas_snps[gtex_id1_match | gtex_id2_match,]
                
                ## set a new column for the actual matching
                ## we use this matching a lot so save it 
                gtex_id1_parsed_match <- this_region_gwas_snps.parsed$gtex_id1 %in% this_eqtl_data$variant_id
                this_region_gwas_snps.parsed$gtex_id <- ifelse(gtex_id1_parsed_match, this_region_gwas_snps.parsed$gtex_id1, this_region_gwas_snps.parsed$gtex_id2)

                ## TODO: if we're using odds ratios, flip the effect directions if necessary
                
                ## also parse the eQTL data
                this_eqtl_data.parsed <- this_eqtl_data[this_eqtl_data$variant_id %in% this_region_gwas_snps.parsed$gtex_id,]
                
                ## check that the number of SNPs are equal
                cat("Analyzing", nrow(this_region_gwas_snps.parsed), "SNPs from GWAS and", nrow(this_eqtl_data.parsed), "matching SNPs tested for eQTL\n")
                if(nrow(this_region_gwas_snps.parsed)!=nrow(this_eqtl_data.parsed)) {
                    cat("Unequal number of SNPs in region", this_tag, "and gene", gene_name, this_gene_id, "\n")
                    next
                }
                
                ## check that everything is in the same order
                if(sum(this_eqtl_data.parsed$variant_id != this_region_gwas_snps.parsed$gtex_id) != 0) {
                    cat("Wrong order in region", this_tag, "and gene", gene_name, this_gene_id, ", attempting to resort\n")
                    this_region_gwas_snps.parsed <- this_region_gwas_snps.parsed[order(this_region_gwas_snps.parsed[,pos_col]),]

                    if(sum(this_eqtl_data.parsed$variant_id != this_region_gwas_snps.parsed$gtex_id) != 0) {
                        cat("Still has out of order variants. Skipping\n")
                        next
                    } else {
                        cat("Successfully resorted variants!\n")
                    }
                }

                this_gtex_samp_size <- gtex_sample_sizes$eqtl_samp_num[gtex_sample_sizes$tissmatch==gtex_tiss_match]
                
                this_pval_coloc_results <- coloc.abf(
                    dataset1=list(pvalues=this_region_gwas_snps.parsed[,pval_col],
                        type="cc", s=case_prop, N=sample_size),
                    dataset2=list(pvalues=this_eqtl_data.parsed$pval_nominal,
                        type="quant", N=this_gtex_samp_size),
                    MAF=this_region_gwas_snps.parsed[,maf_col])
                
                write.table(this_pval_coloc_results[['summary']],
                            file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_summary.txt"),
                            quote=F, sep="\t", col.names=F)
                                
                ## before writing out the full results, need to add a column so that we can
                ## actually match up the SNPs
                ## define the indexing vector
                coloc_snp_idx <- as.numeric(unlist(lapply(strsplit(as.character(this_pval_coloc_results[['results']]$snp), ".", fixed=T), "[[", 2)))
                
                this_pval_coloc_results[['results']]$rsID <- as.character(this_region_gwas_snps.parsed[coloc_snp_idx, rsid_col])
                
                write.table(this_pval_coloc_results[['results']],
                            file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_full_results.txt"),
                            quote=F, sep="\t", col.names=T, row.names=F)

                ## if this comparison meets the threshold and we have a path to locuszoom, make
                ## locuszoom plots:
                if (this_pval_coloc_results[['summary']]['PP.H4.abf'] >= coloc_h4_thresh & exists("locuszoom")) {
                    ## annotate the GTEx data with the matching rsIDs
                    this_eqtl_data.parsed$dbsnp135_rsid <- this_chr_gtex_rsid_match$RS_ID_dbSNP135_original_VCF[match(this_eqtl_data.parsed$variant_id, this_chr_gtex_rsid_match$VariantID)]
                    ## now write out that info
                    write.table(this_eqtl_data.parsed[,c("dbsnp135_rsid", "pval_nominal")],
                                paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_gtex_pvals.txt"),
                                quote=F, sep="\t", col.names=T, row.names=F)

                    ## now make a folder within the plots for this tag region
                    this_tag_lz_out <- paste0(outdir, '/plots/locuszoom_plots/', gsub("/", "_", this_tag, fixed=T), "_region/")
                    dir.create(this_tag_lz_out, F, T)

                    ## do the GWAS analysis, only if we didn't already (this only needs to be run once, even with different GTEx genes and tissues)
                    if(!file.exists(paste0(this_tag_lz_out, gsub("/", "_", this_tag, fixed=T), "_region_gwas_locuszoom_", this_tag_data$rsID, ".pdf"))) {
                        colnames(this_region_gwas_snps.parsed)[rsid_col] <- "MarkerName"
                        colnames(this_region_gwas_snps.parsed)[pval_col] <- "P-value"
                        write.table(this_region_gwas_snps.parsed[,c(rsid_col, pval_col)],
                                    paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_gwas_pvals.txt"),
                                    quote=F, sep="\t", col.names=T, row.names=F)
                        
                        ## make an GWAS locuszoom plot around the tag SNP
                        system2(locuszoom, paste("--metal",
                                                 paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_gwas_pvals.txt"),
                                                 ## "--markercol", rsid_col, "--pvalcol",
                                                 ## pval_col, 
                                                 "--refsnp",
                                                 this_tag_data$rsID, "--flank 500kb --pop EUR --build hg19 --source 1000G_Nov2010 --plotonly --prefix",
                                                 paste0(this_tag_lz_out, gsub("/", "_", this_tag, fixed=T), "_region_gwas_locuszoom"),
                                                 "--no-date legend='right'"))
                    }
                    
                    ## same thing for the eQTL data
                    ## first find the top ABF SNP so that we can label that
                    top_abf_snp <- this_pval_coloc_results[['results']]$rsID[which.max(this_pval_coloc_results[['results']]$SNP.PP.H4)]
                    max_abf <- max(this_pval_coloc_results[['results']]$SNP.PP.H4)

                    ## if it's the same, just label the ABF of the tag variant
                    if(top_abf_snp == this_tag_data$rsID) {
                        write.table(data.frame(snp=this_tag_data$rsID, string=paste0("GWAS tag (", round(max_abf, digits=2), ")"), color="purple"),
                                    file=paste0(this_tiss_outdir, gene_name, "_", this_gene_id, "_lz_labels.txt"),
                                    quote=F, sep="\t", row.names=F)
                        
                    } else {
                        ## find the tag variant ABF
                        tag_abf <- this_pval_coloc_results[['results']]$SNP.PP.H4[this_pval_coloc_results[['results']]$rsID==this_tag_data$rsID]
                        
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
cat("Analyzing high-level colocalization patterns\n")
{
cur_tag <- ""
summ_start <- proc.time()

## i want to read in all the summary data into one unified data frame
all_summ_files <- list.files(paste0(outdir, '/tables/gtex_coloc/'), pattern="*.summary.txt",
                             full.names=T, recursive=T)

all_summary_data <- data.frame(stringsAsFactors = F)
for(f in all_summ_files) {
    ## split the path to get the tag region, tissue, and target gene
    this_info <- strsplit(gsub(outdir, "", gsub("//", "/", f)), "/")[[1]]    

    path_length <- length(this_info)
    
    this_tag <- this_info[path_length-2]
    if(this_tag != cur_tag) {
        cur_tag <- this_tag
        cat("Analyzing", cur_tag, "region\n")
    }
    this_tiss <- this_info[path_length-1]

    this_eqtl_info <- strsplit(this_info[path_length], "_")[[1]]
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
write.table(all_summary_data, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_summaries.txt'), quote=F, sep="\t", row.names=F)

## ## to read this in directly:
## all_summary_data <- read.table(paste0(outdir, '/tables/', outprefix, '_gtex_coloc_summaries.txt'), header=T, sep="\t", quote="")

## make histograms for each of the different hypotheses by tag region
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H0_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H0.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H0 (no causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H0 (no causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H1_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H1.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H1 (GWAS causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H1 (GWAS causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H2_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H2.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H2 (eQTL causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H2 (eQTL causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H3_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H3.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H3 (unshared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H3 (unshared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H4_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H4_highest_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tag_region)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0.5, 1), breaks=seq(0.5, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tag regions, GTEx") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

## also do this split by tissue
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H0_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H0.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H0 (no causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H0 (no causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H1_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H1.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H1 (GWAS causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H1 (GWAS causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H2_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H2.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H2 (eQTL causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H2 (eQTL causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H3_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H3.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H3 (unshared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H3 (unshared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H4_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_H4_highest_prob_hist_by_tissue'), height_ratio=1.5, width_ratio = 2.5)
print(ggplot(all_summary_data, aes(x=PP.H4.abf, fill=tissue)) +
    geom_histogram(binwidth=0.01) +
    scale_x_continuous(limits=c(0.5, 1), breaks=seq(0.5, 1, by=0.05)) +
    facet_wrap(~ tissue, scales="free_y", ncol=5, drop=FALSE) +
    ggtitle("Histograms of H4 (shared causal variant) probability across tissues") + 
    theme_bw() + xlab("Probability of H4 (shared causal variant)") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=20),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

## now make a single histogram looking at the distributions of all hypotheses together
## need to melt the data to do this
melted_summary_data <- melt(all_summary_data, id.vars=1:6, variable.name="hypothesis", value.name="probability")

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_all_hypotheses_prob_density'))
print(ggplot(melted_summary_data, aes(x=probability, fill=hypothesis, color=hypothesis)) +
    geom_density(alpha=0.5) + scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    scale_fill_brewer(palette="Set1") +
    scale_colour_brewer(palette="Set1") + 
    ggtitle("Combined density plot of all hypotheses, GTEx") + 
    theme_bw() + xlab("Probability of hypothesis") + ylab("Density") +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15),
          title=element_text(size=25), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
dev.off()

## -----------------------------------------------------------------------------
## 5. Compare GTEx COLOC hits with top GWAS enhancer overlaps and eQTL effect direction
## -----------------------------------------------------------------------------
## use the colocalization threshold parameter
top_coloc_hits <- all_summary_data[all_summary_data$PP.H4.abf >= coloc_h4_thresh,]
nrow(top_coloc_hits)

## write out the unique genes for pathway analysis
write.table(sort(unique(top_coloc_hits$eqtl_gene_name)), paste0(outdir, '/tables/', outprefix, '_gtex_coloc_all_top_genes.', coloc_h4_thresh, '_thresh.txt'), quote=F, sep="\t", row.names=F, col.names=F)

cat(nrow(unique(top_coloc_hits[,c("tissue", "gtex_tissue_class", "eqtl_gene_name", "eqtl_gene_id")])), "unique tissue-target gene comparisons with P(H_4) >=", coloc_h4_thresh, "\n")

cat(nrow(unique(top_coloc_hits[,c("tissue", "gtex_tissue_class", "eqtl_gene_name", "eqtl_gene_id")])), "unique tissue-target gene comparisons with P(H_4) >=", coloc_h4_thresh, file=summary_file, append=T)

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
## also make one to store more variants, accounting for some amount (coloc_abf_thresh) of the
## individual SNP probability
top_coloc_enh_overlaps.expanded <- data.frame(stringsAsFactors = F)

for(i in seq(nrow(top_coloc_hits))) {
    ## this starts out as just a row of the coloc data but is supplemented
    this_comparison <- top_coloc_hits[i,]
    ## save the original row, for the expanded SNP set analysis
    this_orig_comparison <- this_comparison

    ## to read in the full data, we need to reconstruct the path
    this_data_file <- paste0(outdir, '/tables/gtex_coloc/', this_comparison$tag_region, '/',
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
    snp_set_idx <- indiv_snp_prob_order[1:(sum(snp_cumsum_values <= coloc_abf_thresh)+1)]
    
    top_causal_snp_set <- this_coloc_data$rsID[snp_set_idx]
    expanded_snp_data <- this_coloc_data[snp_set_idx,] 
    
    ## -----------------------------
    ## first do the analysis against the single top SNP
    top_causal_f5_hits <- fantom5_enh_overlaps[fantom5_enh_overlaps$rsID==top_causal_snp_rsid,]
    top_causal_hmm_hits <- roadmap_enh_overlaps[roadmap_enh_overlaps$rsID==top_causal_snp_rsid,]

    ## because we may have one variant in several tag regions, only look at the unique variants
    f5_tag_name_cols <- which(colnames(top_causal_f5_hits) %in% c("tag_name", "tag_no_rsid"))
    top_causal_f5_hits <- unique(top_causal_f5_hits[,-f5_tag_name_cols])
    
    hmm_tag_name_cols <- which(colnames(top_causal_hmm_hits) %in% c("tag_name", "tag_no_rsid"))
    top_causal_hmm_hits <- unique(top_causal_hmm_hits[,-hmm_tag_name_cols])
    
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

        ## because we may have one variant in several tag regions, only look at the unique variants
        f5_tag_name_cols <- which(colnames(this_snp_causal_f5_hits) %in% c("tag_name", "tag_no_rsid"))
        this_snp_causal_f5_hits <- unique(this_snp_causal_hmm_hits[,-f5_tag_name_cols])
        
        hmm_tag_name_cols <- which(colnames(this_snp_causal_hmm_hits) %in% c("tag_name", "tag_no_rsid"))
        this_snp_causal_hmm_hits <- unique(this_snp_causal_hmm_hits[,-hmm_tag_name_cols])

        
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

write.table(top_coloc_enh_overlaps, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_overlaps.txt'), quote=F, sep="\t", row.names=F)

## do the same thing for the expanded data
## add a column for whether there is enhancer overlap or not
top_coloc_enh_overlaps.expanded$any_enh_overlap <- ifelse(top_coloc_enh_overlaps.expanded$num_f5_tissues > 0 | top_coloc_enh_overlaps.expanded$num_hmm_tissues > 0, "yes", "no")
## also add a column for whether the GTEx tissue category is found in roadmap or fantom5
top_coloc_enh_overlaps.expanded$f5_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps.expanded$gtex_tissue_class, x=top_coloc_enh_overlaps.expanded$f5_tissue_classes), "yes", "no")
top_coloc_enh_overlaps.expanded$hmm_category_match <- ifelse(mapply(grepl, pattern=top_coloc_enh_overlaps.expanded$gtex_tissue_class, x=top_coloc_enh_overlaps.expanded$hmm_classes), "yes", "no")
## use these new columns to describe another convenient column for summarizing the tissue matches
top_coloc_enh_overlaps.expanded$enh_levels <- with(top_coloc_enh_overlaps.expanded, ifelse(any_enh_overlap=="no", "no_overlap", ifelse(f5_category_match=="yes" | hmm_category_match=="yes", "tiss_match", "tiss_non_match")))

write.table(top_coloc_enh_overlaps.expanded, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_overlaps.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

cat("Enhancer overlap analysis took", (proc.time() - enh_start_time)[["elapsed"]], "seconds\n")
}

## check how many match up
cat(sum(unique(top_coloc_enh_overlaps$top_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps$top_coloc_snp)),
    "unique max colocalized SNPs match up with INFERNO SNPs\n")
cat(sum(unique(top_coloc_enh_overlaps$top_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps$top_coloc_snp)),
    "unique max colocalized SNPs match up with INFERNO SNPs", file=summary_file, append=T)

cat(sum(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp)),
    "unique high colocalized SNPs in", coloc_abf_thresh, "expanded set match up with INFERNO SNPs\n")
cat(sum(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp) %in% ld_stats_df$rsID), "out of",
    length(unique(top_coloc_enh_overlaps.expanded$high_coloc_snp)),
    "unique high colocalized SNPs in", coloc_abf_thresh, "expanded set match up with INFERNO SNPs", file=summary_file, append=T)

## -----------------------------
## do some visualizations of the single max coloc SNP overlaps
## start with histograms of the maximum colocalization probabilities
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_hist_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(top_coloc_enh_overlaps, aes(x=max_coloc_prob, fill=tag_region)) +
    geom_histogram(binwidth=0.01) +
    scale_x_continuous(limits=c(-0.02, 1.02), breaks=seq(0, 1, by=0.05)) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle("Histograms of maximum SNP causal probability across tag regions, GTEx") + 
    theme_bw() + xlab("Highest probability of SNP being causal") + ylab("Count") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

## now look at the distributions split by annotation overlap
## set the factor levels to get the plots in the right order
top_coloc_enh_overlaps$enh_levels <- factor(top_coloc_enh_overlaps$enh_levels, levels=c("no_overlap", "tiss_non_match", "tiss_match"), ordered=T)

## look at the individual SNP probability distributions by enhancer overlap and tag region
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_overlap_and_tag_region'), height_ratio=1.5, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps, aes(x=enh_levels, y=max_coloc_prob, color=enh_levels)) +
    geom_boxplot() +
    facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
    ggtitle("Distributions of maximum SNP causal probability, GTEx") + 
    theme_bw() + xlab("Type of enhancer overlap") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20)))
dev.off()

## same thing, but collapsed across tag regions
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_max_coloc_snp_prob_distributions_by_enh_overlap'))
print(ggplot(top_coloc_enh_overlaps, aes(x=enh_levels, y=max_coloc_prob, color=enh_levels)) +
    geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
    ggtitle("Distributions of maximum SNP causal probability across tag regions, GTEx") + 
    theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20)))
dev.off()

## next, summarize the counts of highly colocalized SNPs in each region, whether they overlap an enhancer, and whether the category matches or not
tag_region_enh_count_summary <- ddply(top_coloc_enh_overlaps, .(tag_region), summarize,
                                      num_coloc_hits = length(tag_region),
                                      enh_overlap_hits = sum(any_enh_overlap=="yes"),
                                      enh_non_overlap_hits = sum(any_enh_overlap=="no"),
                                      enh_tiss_match = sum(f5_category_match=="yes" | hmm_category_match=="yes"),
                                      enh_tiss_non_match = sum(any_enh_overlap=="yes" & f5_category_match=="no" & hmm_category_match=="no"))

## add in 0 counts for the regions we did not observe (if any)
all_regions <- unique(all_summary_data$tag_region)
if(sum(!(all_regions %in% tag_region_enh_count_summary$tag_region)) > 0) {
    tag_region_enh_count_summary <- rbind(tag_region_enh_count_summary,
                                          data.frame(tag_region=all_regions[!(all_regions %in% tag_region_enh_count_summary$tag_region)], num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))
}
    
## now melt this down into the three categories we care about
melted_enh_count_summ <- melt(tag_region_enh_count_summary, id.vars="tag_region", measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

melted_enh_count_summ$variable <- factor(melted_enh_count_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_overlaps_barplot'))
print(ggplot(melted_enh_count_summ, aes(x=tag_region, y=count, fill=variable)) +
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
          plot.title=element_text(hjust=0.5)))
dev.off()

## eQTL effect direction analysis for the unexpanded set
## number of positive and negative hits
table(ifelse(top_coloc_enh_overlaps$eqtl_beta > 0, "positive", "negative"))

## now plot beta distributions across tag regions
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_beta_dists_across_tag_regions'))
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
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_zscore_dists_across_tag_regions'))
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
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_variance_dists_across_tag_regions'))
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

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_num_snps_to_', coloc_abf_thresh, '_thresh_by_tag_region'), height_ratio=1.5, width_ratio = 1.5)
print(ggplot(unique(top_coloc_enh_overlaps.expanded[,c("tag_region", "tissue", "eqtl_gene_name", "num_prob_expanded_snps")]), aes(x=num_prob_expanded_snps, fill=tag_region)) +
    geom_histogram(binwidth=1) +
    scale_x_continuous(breaks=seq(0, max_snp_set_size, by=ifelse(max_snp_set_size > 50, 10, 2))) +
    facet_wrap(~ tag_region, scales="free_y", ncol=3, drop=FALSE) +
    ggtitle(paste("Histograms of numbers of SNPs to reach", coloc_abf_thresh, "probability")) + 
    theme_bw() + xlab("Number of SNPs to reach probability threshold") + ylab("Number of GTEx - GWAS comparisons") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=25),
          plot.title=element_text(hjust=0.5)))
dev.off()

## set the factor levels to get the plots in the right order
top_coloc_enh_overlaps.expanded$enh_levels <- factor(top_coloc_enh_overlaps.expanded$enh_levels, levels=c("no_overlap", "tiss_non_match", "tiss_match"), ordered=T)

## look at the individual SNP probability distributions by enhancer overlap and tag region
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_overlap_and_tag_region_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
    geom_boxplot() +
    facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
    ggtitle(paste("Distributions of SNP ABFs, expanded to", coloc_abf_thresh, "probability, GTEx")) + 
    theme_bw() + xlab("Type of enhancer overlap") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20)))
dev.off()

## same thing, but collapsed across tag regions
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_overlap_', coloc_abf_thresh, '_prob_thresh'))
print(ggplot(top_coloc_enh_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
    geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
    ggtitle(paste("Distributions of SNP ABFs across tag regions\nexpanded to", coloc_abf_thresh, "probability, GTEx")) + 
    theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          legend.text = element_text(size=15), legend.title=element_text(size=25),
          axis.title=element_text(size=25),
          plot.title=element_text(hjust=0.5, size=20)))
dev.off()

## next, summarize the counts of highly colocalized SNPs in each region, whether they overlap an enhancer, and whether the category matches or not
expanded_tag_region_enh_count_summary <- ddply(top_coloc_enh_overlaps.expanded, .(tag_region), summarize,
                                               num_uniq_snp_comparisons = length(tag_region),
                                               enh_overlap_hits = sum(any_enh_overlap=="yes"),
                                               enh_non_overlap_hits = sum(any_enh_overlap=="no"),
                                               enh_tiss_match = sum(f5_category_match=="yes" | hmm_category_match=="yes"),
                                               enh_tiss_non_match = sum(any_enh_overlap=="yes" & f5_category_match=="no" & hmm_category_match=="no"))

## add in 0 counts for the regions we did not observe
if(sum(!(all_regions %in% expanded_tag_region_enh_count_summary$tag_region)) > 0) {
    expanded_tag_region_enh_count_summary <- rbind(expanded_tag_region_enh_count_summary,
                                                   data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_enh_count_summary$tag_region)], num_uniq_snp_comparisons=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))
}
    
## now melt this down into the three categories we care about
melted_expanded_enh_count_summ <- melt(expanded_tag_region_enh_count_summary, id.vars="tag_region", measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

melted_expanded_enh_count_summ$variable <- factor(melted_expanded_enh_count_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_overlaps_barplot_', coloc_abf_thresh, '_prob_thresh'))
print(ggplot(melted_expanded_enh_count_summ, aes(x=tag_region, y=count, fill=variable)) +
    scale_fill_manual(
        labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
        values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer overlaps of colocalized SNPs expanded to", coloc_abf_thresh, "probability, GTEx")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=3, title="")) + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
dev.off()

## eQTL effect direction analysis for the expanded set
## check whether the eQTLs in each comparison line up or not:
expanded_eqtl_info <- ddply(top_coloc_enh_overlaps.expanded, .(tag_region, tissue, gtex_tissue_class, eqtl_gene_name), summarize, positive_betas = sum(eqtl_beta > 0, na.rm=T), negative_betas = sum(eqtl_beta < 0, na.rm=T))

cat(sum(with(expanded_eqtl_info, positive_betas > 0 & negative_betas > 0)),
    "of", nrow(expanded_eqtl_info), "ABF-expanded sets have variants with different beta directions\n")
cat("Of the rest,", sum(with(expanded_eqtl_info, positive_betas > 0 & negative_betas==0)), "are all positive and",
    sum(with(expanded_eqtl_info, positive_betas==0 & negative_betas > 0)), "are all negative\n")

## now plot beta distributions across comparisons, for the full set of comparisons
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_beta_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_beta, color=eqtl_beta)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL beta values") + theme_bw() +
ggtitle(paste("COLOC eQTL beta distributions, expanded to", coloc_abf_thresh, "probability")) + 
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot z-score distributions across comparisons, for the full set of comparisons
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_zscore_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_z_score, color=eqtl_z_score)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL Z scores") + theme_bw() +
ggtitle(paste("COLOC eQTL Z score distributions, expanded to", coloc_abf_thresh, "probability")) + 
geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## now plot variance distributions across comparisons, for the full set of comparisons
make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_eqtl_variance_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=3.0, width_ratio=1.5)
print(ggplot(top_coloc_enh_overlaps.expanded,
             aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                 y=eqtl_variance, color=eqtl_variance)) +
xlab("Tissue; target gene; tag region") + ylab("eQTL variance") + theme_bw() +
ggtitle(paste("COLOC eQTL variance distributions, expanded to", coloc_abf_thresh, "probability")) + 
geom_point() + coord_flip() + 
facet_grid(tag_region ~ ., scales="free", space="free") +
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
      axis.text.y = element_text(size=15), strip.text=element_blank(),
      title=element_text(size=15), 
      plot.title=element_text(size=20, hjust=0.5)))
dev.off()

## make the same plots for just the relevant tissue categories
if(sum(top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes) > 0) {
    ## plot beta distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_eqtl_beta_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=eqtl_beta, color=eqtl_beta)) +
    xlab("Tissue; target gene; tag region") + ylab("eQTL beta values") + theme_bw() +
    ggtitle(paste("COLOC eQTL beta distributions, expanded to", coloc_abf_thresh, "probability\n",
                  "Only comparisons in", paste(relevant_classes, collapse=", "))) + 
    geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## now plot z-score distributions across comparisons, for the full set of comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_eqtl_zscore_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=eqtl_z_score, color=eqtl_z_score)) +
    xlab("Tissue; target gene; tag region") + ylab("eQTL Z scores") + theme_bw() +
    ggtitle(paste("COLOC eQTL Z score distributions, expanded to", coloc_abf_thresh, "probability\n",
                  "Only comparisons in", paste(relevant_classes, collapse=", "))) +
    geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## now plot variance distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_eqtl_variance_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_overlaps.expanded[top_coloc_enh_overlaps.expanded$gtex_tissue_class %in% relevant_classes,],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=eqtl_variance, color=eqtl_variance)) +
    xlab("Tissue; target gene; tag region") + ylab("eQTL variance") + theme_bw() +
    ggtitle(paste("COLOC eQTL variance distributions, expanded to", coloc_abf_thresh, "probability\n",
                  "Only comparisons in", paste(relevant_classes, collapse=", "))) +
    geom_point() + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()
}

## -----------------------------------------------------------------------------
## 6. Compare GTEx COLOC hits with motif disruptions
## -----------------------------------------------------------------------------
## require PWM disruption calculations
if(check_param(param_ref, "homer_motif_pwm_file") & check_param(param_ref, "homer_motif_seq_file")) {
    cat("Comparing GTEx COLOC hits with HOMER motif disruptions\n")
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
        ## TODO: only use unique SNPs across tag regions here..        
        
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
    write.table(top_coloc_enh_motif_overlaps, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_motif_overlap_summary.txt'), quote=F, sep="\t", row.names=F)

    write.table(top_coloc_motif_overlaps, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_overlaps.txt'), quote=F, sep="\t", row.names=F)
    
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
    if(sum(!(all_regions %in% tag_region_motif_enh_summary$tag_region)) > 0) {
        tag_region_motif_enh_summary <- rbind(tag_region_motif_enh_summary,
                                              data.frame(tag_region=all_regions[!(all_regions %in% tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0),
                                              data.frame(tag_region=all_regions[!(all_regions %in% tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0))
    }
    
    ## write this summary table out
    write.table(tag_region_motif_enh_summary, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_enh_tag_region_summary.txt'), quote=F, sep="\t", row.names=F)
    
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

    ## also add a column for the 'relevant class' enhancer overlaps
    top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels <- with(top_coloc_enh_motif_overlaps.expanded,
        ## check for a concordant hit in a relevant class, otherwise just use the existing enh_levels
        ifelse(enh_levels=="tiss_match", ifelse(gtex_tissue_class %in% relevant_classes, "relevant_tiss_match", "irrelevant_tiss_match"), as.character(enh_levels)))
    
    top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels <- factor(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels, levels=c("no_overlap", "tiss_non_match", "irrelevant_tiss_match", "relevant_tiss_match"), ordered=T)

    ## finally add a column for the 'relevant class' enhancer and motif overlaps
    top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels <- with(top_coloc_enh_motif_overlaps.expanded,
        ## check for a concordant hit in a relevant class, otherwise just use the existing motif levels
        ## need to do this for both motif and non-motif hits
        ifelse(motif_levels=="motif_tiss_match", ifelse(gtex_tissue_class %in% relevant_classes, "motif_relevant_tiss_match", "motif_irrelevant_tiss_match"), ifelse(motif_levels=="no_motif_tiss_match", ifelse(gtex_tissue_class %in% relevant_classes, "no_motif_relevant_tiss_match", "no_motif_irrelevant_tiss_match"), as.character(motif_levels))))   

    top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels <- factor(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels, levels=c("no_motif_no_enh", "no_motif_non_tiss_match", "no_motif_irrelevant_tiss_match", "no_motif_relevant_tiss_match", "motif_no_enh", "motif_non_tiss_match", "motif_irrelevant_tiss_match", "motif_relevant_tiss_match"), ordered=T)
    
    ## write these out
    write.table(top_coloc_enh_motif_overlaps.expanded, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_motif_overlap_summary.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    write.table(top_coloc_motif_overlaps.expanded, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_overlaps.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## now look at the individual SNP probability distributions
    ## set the factor levels to get the plots in the right order
    top_coloc_enh_motif_overlaps.expanded$motif_levels <- factor(top_coloc_enh_motif_overlaps.expanded$motif_levels, levels=c("no_motif_no_enh", "no_motif_non_tiss_match", "no_motif_tiss_match", "motif_no_enh", "motif_non_tiss_match", "motif_tiss_match"), ordered=T)
    ## make another factor so that i can facet on motif ovelrap
    top_coloc_enh_motif_overlaps.expanded$motif_overlap <- factor(ifelse(top_coloc_enh_motif_overlaps.expanded$num_motif_hits > 0, "TFBS", "no TFBS"))

    ## do statistical comparisons of the ABF distributions by enhancer overlap
    ## this is just looking at the enhancer level, without motifs
    enh_levels_with_data <- levels(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels)[table(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels) > 0]

    if(length(enh_levels_with_data) > 1) {
        relevant_class_abf_tests <- do.call(rbind, apply(combn(enh_levels_with_data, m=2),
                                  2, function(x) {
                                      wilcox_res <- with(top_coloc_enh_motif_overlaps.expanded, wilcox.test(coloc_snp_prob[relevant_class_enh_levels==x[1]], coloc_snp_prob[relevant_class_enh_levels==x[2]], alternative="two.sided"))
                                      cond1_size <- sum(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels==x[1])
                                      cond2_size <- sum(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_levels==x[2])

                                      return(data.frame(lev1=x[1], lev2=x[2],
                                                        cond1_size=cond1_size, cond2_size=cond2_size,
                                                        test_stat=wilcox_res[['statistic']],
                                                        p.value=wilcox_res[['p.value']],
                                                        null.value=wilcox_res[['null.value']],
                                                        alternative=wilcox_res[['alternative']],
                                                        method=wilcox_res[['method']]))
                                  }))
        ## write out these statistics
        write.table(relevant_class_abf_tests, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_relevant_class_enh_abf_statistics.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## also write them out as matrices
    write.table(acast(relevant_class_abf_tests, lev1 ~ lev2, value.var="p.value"),
                paste0(outdir, '/tables/', outprefix, '_gtex_coloc_relevant_class_enh_abf_statistics_matrix.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=T, col.names=T)      
    }
    
    ## next do the comparison also considering motif levels
    enh_motif_levels_with_data <- levels(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels)[table(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels) > 0]

    if(length(enh_motif_levels_with_data) > 1) {
        motif_relevant_class_abf_tests <- do.call(rbind,
            apply(combn(enh_motif_levels_with_data, m=2),
                                  2, function(x) {
                                      wilcox_res <- with(top_coloc_enh_motif_overlaps.expanded, wilcox.test(coloc_snp_prob[relevant_class_enh_motif_levels==x[1]], coloc_snp_prob[relevant_class_enh_motif_levels==x[2]], alternative="two.sided"))

                                      cond1_size <- sum(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels==x[1])
                                      cond2_size <- sum(top_coloc_enh_motif_overlaps.expanded$relevant_class_enh_motif_levels==x[2])

                                      return(data.frame(lev1=x[1], lev2=x[2],
                                                        cond1_size=cond1_size, cond2_size=cond2_size,
                                                        test_stat=wilcox_res[['statistic']],
                                                        p.value=wilcox_res[['p.value']],
                                                        null.value=wilcox_res[['null.value']],
                                                        alternative=wilcox_res[['alternative']],
                                                        method=wilcox_res[['method']]))
                                  }))
        
        write.table(motif_relevant_class_abf_tests, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_relevant_class_enh_abf_statistics.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

        write.table(acast(motif_relevant_class_abf_tests, lev1 ~ lev2, value.var="p.value"),
                    paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_relevant_class_enh_abf_statistics_matrix.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=T, col.names=T)

    }
        
    ## just compare having a motif hit or not:
    if(sum(top_coloc_enh_motif_overlaps.expanded$num_motif_hits > 0) > 0 & sum(top_coloc_enh_motif_overlaps.expanded$num_motif_hits==0) > 0) {
        motif_abf_test_result <- with(top_coloc_enh_motif_overlaps.expanded, wilcox.test(coloc_snp_prob[num_motif_hits > 0], coloc_snp_prob[num_motif_hits==0], alternative="two.sided"))
        num_motif_hits <- sum(top_coloc_enh_motif_overlaps.expanded$num_motif_hits > 0)
        num_non_motif_hits <- sum(top_coloc_enh_motif_overlaps.expanded$num_motif_hits==0)
        
        write.table(data.frame(num_motif_hits=num_motif_hits,
                               num_non_motif_hits=num_non_motif_hits,
                               test_stat=motif_abf_test_result[['statistic']],
                               p.value=motif_abf_test_result[['p.value']],
                               null.value=motif_abf_test_result[['null.value']],
                               alternative=motif_abf_test_result[['alternative']],
                               method=motif_abf_test_result[['method']]),
                    paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_relevant_class_motif_vs_no_motif_abf_statistics.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)
    }
        
    ## look at the individual SNP probability distributions by enhancer overlap and tag region
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_motif_overlap_and_tag_region_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_overlaps.expanded, aes(x=motif_levels, y=coloc_snp_prob, color=motif_levels)) +
        geom_boxplot() +
        facet_wrap(~ tag_region, ncol=3, drop=FALSE) + 
        ggtitle(paste("Distributions of SNP ABFs, expanded to", coloc_abf_thresh, "probability, GTEx")) + 
        theme_bw() + xlab("Type of motif and enhancer overlap") + ylab("Individual SNP ABFs") +
        theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, size=15),
              axis.text.y = element_text(size=15), strip.text=element_text(size=25),
              legend.text = element_text(size=15), legend.title=element_text(size=25),
              axis.title=element_text(size=25),
              plot.title=element_text(hjust=0.5, size=20)))
    dev.off()

    ## same thing, but collapsed across tag regions
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_enh_motif_overlap_', coloc_abf_thresh, '_prob_thresh'))
    print(ggplot(top_coloc_enh_motif_overlaps.expanded, aes(x=enh_levels, y=coloc_snp_prob, color=enh_levels)) +
        geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
        ggtitle(paste("Distributions of SNP ABFs across tag regions\nexpanded to", coloc_abf_thresh, "probability, GTEx")) +
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

    ## also make this plot for relevant class overlaps
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_snp_abf_distributions_by_relevant_class_enh_motif_overlap_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_overlaps.expanded, aes(x=relevant_class_enh_levels, y=coloc_snp_prob, color=relevant_class_enh_levels)) +
        geom_boxplot() + geom_jitter(height=0, shape=23, alpha=0.6) +
        ggtitle(paste("Distributions of SNP ABFs across tag regions\nexpanded to", coloc_abf_thresh, "probability, GTEx")) +
          scale_color_manual(values=c("relevant_tiss_match"="red", "irrelevant_tiss_match"="firebrick3", "tiss_non_match"="darkred", "no_overlap"="gray20")) + 
        scale_x_discrete(breaks=c("no_overlap", "tiss_non_match", "irrelevant_tiss_match", "relevant_tiss_match"), labels=c("No enhancer overlap", "Enhancer in inconsistent tissue class", "Enhancer in irrelevant tissue class", paste("Enhancer in", paste0(relevant_classes, collapse=", ")))) + 
        facet_grid(motif_overlap ~ ., scales="free") + 
        theme_bw() + xlab("Type of enhancer overlaps") + ylab("Individual SNP ABFs") +
        scale_y_log10(breaks=c(0, 0.1, 0.25, 0.5, 0.75, 1.0)) + 
        theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=20),
              axis.text.y = element_text(size=20), strip.text=element_text(size=25),
              legend.text = element_text(size=25), legend.title=element_text(size=35),
              axis.title=element_text(size=35),
              plot.title=element_text(hjust=0.5, size=35)))
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
                motif_hits$f5_category_match=="no" & motif_hits$hmm_category_match=="no"),
            ## also get high-ABF variants
            high_abf_variants = sum(motif_hits$coloc_snp_prob >= coloc_abf_thresh),
            high_abf_enh_hit = sum(motif_hits$coloc_snp_prob >= coloc_abf_thresh & motif_hits$any_enh_overlap=="yes"),
            high_abf_non_enh_hit = sum(motif_hits$coloc_snp_prob >= coloc_abf_thresh & motif_hits$any_enh_overlap=="no"),
            high_abf_tiss_match = sum(motif_hits$coloc_snp_prob >= coloc_abf_thresh & (motif_hits$f5_category_match=="yes" | motif_hits$hmm_category_match=="yes")),
            high_abf_tiss_non_match = sum(motif_hits$coloc_snp_prob >= coloc_abf_thresh & motif_hits$any_enh_overlap=="yes" & motif_hits$f5_category_match=="no" & motif_hits$hmm_category_match=="no"))
        
        ## summarize the characteristics of the non-motif hits
        non_motif_hit_summary <- data.frame(
            num_coloc_hits = nrow(non_motif_hits),
            enh_overlap_hits = sum(non_motif_hits$any_enh_overlap=="yes"),
            enh_non_overlap_hits = sum(non_motif_hits$any_enh_overlap=="no"),
            enh_tiss_match = sum(non_motif_hits$f5_category_match=="yes" |
                non_motif_hits$hmm_category_match=="yes"),
            enh_tiss_non_match = sum(non_motif_hits$any_enh_overlap=="yes" &
                non_motif_hits$f5_category_match=="no" & non_motif_hits$hmm_category_match=="no"),
            ## also get high-ABF variants
            high_abf_variants = sum(non_motif_hits$coloc_snp_prob >= coloc_abf_thresh),
            high_abf_enh_hit = sum(non_motif_hits$coloc_snp_prob >= coloc_abf_thresh & non_motif_hits$any_enh_overlap=="yes"),
            high_abf_non_enh_hit = sum(non_motif_hits$coloc_snp_prob >= coloc_abf_thresh & non_motif_hits$any_enh_overlap=="no"),
            high_abf_tiss_match = sum(non_motif_hits$coloc_snp_prob >= coloc_abf_thresh & (non_motif_hits$f5_category_match=="yes" | non_motif_hits$hmm_category_match=="yes")),
            high_abf_tiss_non_match = sum(non_motif_hits$coloc_snp_prob >= coloc_abf_thresh & non_motif_hits$any_enh_overlap=="yes" & non_motif_hits$f5_category_match=="no" & non_motif_hits$hmm_category_match=="no"))

        return(rbind(cbind(hit_motif="motif_overlap", motif_hit_summary),
                     cbind(hit_motif="no_motif_overlap", non_motif_hit_summary)))
    })

    ## add in 0 counts for the regions we did not observe
    if(sum(!(all_regions %in% expanded_tag_region_motif_enh_summary$tag_region)) > 0) {
        expanded_tag_region_motif_enh_summary <- rbind(expanded_tag_region_motif_enh_summary,
                                                       data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0, high_abf_variants=0, high_abf_enh_hit=0, high_abf_non_enh_hit=0, high_abf_tiss_match=0, high_abf_tiss_non_match=0),
                                                       data.frame(tag_region=all_regions[!(all_regions %in% expanded_tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_tiss_match=0, enh_tiss_non_match=0, high_abf_variants=0, high_abf_enh_hit=0, high_abf_non_enh_hit=0, high_abf_tiss_match=0, high_abf_tiss_non_match=0))
    }
        
    ## write this summary table out
    write.table(expanded_tag_region_motif_enh_summary, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_enh_tag_region_summary.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)
    
    ## melt this down
    melted_expanded_enh_motif_summ <- melt(expanded_tag_region_motif_enh_summary, id.vars=c("tag_region", "hit_motif"), measure.vars=c("enh_tiss_non_match", "enh_tiss_match", "enh_non_overlap_hits"), value.name="count")

    melted_expanded_enh_motif_summ$variable <- factor(melted_expanded_enh_motif_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_tiss_match"), ordered=T)

    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_enh_and_motif_overlaps_barplot_', coloc_abf_thresh, '_prob_thresh'), width_ratio=2.0)
    print(ggplot(melted_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_tiss_match"="Enhancer overlaps with consistent tissue class", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_tiss_match"="#B3112E", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_abf_thresh, "probability")) + 
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
    if(sum(!(all_regions %in% expanded_relevant_class_tag_region_motif_enh_summary$tag_region)) > 0) {
        expanded_relevant_class_tag_region_motif_enh_summary <- rbind(expanded_relevant_class_tag_region_motif_enh_summary,
                                                                      data.frame(tag_region=all_regions[!(all_regions %in% expanded_relevant_class_tag_region_motif_enh_summary$tag_region)], hit_motif="motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_relevant_tiss_match=0, enh_irrelevant_tiss_match=0, enh_tiss_non_match=0),
                                                                      data.frame(tag_region=all_regions[!(all_regions %in% expanded_relevant_class_tag_region_motif_enh_summary$tag_region)], hit_motif="no_motif_overlap", num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_relevant_tiss_match=0, enh_irrelevant_tiss_match=0, enh_tiss_non_match=0))
    }
    
    ## write this summary table out
    write.table(expanded_relevant_class_tag_region_motif_enh_summary, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_motif_enh_tag_region_relevant_class_summary.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## melt this down for motif analysis
    melted_relevant_class_expanded_enh_motif_summ <- melt(expanded_relevant_class_tag_region_motif_enh_summary, id.vars=c("tag_region", "hit_motif"), measure.vars=c("enh_tiss_non_match", "enh_relevant_tiss_match", "enh_irrelevant_tiss_match", "enh_non_overlap_hits"), value.name="count")

    melted_relevant_class_expanded_enh_motif_summ$variable <- factor(melted_relevant_class_expanded_enh_motif_summ$variable, levels=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_irrelevant_tiss_match", "enh_relevant_tiss_match"), ordered=T)

    ## make two barplots: one vertical (for paper) and one horizontal (for talks)
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_enh_and_motif_overlaps_barplot_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5)
    print(ggplot(melted_relevant_class_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in irrelevant tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_relevant_tiss_match"="red", "enh_irrelevant_tiss_match"="firebrick3", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_abf_thresh, "probability")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(nrow=4, title="")) +
    facet_wrap(~factor(hit_motif, levels=c("motif_overlap", "no_motif_overlap"),
                       labels=c("Motif Overlap", "No Motif Overlap")),
                       nrow=2, scales="free") + 
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_text(size=25),
          title=element_text(size=15), legend.text=element_text(size=20), 
          plot.title=element_text(hjust=0.5)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_relevant_class_enh_and_motif_overlaps_horizontal_barplot_', coloc_abf_thresh, '_prob_thresh'), width_ratio=1.5)
    print(ggplot(melted_relevant_class_expanded_enh_motif_summ, aes(x=tag_region, y=count, fill=variable)) +
        scale_fill_manual(
            labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in irrelevant tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
            values=c("enh_relevant_tiss_match"="red", "enh_irrelevant_tiss_match"="firebrick3", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of SNP - eQTL comparisons") + theme_bw() +
    ggtitle(paste("Enhancer and motif analysis of GTEx-colocalized SNPs expanded to", coloc_abf_thresh, "probability")) + 
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

    ## now make a master unified plot that stratifies by both motif overlap and ABF level
    ## to do this, make a new summary data frame split by motif and ABF overlaps
    motif_and_abf_summary_df <- ddply(top_coloc_enh_motif_overlaps.expanded, .(tag_region), function(x) {
        ## define a bunch of logical vectors that we use to count various types of variants
        motif_hits <- x$num_motif_hits > 0
        non_motif_hits <- x$num_motif_hits==0
        high_abf <- x$coloc_snp_prob >= coloc_abf_thresh
        low_abf <- x$coloc_snp_prob < coloc_abf_thresh        
        enh_overlap <- x$any_enh_overlap=="yes"
        enh_non_overlap <- x$any_enh_overlap=="no"
        relevant_tiss_match <- (x$f5_category_match=="yes" | x$hmm_category_match=="yes") & x$gtex_tissue_class %in% relevant_classes
        irrelevant_tiss_match <- (x$f5_category_match=="yes" | x$hmm_category_match=="yes") & !(x$gtex_tissue_class %in% relevant_classes)
        tiss_non_match <- x$any_enh_overlap=="yes" & x$f5_category_match=="no" & x$hmm_category_match=="no"
        
        ## summarize the characteristics of the four combinations
        motif_hit_high_abf <- data.frame(
            motif="TFBS overlap", ABF="High ABF",
            num_coloc_hits = sum(motif_hits & high_abf),
            enh_overlap_hits = sum(motif_hits & high_abf & enh_overlap),
            enh_non_overlap_hits = sum(motif_hits & high_abf & enh_non_overlap),
            enh_relevant_tiss_match = sum(motif_hits & high_abf & relevant_tiss_match),
            enh_irrelevant_tiss_match = sum(motif_hits & high_abf & irrelevant_tiss_match),
            enh_tiss_non_match = sum(motif_hits & high_abf & tiss_non_match))

        motif_hit_low_abf <- data.frame(
            motif="TFBS overlap", ABF="Low ABF",
            num_coloc_hits = sum(motif_hits & low_abf),
            enh_overlap_hits = sum(motif_hits & low_abf & enh_overlap),
            enh_non_overlap_hits = sum(motif_hits & low_abf & enh_non_overlap),
            enh_relevant_tiss_match = sum(motif_hits & low_abf & relevant_tiss_match),
            enh_irrelevant_tiss_match = sum(motif_hits & low_abf & irrelevant_tiss_match),
            enh_tiss_non_match = sum(motif_hits & low_abf & tiss_non_match))

        non_motif_hit_high_abf <- data.frame(
            motif="No TFBS overlap", ABF="High ABF",
            num_coloc_hits = sum(non_motif_hits & high_abf),
            enh_overlap_hits = sum(non_motif_hits & high_abf & enh_overlap),
            enh_non_overlap_hits = sum(non_motif_hits & high_abf & enh_non_overlap),
            enh_relevant_tiss_match = sum(non_motif_hits & high_abf & relevant_tiss_match),
            enh_irrelevant_tiss_match = sum(non_motif_hits & high_abf & irrelevant_tiss_match),
            enh_tiss_non_match = sum(non_motif_hits & high_abf & tiss_non_match))

        non_motif_hit_low_abf <- data.frame(
            motif="No TFBS overlap", ABF="Low ABF",
            num_coloc_hits = sum(non_motif_hits & low_abf),
            enh_overlap_hits = sum(non_motif_hits & low_abf & enh_overlap),
            enh_non_overlap_hits = sum(non_motif_hits & low_abf & enh_non_overlap),
            enh_relevant_tiss_match = sum(non_motif_hits & low_abf & relevant_tiss_match),
            enh_irrelevant_tiss_match = sum(non_motif_hits & low_abf & irrelevant_tiss_match),
            enh_tiss_non_match = sum(non_motif_hits & low_abf & tiss_non_match))        
        
        return(rbind(motif_hit_high_abf, motif_hit_low_abf, non_motif_hit_high_abf, non_motif_hit_low_abf))
    })

    ## add in 0 counts for the regions we did not observe
    if(sum(!(all_regions %in% motif_and_abf_summary_df$tag_region)) > 0) {
        motif_and_abf_summary_df <- rbind(motif_and_abf_summary_df,
                                          data.frame(tag_region=rep(all_regions[!(all_regions %in% motif_and_abf_summary_df$tag_region)], each=4), motif=c("TFBS overlap", "TFBS overlap", "No TFBS overlap", "No TFBS overlap"), ABF=c("High ABF", "Low ABF", "High ABF", "Low ABF"), num_coloc_hits=0, enh_overlap_hits=0, enh_non_overlap_hits=0, enh_relevant_tiss_match=0, enh_irrelevant_tiss_match=0, enh_tiss_non_match=0))
    }

    ## write this summary table out
    write.table(motif_and_abf_summary_df, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_ABF_and_motif_summary.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## now analyze the enrichment patterns
    melted_motif_and_abf_summary <- melt(motif_and_abf_summary_df, id.vars=c("tag_region", "motif", "ABF"), measure.vars=c("enh_non_overlap_hits", "enh_tiss_non_match", "enh_irrelevant_tiss_match", "enh_relevant_tiss_match"), value.name="count")

    melted_motif_and_abf_summary$motif <- factor(melted_motif_and_abf_summary$motif, levels=c("TFBS overlap", "No TFBS overlap"), ordered=T)
    
    ## make the master barplot
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_ABF_and_motif_summary_barplot_', coloc_abf_thresh, '_prob_thresh'), width_ratio=1.75, height_ratio=1.5)
    print(ggplot(melted_motif_and_abf_summary, aes(x=tag_region, y=count, fill=variable)) +
          scale_fill_manual(
              labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in irrelevant tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
              values=c("enh_relevant_tiss_match"="cyan3", "enh_irrelevant_tiss_match"="darkblue", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of variants") + theme_bw() +
    ggtitle(paste("Enhancer overlap analysis of GTEx-colocalized SNPs\n expanded to", coloc_abf_thresh, "probability")) + 
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(ncol=2, title="")) +
    facet_wrap(~ motif + ABF, scales="free", nrow=2, ncol=2) + 
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=60, hjust=1, size=20*text_mult),
          axis.text.y = element_text(size=20*text_mult),
          strip.text=element_text(size=25*text_mult),
          title=element_text(size=35*text_mult),
          legend.text=element_text(size=25*text_mult),
          plot.title=element_text(hjust=0.5)))          
    dev.off()        

    ## no title
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_ABF_and_motif_summary_barplot_', coloc_abf_thresh, '_prob_thresh_no_title'), width_ratio=1.75, height_ratio=1.5)
    print(ggplot(melted_motif_and_abf_summary, aes(x=tag_region, y=count, fill=variable)) +
          scale_fill_manual(
              labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in irrelevant tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
              values=c("enh_relevant_tiss_match"="cyan3", "enh_irrelevant_tiss_match"="darkblue", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
    xlab("Tag region") + ylab("Number of variants") + theme_bw() +
    geom_bar(stat="identity", position="stack") +
    guides(fill=guide_legend(ncol=2, title="")) +
    facet_wrap(~ motif + ABF, scales="free", nrow=2, ncol=2) + 
    theme(legend.position="bottom",
          axis.text.x=element_text(angle=60, hjust=1, size=20*text_mult),
          axis.text.y = element_text(size=20*text_mult),
          strip.text=element_text(size=25*text_mult),
          title=element_text(size=35*text_mult),
          legend.text=element_text(size=25*text_mult)))
    dev.off()        

    ## now make barplots that only include the tag regions with hits
    nonzero_motif_and_abf_summary <- melted_motif_and_abf_summary[melted_motif_and_abf_summary$count > 0,]
    ## only make it if we have any hits
    if(nrow(nonzero_motif_and_abf_summary) > 0) {
        make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_ABF_and_motif_nonzero_summary_barplot_', coloc_abf_thresh, '_prob_thresh_no_title'), width_ratio=1.75, height_ratio=1.5)
        print(ggplot(nonzero_motif_and_abf_summary, aes(x=tag_region, y=count, fill=variable)) +
              scale_fill_manual(
                  labels=c("enh_relevant_tiss_match"=paste("Enhancer overlaps in", paste0(relevant_classes, collapse=", ")), "enh_irrelevant_tiss_match"="Enhancer overlaps in irrelevant tissue classes", "enh_tiss_non_match"="Enhancer overlaps without consistent tissue class", "enh_non_overlap_hits"="Colocalization hits without enhancer overlap"),
                  values=c("enh_relevant_tiss_match"="cyan3", "enh_irrelevant_tiss_match"="darkblue", "enh_tiss_non_match"="darkred", "enh_non_overlap_hits"="gray20")) + 
              xlab("Tag region") + ylab("Number of variants") + theme_bw() +
              geom_bar(stat="identity", position="stack") +
              guides(fill=guide_legend(ncol=2, title="")) +
              facet_wrap(~ motif + ABF, scales="free", nrow=2, ncol=2) + 
              theme(legend.position="bottom",
                    axis.text.x=element_text(angle=60, hjust=1, size=20*text_mult),
                    axis.text.y = element_text(size=20*text_mult),
                    strip.text=element_text(size=25*text_mult),
                    title=element_text(size=35*text_mult),
                    legend.text=element_text(size=25*text_mult)))
        dev.off()                
    }    
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
    cat("Comparing GTEx direct eQTL overlaps with HOMER motif disruptions\n")
    eqtl_overlap_file <- paste0(param_ref[['outdir']], '/analysis_results/gtex_eqtl_overlap/tables/',
                                param_ref[['outprefix']], "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh,
                                "_uniq_snp_eqtl_overlaps_no_tagsnp_info.txt")
    
    uniq_snp_eqtl_overlaps <- read.table(eqtl_overlap_file, header=T, sep="\t", quote="", as.is=T)

    if(nrow(uniq_snp_eqtl_overlaps)==0) {
        cat("No direct eQTL overlaps found!\n")
    } else {
    
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
    write.table(top_coloc_enh_motif_eqtl_overlaps, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_motif_eqtl_info.txt'), quote=F, sep="\t", row.names=F)
    
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
            cat("There were multiple hits for SNP", as.character(this_data$high_coloc_snp), "tissue", as.character(this_data$tissue), "and gene", as.character(this_data$eqtl_gene_name), "\n")
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
    write.table(top_coloc_enh_motif_eqtl_overlaps.expanded, paste0(outdir, '/tables/', outprefix, '_gtex_coloc_enh_motif_eqtl_info.', coloc_abf_thresh, '_thresh_expanded.txt'), quote=F, sep="\t", row.names=F)

    ## check whether the eQTLs in each comparison line up or not:
    expanded_direct_eqtl_info <- ddply(top_coloc_enh_motif_eqtl_overlaps.expanded, .(tag_region, tissue, gtex_tissue_class, eqtl_gene_name), summarize, eqtl_overlaps = sum(!is.na(beta)), eqtl_overlap_prop = sum(!is.na(beta)) / length(beta), positive_betas = sum(beta > 0, na.rm=T), negative_betas = sum(beta < 0, na.rm=T))

    cat(sum(expanded_direct_eqtl_info$eqtl_overlaps!=0), "comparisons have eQTL overlap and", sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas > 0 & negative_betas > 0)),
        "of these have variants with different beta directions\n")
    cat("Of the rest,", sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas > 0 & negative_betas==0)), "are all positive and",
        sum(with(expanded_direct_eqtl_info[expanded_direct_eqtl_info$eqtl_overlaps!=0,], positive_betas==0 & negative_betas > 0)), "are all negative\n")

    ## make a plot of the beta distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_beta_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$beta),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=beta, color=beta)) +
    xlab("Tissue; target gene; tag region") + ylab("Beta value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap beta distributions, expanded to", coloc_abf_thresh, "probability")) + 
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
    
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_t_stat_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$t_stat),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=t_stat, color=t_stat)) +
    xlab("Tissue; target gene; tag region") + ylab("T statistics value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap T-statistic distributions, expanded to", coloc_abf_thresh, "probability")) +
    scale_y_continuous(breaks=seq(tstat_min, tstat_max, by=1)) + 
    geom_point() + geom_hline(color="red", yintercept=0, linetype=3) + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## make a plot of the se distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_se_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$se),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=se, color=se)) +
    xlab("Tissue; target gene; tag region") + ylab("Standard error value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap SE distributions, expanded to", coloc_abf_thresh, "probability")) + 
    geom_point() + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()

    ## make a plot of the p value distributions across comparisons
    make_graphic(paste0(outdir, '/plots/', outprefix, '_gtex_coloc_direct_eqtl_p_value_dists_', coloc_abf_thresh, '_prob_thresh'), height_ratio=1.5, width_ratio=1.5)
    print(ggplot(top_coloc_enh_motif_eqtl_overlaps.expanded[!is.na(top_coloc_enh_motif_eqtl_overlaps.expanded$p_value),],
                 aes(x=paste(tissue, eqtl_gene_name, tag_region, sep=";"),
                     y=p_value, color=p_value)) +
    xlab("Tissue; target gene; tag region") + ylab("P-value from direct eQTL overlap") + theme_bw() +
    ggtitle(paste("Direct eQTL overlap P value distributions, expanded to", coloc_abf_thresh, "probability")) + 
    geom_point() + coord_flip() + 
    facet_grid(tag_region ~ ., scales="free", space="free") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.text.y = element_text(size=15), strip.text=element_blank(),
          title=element_text(size=15), 
          plot.title=element_text(size=20, hjust=0.5)))
    dev.off()    

    }
}
