## all_gene_GTEx_correlations.R
## alex amlie-wolf 06/13/18
## analyzing correlation distribution of all genes in the genome

library(ggplot2)
library(reshape2)
library(plyr)
library(psych)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------

## 0. Table of Contents
## 1. Function Definitions
## 2. Set up data paths
## 3. Read in expression data, calculate correlations

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

## -----------------------------------------------------------------------------
## 2. Set up data paths
## -----------------------------------------------------------------------------
gtex_expr_dir <- "/home/alexaml/data/GTEx/GTEx_v6p/rnaseq/all_chr_gtex_v6p_expression/"
sample_info_file <- "/home/alexaml/data/GTEx/subject_sample_info/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
outdir <- "~/data/GTEx/GTEx_v6p/all_gene_correlation_analysis/"
dir.create(paste0(outdir, "/plots/"), F, T)
dir.create(paste0(outdir, "/tables/"), F, T)
dir.create(paste0(outdir, "/full_correlation_tables/"), F, T)

## -----------------------------------------------------------------------------
## 3. Read in expression data, calculate correlations
## -----------------------------------------------------------------------------
{
    expr_reading_start <- proc.time()
    ## just make a huge matrix of all the gene expression
    ## rows are samples, columns are genes
    ## hard-coded sample number for GTEx v6p
    all_chr_gene_expression <- matrix(data=NA, nrow=8555, ncol=0)

    for(chr_file in sort(list.files(gtex_expr_dir, 'chr.*_gtex_expression.txt', full.names=T))) {    
        this_chr <- strsplit(basename(chr_file), "_")[[1]][1]
        cat("Reading gene expression data from", this_chr, "\n")
        this_chr_expression <- read.table(chr_file, header=T, sep="\t", as.is=T)
        
        ## make this into a matrix of gene expression levels (genes are columns, samples are rows)
        gene_expression_mat <- t(data.matrix(this_chr_expression[,3:ncol(this_chr_expression)]))
        colnames(gene_expression_mat) <- this_chr_expression$Description

        ## filter out genes that are never expressed
        nonexpressed_genes <- apply(gene_expression_mat, 2, function(x) {all(x==0)})
        cat(sum(nonexpressed_genes), "genes are never expressed in", this_chr, "and are filtered out.\n")

        ## remove duplicated genes
        gene_expression_mat <- gene_expression_mat[,!nonexpressed_genes & !duplicated(colnames(gene_expression_mat))]

        ## add this to the full matrix
        all_chr_gene_expression <- cbind(all_chr_gene_expression, gene_expression_mat)
    }

    cat("Reading in all expression data took", (proc.time() - expr_reading_start)[['elapsed']], 'seconds\n')
}

{
    pearson_start <- proc.time()
    ## pearson_cor <- corr.test(all_chr_gene_expression, method="pearson", ci=FALSE)
    pearson_cor <- cor(all_chr_gene_expression, method="pearson")
    cat("Computing Pearson took", (proc.time() - pearson_start)[['elapsed']], 'seconds\n')

    spearman_start <- proc.time()
    ## spearman_cor <- corr.test(all_chr_gene_expression, method="spearman", ci=FALSE)
    spearman_cor <- cor(all_chr_gene_expression, method="spearman")
    cat("Computing Spearman took", (proc.time() - spearman_start)[['elapsed']], 'seconds\n')

}
