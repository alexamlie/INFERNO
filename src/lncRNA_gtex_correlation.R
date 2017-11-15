## lncRNA_gtex_correlation.R
## alex amlie-wolf 06/16/17
## after colocalization analysis, this script takes any highly colocalized lncRNAs and finds
## their targets by using expression correlation across all genes and samples in GTEx

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
## 2. Read in data
## 3. Identify colocalized lncRNAs, pull out their expression vectors
## 4. Compute correlation of all lncRNAs with all GTEx genes
## 5. Analyze correlation patterns and GTEx expression of individual lncRNAs
## 6. Summarize correlation patterns across all lncRNAs

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
## 2. Read in data
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
cat(args, '\n')

if(length(args)==10) {
    outdir <- args[1]
    coloc_result_file <- args[2]   
    gtex_expr_dir <- args[3] 
    sample_info_file <- args[4]
    gencode_lncrna_file <- args[5]    
    fantom5_class_file <- args[6]
    gtex_class_file <- args[7]
    roadmap_class_file <- args[8]
    coloc_h4_thresh <- as.numeric(args[9])
    cor_thresh <- as.numeric(args[10])
} else {
    stop("Incorrect number of arguments. Should be: <outdir> <coloc_result_file> <gtex_expression_dir> <gtex_sample_info> <gencode lncRNA annotation table> <fantom5_class_file> <gtex_class_file> <roadmap_class_file> <P(H_4) threshold> <correlation value threshold>\n")
}

dir.create(paste0(outdir, "/plots/"), F, T)
dir.create(paste0(outdir, "/tables/"), F, T)
dir.create(paste0(outdir, "/full_correlation_tables/"), F, T)

summary_file <- paste0(outdir, '/tables/', 'lncRNA_corr_summary.txt')

## read in the lncRNA annotations (genes only)
gencode_lncrnas <- read.table(pipe(paste0("awk -F$'\t' '$3==\"gene\"' ", gencode_lncrna_file)), header=F, sep="\t", quote="", as.is=T)
## parse out the gene_id information in a new column
all_gencode_info <- strsplit(gsub("\"", "", gencode_lncrnas$V9), "; ")
gencode_lncrnas$gene_id <- unlist(lapply(all_gencode_info, function(x) {
    gsub("gene_id ", "", x[grep("gene_id", x)])
}))

## to get the right colors, read in all the categorizations from INFERNO
fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
## tweak the GTEx categories to match up with the sample data frame
gtex_category_df$SMTSD <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data)
gtex_category_df$SMTSD <- gsub("_", " ", gtex_category_df$SMTSD)

roadmap_category_df <- read.table(roadmap_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")

## find the full list of classes
all_classes <- sort(union(fantom5_category_df$Class, union(gtex_category_df$Class, roadmap_category_df$Class)))

## now define the color palette
## generated from http://tools.medialab.sciences-po.fr/iwanthue/
## parameters: H 0-360, C 0.46 - 3, L 0.5-1.5
category_colors <- c(rgb(227,191,35, maxColorValue=255), rgb(224,98,247, maxColorValue=255), rgb(36,153,165, maxColorValue=255), rgb(232,2,58, maxColorValue=255), rgb(30,181,54, maxColorValue=255), rgb(253,164,185, maxColorValue=255), rgb(44,115,232, maxColorValue=255), rgb(83,118,43, maxColorValue=255), rgb(203,50,133, maxColorValue=255), rgb(194,196,252, maxColorValue=255), rgb(157,236,192, maxColorValue=255), rgb(141,86,174, maxColorValue=255), rgb(253,133,114, maxColorValue=255), rgb(175,242,253, maxColorValue=255), rgb(149,87,122, maxColorValue=255), rgb(131,233,24, maxColorValue=255), rgb(182,241,140, maxColorValue=255), rgb(252,50,7, maxColorValue=255), rgb(244,148,221, maxColorValue=255), rgb(40,169,123, maxColorValue=255), rgb(247,183,144, maxColorValue=255), rgb(242,184,94, maxColorValue=255), rgb(53,108,145, maxColorValue=255), rgb(198,8,78, maxColorValue=255), rgb(61,115,80, maxColorValue=255), rgb(41,232,215, maxColorValue=255), rgb(122,107,24, maxColorValue=255), rgb(79,153,241, maxColorValue=255), rgb(100,130,128, maxColorValue=255), rgb(166,169,37, maxColorValue=255), rgb(203,137,237, maxColorValue=255), rgb(178,204,231, maxColorValue=255))
names(category_colors) <- all_classes

cat_col_scale <- scale_fill_manual(name="Tissue Category", values=category_colors)

## for the tissue-specific analysis, we also want to get the sample attributes
sample_info_df <- read.table(sample_info_file, header=T, sep="\t", quote="", as.is=T)
## also add a new column to the sample info DF to match (get rid of the stupid naming in the
## GTEx data)
sample_info_df$inferno_match <- gsub("- ", "", gsub("\\(|\\)", "", sample_info_df$SMTSD))
## and a column for the actual class
sample_info_df$inferno_class <- gtex_category_df$Class[match(sample_info_df$inferno_match, gtex_category_df$SMTSD)]

## get the sample IDs of all samples that we actually have an INFERNO category for
expression_samples <- gsub("-", ".", sample_info_df$SAMPID[sample_info_df$inferno_match %in% gtex_category_df$SMTSD])

## read in the colocalization results
all_coloc_results <- read.table(coloc_result_file, header=T, sep="\t", quote="", as.is=T)
## also figure out the prefix from this
outprefix <- gsub("_gtex_coloc_summaries.txt", "", basename(coloc_result_file), fixed=T)

## -----------------------------------------------------------------------------
## 3. Identify colocalized lncRNAs, pull out their expression vectors
## -----------------------------------------------------------------------------
top_coloc_hits <- all_coloc_results[all_coloc_results$PP.H4.abf >= coloc_h4_thresh,]
rm(all_coloc_results)

## find lncRNAs
## skip the information column from GENCODE, for now
top_lncrna_hits <- merge(top_coloc_hits, gencode_lncrnas[,-9], by.x="eqtl_gene_id", by.y="gene_id")

if(nrow(top_lncrna_hits)==0) {
    stop("No strongly colocalized lncRNA signals found!\n")
} else {
    cat(length(unique(top_lncrna_hits$eqtl_gene_name)), "unique lncRNAs across",
        length(unique(top_lncrna_hits$tissue)), "unique GTEx tissues and",
        length(unique(top_lncrna_hits$gtex_tissue_class)), "tissue classes found in data\n")
}
cat(length(unique(top_lncrna_hits$eqtl_gene_name)), "unique lncRNAs across",
                 length(unique(top_lncrna_hits$tissue)), "unique GTEx tissues and",
                 length(unique(top_lncrna_hits$gtex_tissue_class)), "tissue classes found in data\n", file=summary_file, append=T)

## save the expression vectors as a named list
{
expr_vec_start <- proc.time()
lncrna_expression_vectors <- list()

for(this_chr in sort(unique(top_lncrna_hits$V1))) {
    cat("Reading in lncRNA expression levels from", this_chr, "\n")

    chr_file <- paste0(gtex_expr_dir, '/', this_chr, '_gtex_expression.txt')

    this_chr_expression <- read.table(chr_file, header=T, sep="\t", as.is=T)
    
    ## make this into a matrix of gene expression levels (genes are columns, samples are rows)
    gene_expression_mat <- t(data.matrix(this_chr_expression[,3:ncol(this_chr_expression)]))
    colnames(gene_expression_mat) <- this_chr_expression$Description
    
    ## filter out genes that are never expressed
    nonexpressed_genes <- apply(gene_expression_mat, 2, function(x) {all(x==0)})

    gene_expression_mat <- gene_expression_mat[,!nonexpressed_genes & !duplicated(colnames(gene_expression_mat))]

    for(lncrna in unique(top_lncrna_hits$eqtl_gene_name[top_lncrna_hits$V1==this_chr])) {
        if(lncrna %in% colnames(gene_expression_mat)) {
            lncrna_expression_vectors[[lncrna]] <- as.matrix(gene_expression_mat[,which(colnames(gene_expression_mat)==lncrna)])
        } else {
            cat("lncRNA", lncrna, "was not found in GTEx expression data by name!\n")
        }
    }
}

cat("Reading in all lncRNA expression vectors took", (proc.time() - expr_vec_start)[['elapsed']], 'seconds\n')
}

## -----------------------------------------------------------------------------
## 4. Compute correlation of all lncRNAs with all GTEx genes
## -----------------------------------------------------------------------------
{
corr_start <- proc.time()
## make another list to store the correlation results for each lncRNA
lncrna_correlation_dfs <- list()

for(chr_file in sort(list.files(gtex_expr_dir, 'chr.*_gtex_expression.txt', full.names=T))) {    
    this_chr <- strsplit(basename(chr_file), "_")[[1]][1]
    cat("Reading gene expression data from", this_chr, "\n")
    this_chr_expression <- read.table(chr_file, header=T, sep="\t", as.is=T)
    
    ## make this into a matrix of gene expression levels (genes are columns, samples are rows)
    gene_expression_mat <- t(data.matrix(this_chr_expression[,3:ncol(this_chr_expression)]))
    colnames(gene_expression_mat) <- this_chr_expression$Description
    
    ## filter out genes that are never expressed
    nonexpressed_genes <- apply(gene_expression_mat, 2, function(x) {all(x==0)})
    cat(sum(nonexpressed_genes), "genes are never expressed and are filtered out.\n")

    gene_expression_mat <- gene_expression_mat[,!nonexpressed_genes & !duplicated(colnames(gene_expression_mat))]

    ## for partial correlation analysis against the tissues of each sample
    sample_tissues <- sample_info_df$SMTS[match(rownames(gene_expression_mat),
                                                gsub("-", ".", sample_info_df$SAMPID, fixed=T))]

    ## now go through and compute the correlation of each lncRNA with the genes on this chromosome
    cat("Analyzing lncRNA correlations: ")
    for(lncrna in names(lncrna_expression_vectors)) {
        cat(lncrna, " ")
        this_lncrna_pearson_cor <- corr.test(lncrna_expression_vectors[[lncrna]], gene_expression_mat, method="pearson", adjust="none")
        this_lncrna_spearman_cor <- corr.test(lncrna_expression_vectors[[lncrna]], gene_expression_mat, method="spearman", adjust="none")

        ## extract the residuals of the regression of each gene's expression against these tissues
        this_lncrna_resids <- residuals(lm(lncrna_expression_vectors[[lncrna]] ~ sample_tissues))
        gene_resids <- residuals(lm(gene_expression_mat ~ sample_tissues))
        
        ## now do correlation with these residuals (skip the p-values for these)
        this_lncrna_pearson_partial_cor <- cor(this_lncrna_resids, gene_resids, method="pearson")
        this_lncrna_spearman_partial_cor <- cor(this_lncrna_resids, gene_resids, method="spearman")

        this_lncrna_corr_df <- data.frame(chr=this_chr,
                                          gene=colnames(this_lncrna_pearson_cor$r),
                                          pearson_cor=this_lncrna_pearson_cor$r[1,],
                                          spearman_cor=this_lncrna_spearman_cor$r[1,],
                                          pearson_pval=this_lncrna_pearson_cor$p[1,],
                                          spearman_pval=this_lncrna_spearman_cor$p[1,],
                                          pearson_partial_cor=this_lncrna_pearson_partial_cor[1,],
                                          spearman_partial_cor=this_lncrna_spearman_partial_cor[1,],
                                          stringsAsFactors = F)

        ## add this to the data frame for this lncRNA
        lncrna_correlation_dfs[[lncrna]] <- rbind(lncrna_correlation_dfs[[lncrna]], this_lncrna_corr_df)
        
    }
    cat("\n")
}

cat("Computing all expression correlations took", (proc.time() - corr_start)[['elapsed']], 'seconds\n')
}

## -----------------------------------------------------------------------------
## 5. Analyze correlation patterns and GTEx expression of individual lncRNAs
## -----------------------------------------------------------------------------
## now loop through all the lncRNAs and do individual analysis on them
## for convenience, make a data frame storing the correlated genes and which lncRNA they came from
all_lncrna_targets <- data.frame(stringsAsFactors = FALSE)

for(lncrna in names(lncrna_correlation_dfs)) {
    cat("Analyzing", lncrna, "correlation patterns\n")    
    cat("Analyzing", lncrna, "correlation patterns\n", file=summary_file, append=T)    
    
    ## add mean correlation values
    lncrna_correlation_dfs[[lncrna]]$mean_cor <- rowMeans(lncrna_correlation_dfs[[lncrna]][,c("pearson_cor", "spearman_cor")])
    lncrna_correlation_dfs[[lncrna]]$mean_partial_cor <- rowMeans(lncrna_correlation_dfs[[lncrna]][,c("pearson_partial_cor", "spearman_partial_cor")])
    ## sort by that
    lncrna_correlation_dfs[[lncrna]] <- lncrna_correlation_dfs[[lncrna]][order(abs(lncrna_correlation_dfs[[lncrna]]$mean_cor), decreasing=T),]

    ## correct the p-values
    lncrna_correlation_dfs[[lncrna]]$pearson_padj <- p.adjust(lncrna_correlation_dfs[[lncrna]]$pearson_pval, method="bonferroni")
    lncrna_correlation_dfs[[lncrna]]$spearman_padj <- p.adjust(lncrna_correlation_dfs[[lncrna]]$spearman_pval, method="bonferroni")

    ## write out the table with full correlations, for posterity
    write.table(lncrna_correlation_dfs[[lncrna]], 
                paste0(outdir, '/full_correlation_tables/', lncrna, '_correlations.txt'),
                quote=F, sep="\t", row.names=F, col.names=T)

    ## now get rid of the entry for the lncRNA itself
    lncrna_correlation_dfs[[lncrna]] <- lncrna_correlation_dfs[[lncrna]][lncrna_correlation_dfs[[lncrna]]$gene!=lncrna,]   
    
    ## also write out all the genes meeting a 0.5 correlation threshold
    corr_thresh_vec <- abs(lncrna_correlation_dfs[[lncrna]]$pearson_cor) > cor_thresh & abs(lncrna_correlation_dfs[[lncrna]]$spearman_cor) > cor_thresh
    cat(sum(corr_thresh_vec), 'genes met a', cor_thresh, 'correlation threshold\n')
    cat(sum(corr_thresh_vec), 'genes met a', cor_thresh, 'correlation threshold\n', file=summary_file, append=T)
    
    if(sum(corr_thresh_vec) > 0) {
        write.table(lncrna_correlation_dfs[[lncrna]][corr_thresh_vec,],
                    paste0(outdir, '/tables/', lncrna, '_correlated_genes_over_', cor_thresh, '_spearman_pearson.txt'),
                    quote=F, sep="\t", row.names=F, col.names=T)

        ## save this to the master data frame with only the main correlation numbers
        ## also save the tissues
        all_lncrna_targets <- rbind(all_lncrna_targets,
                                    data.frame(
                                        lncRNA = lncrna,
                                        tissues = paste0(unique(top_lncrna_hits$tissue[top_lncrna_hits$eqtl_gene_name==lncrna]), collapse=";"),
                                        classes = paste0(unique(top_lncrna_hits$gtex_tissue_class[top_lncrna_hits$eqtl_gene_name==lncrna]), collapse=";"),
                                        target_gene = lncrna_correlation_dfs[[lncrna]]$gene[corr_thresh_vec],
                                        lncrna_correlation_dfs[[lncrna]][corr_thresh_vec, c("pearson_cor", "spearman_cor")], stringsAsFactors = F))
    }

    lncrna_out <- gsub("/", "_", lncrna, fixed=T)
    dir.create(paste0(outdir, '/plots/', lncrna_out), F, T)

    ## next, make a scatterplot to compare the two types of correlations
    heatmap_col_lims <- c("blue", "red")

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_correlation_scatterplot'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=pearson_cor, y=spearman_cor)) +
          geom_point(aes(color=mean_cor)) + 
          scale_colour_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2],
                                guide="colorbar", limits=range(lncrna_correlation_dfs[[lncrna]]$mean_cor)) +
          theme_bw() +
          scale_x_continuous(limits=range(lncrna_correlation_dfs[[lncrna]]$pearson_cor)*1.1, expand = c(0.05, 0),
                             breaks=seq(-1, 1, by=0.1)) +
          scale_y_continuous(limits=range(lncrna_correlation_dfs[[lncrna]]$spearman_cor)*1.1, expand = c(0.05, 0),
                             breaks=seq(-1, 1, by=0.1)) +
          geom_hline(yintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_hline(yintercept=-cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +
          xlab("Pearson correlation") + ylab("Spearman correlation") +
          ggtitle(paste("Correlation values of GTEx genes with", lncrna)) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=20),
                plot.title = element_text(hjust = 0.5)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_partial_correlation_scatterplot'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=pearson_partial_cor, y=spearman_partial_cor)) +
          geom_point(aes(color=mean_partial_cor)) + 
          scale_colour_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2],
                                guide="colorbar", limits=range(lncrna_correlation_dfs[[lncrna]]$mean_partial_cor)) +
          theme_bw() +
          scale_x_continuous(limits=range(lncrna_correlation_dfs[[lncrna]]$pearson_partial_cor), expand = c(0.05, 0),
                             breaks=seq(-1, 1, by=0.1)) +
          scale_y_continuous(limits=range(lncrna_correlation_dfs[[lncrna]]$spearman_partial_cor), expand = c(0.05, 0),
                             breaks=seq(-1, 1, by=0.1)) +
          geom_hline(yintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_hline(yintercept=-cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +          
          xlab("Pearson partial correlation") + ylab("Spearman partial correlation") +
          ggtitle(paste("Partial correlation values of all genes with", lncrna, "\nAdjusted for tissues")) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=20)))
    dev.off()

    ## now make histograms for each type of correlation value
    hist_binwidth <- 0.01

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_pearson_cor_histogram'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=pearson_cor, colour=pearson_cor)) +
          geom_histogram(colour="black", fill="lightblue", binwidth=hist_binwidth) +
          scale_y_continuous(trans="log1p", breaks=c(1, 10, 100, 500, 1000, 2000, 4000, 6000)) + 
          scale_x_continuous(breaks=seq(-1, 1, by=0.1)) + 
          theme_bw() +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +          
          xlab("Pearson correlation") + ylab("Number of genes") + 
          ggtitle(paste("Histogram of Pearson correlation of all genes with", lncrna)) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=20)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_pearson_partial_cor_histogram'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=pearson_partial_cor, colour=pearson_partial_cor)) +
          geom_histogram(colour="black", fill="lightblue", binwidth=hist_binwidth) +
          scale_y_continuous(trans="log1p", breaks=c(1, 10, 100, 500, 1000, 2000, 4000, 6000)) + 
          scale_x_continuous(breaks=seq(-1, 1, by=0.1)) + 
          theme_bw() +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +                    
          xlab("Pearson partial correlation") + ylab("Number of genes") + 
          ggtitle(paste("Histogram of Pearson partial correlation of all genes with", lncrna, "\nAdjusted for tissue")) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=20)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_spearman_cor_histogram'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=spearman_cor, colour=spearman_cor)) +
          geom_histogram(colour="black", fill="lightblue", binwidth=hist_binwidth) +
          scale_y_continuous(trans="log1p", breaks=c(1, 10, 100, 500, 1000, 2000, 4000, 6000)) + 
          scale_x_continuous(breaks=seq(-1, 1, by=0.1)) + 
          theme_bw() +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +                    
          xlab("Spearman correlation") + ylab("Number of genes") + 
          ggtitle(paste("Histogram of Spearman correlation of all genes with", lncrna)) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=19)))
    dev.off()

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_spearman_partial_cor_histogram'))
    print(ggplot(lncrna_correlation_dfs[[lncrna]], aes(x=spearman_partial_cor, colour=spearman_partial_cor)) +
          geom_histogram(colour="black", fill="lightblue", binwidth=hist_binwidth) +
          scale_y_continuous(trans="log1p", breaks=c(1, 10, 100, 500, 1000, 2000, 4000, 6000)) + 
          scale_x_continuous(breaks=seq(-1, 1, by=0.1)) + 
          theme_bw() +
          geom_vline(xintercept=cor_thresh, linetype=3) +
          geom_vline(xintercept=-cor_thresh, linetype=3) +                    
          xlab("Spearman partial correlation") + ylab("Number of genes") + 
          ggtitle(paste("Histogram of Spearman partial correlation of all genes with", lncrna, "\nAdjusted for tissue")) +
          theme(legend.position="none", axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=20), title=element_text(size=19)))
    dev.off()

    this_expression_df <- data.frame(lncrna_expression_vectors[[lncrna]])

    ## get the subjects as their own column
    this_expression_df$subject <- rownames(this_expression_df)
    rownames(this_expression_df) <- NULL
    this_expression_df <- this_expression_df[this_expression_df$subject %in% expression_samples,]
    colnames(this_expression_df)[1] <- "lncrna_expression"
    ## now add tissue categorization to this
    this_expression_df$tissue_category <- sample_info_df$inferno_class[match(this_expression_df$subject, gsub("-", ".", sample_info_df$SAMPID))]

    make_graphic(paste0(outdir, '/plots/', lncrna_out, '/', lncrna_out, '_gtex_expression_plot'))
    print(ggplot(this_expression_df, aes(x=tissue_category, y=lncrna_expression, fill=tissue_category)) +
          geom_boxplot() + theme_bw() + xlab("Tissue category") + ylab(paste(lncrna, "RPKM expression")) +
          ggtitle(paste(lncrna, "expression across GTEx samples")) +
          cat_col_scale + 
          theme(legend.position="none", axis.text.x = element_text(size=25, angle=45, hjust=1),
                axis.text.y = element_text(size=25), title=element_text(size=20),
                plot.title = element_text(hjust = 0.5)))
    dev.off()    
}

cat(nrow(all_lncrna_targets), "total genes met the", cor_thresh, "correlation threshold, representing", length(unique(all_lncrna_targets$target_gene)), "unique genes and", length(unique(all_lncrna_targets$lncRNA)), "unique lncRNAs across", length(unique(unlist(strsplit(as.character(all_lncrna_targets$tissues), ";")))), "GTEx tissues and", length(unique(unlist(strsplit(as.character(all_lncrna_targets$classes), ";")))), "tissue classes\n")
cat(nrow(all_lncrna_targets), "total genes met the", cor_thresh, "correlation threshold, representing", length(unique(all_lncrna_targets$target_gene)), "unique genes and", length(unique(all_lncrna_targets$lncRNA)), "unique lncRNAs across", length(unique(unlist(strsplit(as.character(all_lncrna_targets$tissues), ";")))), "GTEx tissues and", length(unique(unlist(strsplit(as.character(all_lncrna_targets$classes), ";")))), "tissue classes\n", file=summary_file, append=T)

## now write the full target file
write.table(all_lncrna_targets, paste0(outdir, '/tables/all_lncRNA_targets_', cor_thresh, '_correlation_threshold.txt'), quote=F, sep="\t", row.names=F, col.names=T)
## also write all the target genes together, for pathway analysis
write.table(sort(unique(all_lncrna_targets$target_gene)), paste0(outdir, '/tables/all_lncRNA_genes_', cor_thresh, '_correlation_threshold.txt'), quote=F, sep="\t", row.names=F, col.names=F)
## save an Rdata file of the list of correlation DFs
save(lncrna_correlation_dfs, file=paste0(outdir, '/full_correlation_tables/full_lncRNA_correlation_df_list.Rdata'))

## now write out tissue-specific lists of target genes
dir.create(paste0(outdir, '/tables/tissue_specific_gene_lists/'), F, T)
for(gtex_tissue in unique(unlist(strsplit(as.character(all_lncrna_targets$tissues), ";")))) {
    write.table(sort(unique(all_lncrna_targets$target_gene[grep(gtex_tissue, all_lncrna_targets$tissues)])), paste0(outdir, '/tables/tissue_specific_gene_lists/', gtex_tissue, '_lncRNA_genes_', cor_thresh, '_correlation_threshold.txt'), quote=F, sep="\t", row.names=F, col.names=F)
}

dir.create(paste0(outdir, '/tables/class_specific_gene_lists/'), F, T)
for(gtex_class in unique(unlist(strsplit(as.character(all_lncrna_targets$classes), ";")))) {
    write.table(sort(unique(all_lncrna_targets$target_gene[grep(gtex_class, all_lncrna_targets$classes)])), paste0(outdir, '/tables/class_specific_gene_lists/', gtex_class, '_lncRNA_genes_', cor_thresh, '_correlation_threshold.txt'), quote=F, sep="\t", row.names=F, col.names=F)
}

## ## to read the target file back in
## all_lncrna_targets <- read.table(paste0(outdir, '/tables/all_lncRNA_targets_', cor_thresh, '_correlation_threshold.txt'), header=T, sep="\t", quote="", as.is=T)
## ## to read the correlations back in
## load(file=paste0(outdir, '/full_correlation_tables/full_lncRNA_correlation_df_list.Rdata'))

## -----------------------------------------------------------------------------
## 6. Summarize correlation patterns across all lncRNAs
## -----------------------------------------------------------------------------
## make this a factor so that we have counts for all lncRNAs
all_lncrna_targets$lncRNA <- factor(all_lncrna_targets$lncRNA, levels=names(lncrna_correlation_dfs))

## count number of targets per lncRNA:
lncrna_target_counts <- ddply(all_lncrna_targets, .(lncRNA), summarize, num_targets=length(unique(target_gene)), .drop=FALSE)

## plot the histogram of number of targets per lncRNA
if(max(lncrna_target_counts$num_targets) > 100) {
    scale_amt <- 50
} else {
    scale_amt <- 10
}

max_y_count <- max(count(lncrna_target_counts$num_targets)$freq)

make_graphic(paste0(outdir, '/plots/', outprefix, '_all_lncRNA_target_num_histogram'))
print(ggplot(lncrna_target_counts, aes(x=num_targets)) + 
      geom_histogram(binwidth=1) +
      scale_x_continuous(breaks=seq(0, round_any(max(lncrna_target_counts$num_targets), accuracy=scale_amt, f=ceiling), by=scale_amt)) +
      scale_y_continuous(breaks=seq(0, max_y_count, by=1)) + 
      ggtitle("Histogram of number of gene targets per lncRNA") + 
      theme_bw() + xlab("Number of unique target gene(s) per lncRNA") +
      ylab("Number of lncRNAs") + 
      theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
            axis.text.y = element_text(size=15), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## count number of lncRNAs per gene target (i.e. how many lncRNAs target a given gene)
target_lncrna_counts <- ddply(all_lncrna_targets, .(target_gene), summarize, num_lncrnas=length(unique(lncRNA)))

## plot the histogram of number of lncRNAs per target gene
max_y_count <- max(count(target_lncrna_counts$num_lncrnas)$freq)
if(max_y_count > 1000) {
    y_scale_amt <- 500
} else {
    y_scale_amt <- 100
}

make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_lncRNA_num_histogram'))
print(ggplot(target_lncrna_counts, aes(x=num_lncrnas)) + 
      geom_histogram(binwidth=1) + 
      scale_x_continuous(breaks=seq(0, max(target_lncrna_counts$num_lncrnas), by=1)) +
      scale_y_continuous(breaks=seq(0, max_y_count, by=y_scale_amt)) + 
      ggtitle("Histogram of number of lncRNAs per gene") + 
      theme_bw() + xlab("Number of unique lncRNAs targeting gene") +
      ylab("Number of target genes") + 
      theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=15),
            axis.text.y = element_text(size=15), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## make a histogram of the pearson correlation values across all targets
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_pearson_histogram'))
print(ggplot(all_lncrna_targets, aes(x=pearson_cor)) + 
      geom_histogram(binwidth=0.01) + 
      ggtitle("Histogram of lncRNA-mRNA Pearson correlations") + 
      theme_bw() + xlab("lncRNA - mRNA expression Pearson correlation") +
      ylab("Number of lncRNA - target gene pairs") + 
      theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=20),
            axis.text.y = element_text(size=20), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## make a histogram of the spearman correlation values across all targets
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_spearman_histogram'))
print(ggplot(all_lncrna_targets, aes(x=spearman_cor)) + 
      geom_histogram(binwidth=0.01) + 
      ggtitle("Histogram of lncRNA-mRNA Spearman correlations") + 
      theme_bw() + xlab("lncRNA - mRNA expression Spearman correlation") +
      ylab("Number of lncRNA - target gene pairs") + 
      theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=20),
            axis.text.y = element_text(size=20), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## make a combined density plot for this
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_pearson_and_spearman_dist'))
print(ggplot(melt(all_lncrna_targets, id.vars=c("lncRNA", "target_gene"), variable.name="cor"),
             aes(x=value, fill=cor)) +
      geom_density(alpha=0.5) + scale_fill_hue(h=c(270, 360)) + 
      ggtitle("Spearman and Pearson correlation distributions") + 
      theme_bw() + xlab("Correlation value") +
      ylab("Density") + 
      guides(fill=guide_legend(title="Type of correlation")) + 
      theme(legend.position="bottom", axis.text.x=element_text(angle=60, hjust=1, size=15),
            axis.text.y = element_text(size=20), legend.text = element_text(size=20), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## plot distributions of correlation values split by lncRNA
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_pearson_dist_per_lncRNA'))
print(ggplot(all_lncrna_targets, aes(x=lncRNA, y=pearson_cor, fill=lncRNA)) +
      geom_boxplot() + scale_fill_hue(h=c(270, 360)) + 
      ggtitle("Pearson correlation distributions per lncRNA") + 
      theme_bw() + xlab("lncRNA") +
      ylab("Pearson correlations of lncRNA - target gene pairs") + 
      theme(legend.position="none", axis.text.x=element_text(angle=60, hjust=1, size=15),
            axis.text.y = element_text(size=20), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_spearman_dist_per_lncRNA'))
print(ggplot(all_lncrna_targets, aes(x=lncRNA, y=spearman_cor, fill=lncRNA)) +
      geom_boxplot() + scale_fill_hue(h=c(270, 360)) + 
      ggtitle("Spearman correlation distributions per lncRNA") + 
      theme_bw() + xlab("lncRNA") +
      ylab("Spearman correlations of lncRNA - target gene pairs") + 
      theme(legend.position="none", axis.text.x=element_text(angle=60, hjust=1, size=15),
            axis.text.y = element_text(size=20), 
            title=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## make a boxplot with combined distributions
melted_lncrna_targets <- melt(all_lncrna_targets, measure.vars=c("pearson_cor", "spearman_cor"), variable="cor_type")
melted_lncrna_targets$lncRNA <- factor(melted_lncrna_targets$lncRNA, levels=sort(unique(as.character(melted_lncrna_targets$lncRNA))), ordered=T)

make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_spearman_and_pearson_dist_per_lncRNA'), width_ratio = 2.0)
print(ggplot(melted_lncrna_targets, aes(x=lncRNA, y=value, fill=cor_type)) +
      geom_boxplot() + 
      scale_fill_manual(values=c("#FF6C91", "#9F8CFF"), labels=c("Pearson", "Spearman"),
                        name="Type of correlation") + 
      ## ggtitle("Spearman and Pearson correlation distributions per lncRNA") + 
      theme_bw() + xlab("lncRNA") +
      ylab("Correlation values") + 
      theme(legend.position="bottom", axis.text.x=element_text(angle=60, hjust=1, size=18),
            axis.text.y = element_text(size=25), 
            title=element_text(size=25),
            legend.text=element_text(size=25),
            plot.title=element_text(hjust=0.5)))
dev.off()

## finally, make a big scatterplot of the correlations
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_correlation_scatterplot'),
             height_ratio = 1.5)
print(ggplot(all_lncrna_targets, aes(x=pearson_cor, y=spearman_cor, fill=lncRNA)) +
      geom_point(aes(color=lncRNA), alpha=0.5, shape=5) + scale_color_hue(h=c(0, 360)) + 
      theme_bw() +
      ## geom_hline(yintercept=cor_thresh, linetype=3) +
      ## geom_vline(xintercept=cor_thresh, linetype=3) +
      ## geom_hline(yintercept=-cor_thresh, linetype=3) +
      ## geom_vline(xintercept=-cor_thresh, linetype=3) +
      xlab("Pearson correlation") + ylab("Spearman correlation") +
      ggtitle("Correlation values of top lncRNA - mRNA pairs") +
      theme(legend.position="bottom", axis.text.x = element_text(size=20),
            axis.text.y = element_text(size=20), title=element_text(size=20),
            plot.title = element_text(hjust = 0.5), legend.text = element_text(size=15)))
dev.off()

## make this scatterplot without splitting by lncRNA
make_graphic(paste0(outdir, '/plots/', outprefix, '_all_target_correlation_scatterplot_across_lncRNAs')## , width_ratio = 2.0, height_ratio = 2.0
             )
print(ggplot(all_lncrna_targets, aes(x=pearson_cor, y=spearman_cor)) +
      stat_binhex(bins=50) + 
#      geom_point(aes(color=lncRNA), alpha=0.5, shape=5) + scale_color_hue(h=c(0, 360)) + 
      theme_bw() +
      ## geom_hline(yintercept=cor_thresh, linetype=3) +
      ## geom_vline(xintercept=cor_thresh, linetype=3) +
      ## geom_hline(yintercept=-cor_thresh, linetype=3) +
      ## geom_vline(xintercept=-cor_thresh, linetype=3) +
      xlab("Pearson correlation") + ylab("Spearman correlation") +
      ggtitle("Correlation values of top lncRNA - mRNA pairs") +
      theme(legend.position="bottom", axis.text.x = element_text(size=20),
            legend.key.width=unit(3,"cm"),
            axis.text.y = element_text(size=20), title=element_text(size=20),
            plot.title = element_text(hjust = 0.5), legend.text = element_text(size=15)))
dev.off()

## make expression distribution plots for all target genes (?)
