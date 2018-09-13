## all_gene_GTEx_correlations.R
## alex amlie-wolf 06/13/18
## analyzing correlation distribution of all genes in the genome

library(ggplot2)
library(scales)
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

save(pearson_cor, spearman_cor, file=paste0(outdir, "/full_correlation_tables/pearson_and_spearman_matrices.Rdata"))

## load(paste0(outdir, "/full_correlation_tables/pearson_and_spearman_matrices.Rdata"))

summary(pearson_cor[upper.tri(pearson_cor)])
summary(spearman_cor[upper.tri(spearman_cor)])

## {
## pearson_cor[lower.tri(pearson_cor)] <- NA
## spearman_cor[lower.tri(spearman_cor)] <- NA
## }

## {
##     melt_start <- proc.time()
##     ## visualize the correlation patterns
##     pearson_correlation_melt <- melt(pearson_cor, na.rm=TRUE, as.is=TRUE)
##     cat("Melting Pearson took", (proc.time() - melt_start)[['elapsed']], 'seconds\n')

##     melt_start <- proc.time()
##     spearman_correlation_melt <- melt(spearman_cor, na.rm=TRUE, as.is=TRUE)
##     cat("Melting Spearman took", (proc.time() - melt_start)[['elapsed']], 'seconds\n')    
## }

{
    combn_start <- proc.time()
    pearson_vec <- c(pearson_cor[upper.tri(pearson_cor)])
    rm(pearson_cor)
    cat("Pearson vectorization took", (proc.time() - combn_start)[['elapsed']], 'seconds\n')

    combn_start <- proc.time()
    spearman_vec <- c(spearman_cor[upper.tri(spearman_cor)])
    rm(spearman_cor)
    cat("Spearman vectorization took", (proc.time() - combn_start)[['elapsed']], 'seconds\n')
}

save(pearson_vec, spearman_vec, file=paste0(outdir, "/full_correlation_tables/pearson_and_spearman_vectors.Rdata"))

## load(paste0(outdir, "/full_correlation_tables/pearson_and_spearman_vectors.Rdata"))

## check the different quantiles
tab_start <- proc.time()

table(pearson_pos=pearson_vec >= 0.5, spearman_pos=spearman_vec >= 0.5, pearson_neg=pearson_vec <= -0.5, spearman_neg=spearman_vec <= -0.5)

cat("Tables took", (proc.time() - tab_start)[['elapsed']], "seconds\n")

sum(rownames(pearson_cor) != rownames(spearman_cor))
sum(colnames(pearson_cor) != colnames(spearman_cor))

## {
##     merge_start <- proc.time()
##     combined_correlation_melt <- merge(pearson_correlation_melt, spearman_correlation_melt, by=c("gene_a", "gene_b"), suffixes=c(".pearson", ".spearman"))
##     cat("Merging took", (proc.time() - merge_start)[['elapsed']], "seconds\n")
## }

summary(spearman_vec)

## get an index for subsampling the vectors for this
{
## num_pairs <- length(pearson_vec)    
## vec_sample <- seq(1, num_pairs, by=2)
## correlation_df <- data.frame(pearson=pearson_vec[vec_sample], spearman=spearman_vec[vec_sample])
correlation_df <- data.frame(pearson=pearson_vec, spearman=spearman_vec)
}

{
make_graphic(paste0(outdir, '/plots/GTEx_v6p_all_correlation_distribution'), height_ratio=2.0, width_ratio=2.0)
print(ggplot(correlation_df, aes(x=pearson, y=spearman)) + 
      stat_binhex(bins=75) +
      scale_fill_gradientn(colors=c("lightgray", "cyan", "cyan1", "cyan2", "cyan3", "cyan4",
                               "darkcyan", "blue", "blue1", "blue2", "blue3", "blue4", "darkblue"),
                          name="Number of gene-gene pairs", labels=comma) +
      theme_bw() +
      geom_hline(yintercept=0.5, linetype=3) +
      geom_vline(xintercept=0.5, linetype=3) +
      geom_hline(yintercept=-0.5, linetype=3) +
      geom_vline(xintercept=-0.5, linetype=3) +    
      scale_x_continuous(limits=c(-1, 1), expand = c(0.05, 0),
                         breaks=seq(-1, 1, by=0.2)) +
      scale_y_continuous(limits=c(-1, 1), expand = c(0.05, 0),
                         breaks=seq(-1, 1, by=0.2)) +
      xlab("Pearson correlation") + ylab("Spearman correlation") +
      ggtitle("Correlation distribution of all transcript pairs") +
      theme(legend.position="bottom", axis.text.x = element_text(size=20),
            legend.key.width=unit(3,"cm"),
            axis.text.y = element_text(size=20), title=element_text(size=20),
            plot.title = element_text(hjust = 0.5), legend.text = element_text(size=15)))
dev.off()

}

{
make_graphic(paste0(outdir, '/plots/GTEx_v6p_all_correlation_distribution_logscale'), height_ratio=2.0, width_ratio=2.0)
print(ggplot(correlation_df, aes(x=pearson, y=spearman)) + 
      stat_binhex(bins=75, aes(fill=log(..count..))) +
      scale_fill_gradientn(colors=c("lightgray", "cyan", "cyan1", "cyan2", "cyan3", "cyan4",
                               "darkcyan", "blue", "blue1", "blue2", "blue3", "blue4", "darkblue"),
                          name="Number of gene-gene pairs", labels=comma) +
      theme_bw() +
      geom_hline(yintercept=0.5, linetype=3) +
      geom_vline(xintercept=0.5, linetype=3) +
      geom_hline(yintercept=-0.5, linetype=3) +
      geom_vline(xintercept=-0.5, linetype=3) +    
      scale_x_continuous(limits=c(-1, 1), expand = c(0.05, 0),
                         breaks=seq(-1, 1, by=0.2)) +
      scale_y_continuous(limits=c(-1, 1), expand = c(0.05, 0),
                         breaks=seq(-1, 1, by=0.2)) +
      xlab("Pearson correlation") + ylab("Spearman correlation") +
      ggtitle("Correlation distribution of all transcript pairs") +
      theme(legend.position="bottom", axis.text.x = element_text(size=20),
            legend.key.width=unit(3,"cm"),
            axis.text.y = element_text(size=20), title=element_text(size=20),
            plot.title = element_text(hjust = 0.5), legend.text = element_text(size=15)))
dev.off()

}
