## analyze_GTEx_v7_metaXcan.R
## alex amlie-wolf 06/13/18
## a script to analyze the results of the metaXcan GTEx v7 analysis

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
## 3. Analyze MetaXcan results

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------
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
## 2. Read in data, set data directories
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==1) {
    datadir <- args[1]
} else {
    stop("Requires MetaXcan output directory")
}

dir.create(paste0(datadir, '/metaxcan_analysis/plots/'), F, T)
dir.create(paste0(datadir, '/metaxcan_analysis/tables/'), F, T)

## read in all the datasets, save in a single data frame
full_metaxcan_results <- data.frame(stringsAsFactors=FALSE)
for(meta_file in list.files(datadir, pattern="*csv")) {
    this_tissue <- gsub("results-gtex_v7_|_imputed_eur.csv", "", meta_file)
    cat("Reading data from", this_tissue, "\n")
    this_data <- read.table(paste0(datadir, '/', meta_file), header=T, sep=",", quote="", as.is=T)

    full_metaxcan_results <- rbind(full_metaxcan_results,
                                   cbind(tissue=this_tissue, this_data))
}

## correct all the p-values
full_metaxcan_results$bonf_p <- p.adjust(full_metaxcan_results$pvalue, method="bonferroni")

write.table(full_metaxcan_results,
            paste0(datadir, '/metaxcan_analysis/tables/full_metaxcan_results.txt'),
            quote=F, sep="\t", row.names=F)

cat(sum(full_metaxcan_results$bonf_p <= 0.05), "significant tissue-gene pairs found with Bonferroni cutoff of 0.05 across all tissues and genes\n")

write.table(full_metaxcan_results[full_metaxcan_results$bonf_p <= 0.05,],
            paste0(datadir, '/metaxcan_analysis/tables/full_metaxcan_results.0.05_bonf_signif.txt'),
            quote=F, sep="\t", row.names=F)

## -----------------------------------------------------------------------------
## 3. Analyze MetaXcan results
## -----------------------------------------------------------------------------
for(var in colnames(full_metaxcan_results)[!(colnames(full_metaxcan_results) %in% c("tissue", "gene", "gene_name"))]) {
    make_graphic(paste0(datadir, '/metaxcan_analysis/plots/', var, '_distribution_across_tissues'),
                 width_ratio = 1.5)
    print(ggplot(full_metaxcan_results,
                 aes_string(x="tissue", y=var, color="tissue")) +
          geom_boxplot() +
          ggtitle(paste("Distributions of", var, "across tissues")) +
          theme_bw() + xlab("Tissue") + ylab(var) +
          theme(legend.position="none", axis.text.x=element_text(angle=60, hjust=1, size=15),
                axis.text.y = element_text(size=15), strip.text=element_text(size=25),
                legend.text = element_text(size=15), legend.title=element_text(size=25),
                axis.title=element_text(size=25),
                plot.title=element_text(hjust=0.5, size=20)))
    dev.off()

    if(sum(full_metaxcan_results$bonf_p <= 0.05) > 0) {
        make_graphic(paste0(datadir, '/metaxcan_analysis/plots/', var, '_distribution_across_tissues.signif'),
                     width_ratio = 1.5)
        print(ggplot(full_metaxcan_results[full_metaxcan_results$bonf_p <= 0.05,],
                     aes_string(x="tissue", y=var, color="tissue")) +
              geom_boxplot() +
              ggtitle(paste("Distributions of", var, "across tissues for significant results")) +
              theme_bw() + xlab("Tissue") + ylab(var) +
              theme(legend.position="none", axis.text.x=element_text(angle=60, hjust=1, size=15),
                    axis.text.y = element_text(size=15), strip.text=element_text(size=25),
                    legend.text = element_text(size=15), legend.title=element_text(size=25),
                    axis.title=element_text(size=25),
                    plot.title=element_text(hjust=0.5, size=20)))
        dev.off()
    }    
}
    
