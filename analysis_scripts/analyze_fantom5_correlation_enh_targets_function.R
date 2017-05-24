## analyze_fantom5_correlation_enh_targets.R
## alex amlie-wolf 01/11/16
## analysis of the enhancer targets defined by enh:promoter correlation

library(ggplot2)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------

## 0. Table of Contents
## 1. Function Definitions
## 2. Read in data
## 3. Analyze enhancer overlap targets
## 4. Analyze closest enhancer targets

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------
## the function library should be in the same place
source("enhancer_analysis_functions.R")

## -----------------------------------------------------------------------------
## 2. Read in data
## -----------------------------------------------------------------------------
outdir <- "/home/alexaml/data/enhancer_snp_pipeline/results/correlation_enh_targets/"

outprefix <- "IGAP_top_hits"
dir.create(paste0(outdir, 'plots/', outprefix), F, T)
dir.create(paste0(outdir, 'tables/', outprefix), F, T)

param_list <- get_params(outprefix)
inprefix <- param_list$inprefix
datadir <- paste0(param_list$datadir, '/correlation_enh_targets/')
out_subtitle <- param_list$out_subtitle
## get the cross-referencing description file
description_xref <- param_list$description_xref

## select the values to analyze
r2_thresh <- 0.7
dist_thresh <- 500000

## read in the overlap data file
enh_overlap_target_file <- paste0(datadir, inprefix, "_", format(r2_thresh, nsmall=1), "_ld_cutoff_snps_within_", format(dist_thresh, scientific=F), "_fantom5_overlap_target_genes.txt")
enh_overlap_target_df <- read.table(enh_overlap_target_file, header=T, sep="\t", quote="", as.is=T)
## fix the R and FDR columns
enh_overlap_target_df$r <- as.numeric(gsub("R:", "", enh_overlap_target_df$r, fixed=T))
enh_overlap_target_df$fdr <- as.numeric(gsub("FDR:", "", enh_overlap_target_df$fdr, fixed=T))

## read in the closest enhancer data file
closest_enh_target_file <- paste0(datadir, inprefix, "_", format(r2_thresh, nsmall=1), "_ld_cutoff_snps_within_", format(dist_thresh, scientific=F), "_closest_enhancer_target_genes.txt")
closest_enh_target_df <- read.table(closest_enh_target_file, header=T, sep="\t", quote="", as.is=T)
## fix the R and FDR columns
closest_enh_target_df$r <- as.numeric(gsub("R:", "", closest_enh_target_df$r, fixed=T))
closest_enh_target_df$fdr <- as.numeric(gsub("FDR:", "", closest_enh_target_df$fdr, fixed=T))

## get the linkage information
ld_stats_file <- paste0(param_list$datadir, '/ld_expansion/', inprefix, "_", format(r2_thresh, nsmall=1),
                         "_ld_cutoff_snps_within_", format(dist_thresh, scientific=F), ".txt")
ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)
## merge the linkage info with the others
enh_overlap_merge_df <- merge(ld_stats_df, enh_overlap_target_df[,-3],
                              by.x=c("chr", "pos", "rsID", "tag_rsID"),
                              by.y=c("chr", "start", "rsID", "tag_rsID"))
closest_enh_merge_df <- merge(ld_stats_df, closest_enh_target_df[,-3],
                              by.x=c("chr", "pos", "rsID", "tag_rsID"),
                              by.y=c("chr", "start", "rsID", "tag_rsID"))

## if the description exists, read it in and augment the input datasets
if(!is.na(description_xref)) {
    description_tab <- read.table(description_xref, header=F, sep="\t", quote="", as.is=T,
                                  col.names=c("snp_rsID", "description"))
    rownames(description_tab) <- description_tab$snp_rsID
    ## paste the description with the tag rsID for each dataset
    enh_overlap_target_df$tag_rsID <- paste(enh_overlap_target_df$tag_rsID,
                                  paste0("(", description_tab[enh_overlap_target_df$tag_rsID,"description"],")"),
                                  sep="\n")

    closest_enh_target_df$tag_rsID <- paste(closest_enh_target_df$tag_rsID,
                                  paste0("(", description_tab[closest_enh_target_df$tag_rsID,"description"],")"),
                                  sep="\n")

    ## also do the merged tables, but use a space separator
    enh_overlap_merge_df$tag_rsID <- paste(enh_overlap_merge_df$tag_rsID,
                                      paste0("(", description_tab[enh_overlap_merge_df$tag_rsID,"description"],")"),
                                      sep=" ")    
    closest_enh_merge_df$tag_rsID <- paste(closest_enh_merge_df$tag_rsID,
                                      paste0("(", description_tab[closest_enh_merge_df$tag_rsID,"description"],")"),
                                      sep=" ")    
}

## write out the merged tables:
full_overlap_outf <- paste0(outdir, 'tables/', outprefix, '/', outprefix, "_", format(r2_thresh, nsmall=1),
                         "_ld_cutoff_snps_within_", format(dist_thresh, scientific=F),
                        "_fantom5_overlap_targets_with_linkage.txt")
write.table(enh_overlap_merge_df, full_overlap_outf, quote=F, sep="\t", col.names=T, row.names=F)

full_closest_enh_outf <- paste0(outdir, 'tables/', outprefix, '/', outprefix, "_", format(r2_thresh, nsmall=1),
                         "_ld_cutoff_snps_within_", format(dist_thresh, scientific=F),
                        "_closest_enh_targets_with_linkage.txt")
write.table(closest_enh_merge_df, full_closest_enh_outf, quote=F, sep="\t", col.names=T, row.names=F)
rm(enh_overlap_merge_df, closest_enh_merge_df, ld_stats_df)


## -----------------------------------------------------------------------------
## 3. Analyze enhancer overlap targets
## -----------------------------------------------------------------------------
## plot the distribution of correlation values
make_graphic(paste0(outdir, 'plots/', outprefix, '/', outprefix, '_enh_overlap_target_r_vals_',
                    format(r2_thresh, nsmall=1), "_ld_", format(dist_thresh, scientific=F),
                    "_dist"), type='png')
ggplot(enh_overlap_target_df[enh_overlap_target_df$num_assocs > 0,],
       aes(x=tag_rsID, y=r, colour=r)) +
    scale_fill_hue(h=c(180, 270)) + ylab("R value") + xlab("Tag SNP rsID / info") +
    theme_bw() + geom_boxplot(fill=NA) +
    geom_point(position=position_jitter(width=0.2), alpha=0.4) +
    theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(paste0("Correlation values of target promoters for enhancers overlapping SNPs\nR2 >= ", format(r2_thresh, nsmall=1), ", Distance <= ", format(dist_thresh, scientific=F), " bp\n", out_subtitle)) 
dev.off()

## plot the number of unique LD SNPs per tag SNP
num_ld_snps <- ddply(enh_overlap_target_df, .(tag_rsID), summarize,
                     num_enh_overlap_snps = length(unique(rsID)))
## add 0-counts for non-overlapping ones:
snp_descriptions <- paste(description_tab$snp_rsID,
                          paste0("(", description_tab$description,")"), sep="\n")
if(sum(!(snp_descriptions %in% enh_overlap_target_df$tag_rsID)) > 0) {
    num_ld_snps <- rbind(num_ld_snps,
                         data.frame(tag_rsID = snp_descriptions[which(!(snp_descriptions %in% enh_overlap_target_df$tag_rsID))],
                                    num_enh_overlap_snps = 0, stringsAsFactors = F))
}

## re-order to match the input order
num_ld_snps <- num_ld_snps[match(description_tab$snp_rsID, gsub("\n.*", "", num_ld_snps$tag_rsID)),]

## look at number of linked SNPs per tag SNP with enhancer overlap + target
make_graphic(paste0(outdir, 'plots/', outprefix, '/', outprefix, '_num_linked_snps_with_corr_per_tagsnp_',
                    format(r2_thresh, nsmall=1), "_ld_", format(dist_thresh, scientific=F),
                    "_dist"))  
ggplot(num_ld_snps, aes(x=tag_rsID, y=num_enh_overlap_snps, fill=tag_rsID)) +
    xlab("Tag SNP rsID") + ylab("Number of SNPs overlapping any FANTOM5 enhancers") +
    scale_fill_hue(h=c(180, 270)) +
    scale_y_discrete(breaks=seq(0, max(num_ld_snps$num_enh_overlap_snps)*1.1, by=5)) + 
    theme_bw() + geom_bar(stat="identity", position="stack") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) + 
    ggtitle(paste0("Numbers of linked SNPs with enhancer overlap and correlation-based target per tag SNP\nR2 >= ", format(r2_thresh, nsmall=1), ", Distance <= ", format(dist_thresh, scientific=F), " bp\n", out_subtitle))
dev.off()


## -----------------------------------------------------------------------------
## 4. Analyze closest enhancer targets
## -----------------------------------------------------------------------------
## plot the distribution of correlation values
make_graphic(paste0(outdir, 'plots/', outprefix, '/', outprefix, '_closest_enh_target_r_vals_',
                    format(r2_thresh, nsmall=1), "_ld_", format(dist_thresh, scientific=F),
                    "_dist"), type='png')
ggplot(closest_enh_target_df[closest_enh_target_df$num_assocs > 0,],
       aes(x=tag_rsID, y=r, fill=tag_rsID)) +
    scale_fill_hue(h=c(180, 270)) + ylab("R value") + xlab("Tag SNP rsID / info") +
    theme_bw() + geom_boxplot() +
    geom_point(position=position_jitter(width=0.2), alpha=0.4) +
    theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(paste0("Correlation values of target promoters for closest enhancers to SNPs\nR2 >= ", format(r2_thresh, nsmall=1), ", Distance <= ", format(dist_thresh, scientific=F), " bp\n", out_subtitle)) 
dev.off()

## compare distance to enhancer with correlation
make_graphic(paste0(outdir, 'plots/', outprefix, '/', outprefix, '_closest_enh_r_vs_dist_to_enh_',
                    format(r2_thresh, nsmall=1), "_ld_", format(dist_thresh, scientific=F),
                    "_dist"), type='png', height_ratio = 1.5)
ggplot(closest_enh_target_df[closest_enh_target_df$num_assocs>0,],
       aes(x=abs(dist_to_enh), y=r)) +
    facet_grid(tag_rsID ~ .) +
    ylab("R value") + xlab("Absolute value of distance to closest enhancer") + 
    theme_bw() + geom_point(aes(colour=abs(dist_to_enh))) +
    theme(legend.position="none") + # , axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(paste0("Correlation value vs distance of closest enhancers to SNPs\nR2 >= ", format(r2_thresh, nsmall=1), ", Distance <= ", format(dist_thresh, scientific=F), " bp\n", out_subtitle)) 
dev.off()
