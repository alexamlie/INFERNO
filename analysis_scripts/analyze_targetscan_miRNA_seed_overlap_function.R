## analyze_targetscan_miRNA_seed_overlap_function.R
## alex amlie-wolf 06/22/17
## analyzes targetscan-predicted miRNA seed site overlaps

analyze_targetscan_seed_overlaps <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    targetscan_file <- paste0(datadir, '/targetscan_miRNA_overlap/', prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh, "_targetscan_overlap.txt")
    ## read in this data
    targetscan_overlap_df <- read.table(targetscan_file, header=T, sep="\t", quote="", as.is=T)
    
    if(nrow(targetscan_overlap_df)==0) {
        ## still write out the table so that other stuff doesn't break
        targetscan_overlap_df <- cbind(targetscan_overlap_df, tag_no_rsid=numeric(0))
        ## write this table out
        targetscan_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                             dist_thresh, '_targetscan_miRNA_seed_overlap.txt')
        write.table(targetscan_overlap_df, targetscan_outf, quote=F, sep="\t", col.names=T, row.names=F)
        cat("No TargetScan miRNA seed overlaps found in this dataset!\n")
        return("No TargetScan miRNA seed overlaps found in this dataset!")
    }

    ## add annotation for the tag regions without tag rsIDs
    targetscan_overlap_df$tag_no_rsid <- gsub(":rs.*", "", targetscan_overlap_df$tag_name)

    ## write this table out
    targetscan_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                         dist_thresh, '_targetscan_miRNA_seed_overlap.txt')
    write.table(targetscan_overlap_df, targetscan_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## make a barplot of the different targetScan conditions
    make_graphic(paste0(outdir, 'plots/', prefix, '_targetscan_condition_counts_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(count(targetscan_overlap_df, vars="targetscan_condition"),
                 aes(x=targetscan_condition, y=freq, fill=targetscan_condition)) +
          xlab("TargetScan Condition") + ylab("Number of Overlaps") +
          theme_bw() + geom_bar(stat="identity", position="dodge") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Counts of overlaps from each TargetScan condition", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## distributions of context score percentiles by condition
    make_graphic(paste0(outdir, 'plots/', prefix, '_targetscan_scores_by_condition_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(targetscan_overlap_df, aes(x=targetscan_condition, y=context_score_percentile, fill=targetscan_condition)) +
          xlab("TargetScan Condition") + ylab("Context Score Percentile") +
          theme_bw() + geom_violin() +
          geom_jitter(aes(fill=targetscan_condition), colour='black', shape=23, alpha=0.5, height=0) +          
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of context score percentiles from each TargetScan condition", r2_thresh, dist_thresh, out_subtitle))
    dev.off()
    
}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/targetscan_miRNA_overlap/")

